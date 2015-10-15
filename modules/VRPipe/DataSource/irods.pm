
=head1 NAME

VRPipe::DataSource::irods - get pipeline inputs from an iRods data store

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013-2015 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see L<http://www.gnu.org/licenses/>.

=cut

use VRPipe::Base;

class VRPipe::DataSource::irods with VRPipe::DataSourceFilterRole {
    use Digest::MD5 qw(md5_hex);
    use File::Spec::Functions;
    use Path::Class;
    use VRPipe::Persistent::InMemory;
    use VRPipe::Schema;
    use File::Which qw(which);
    use JSON::XS;
    use IPC::Run;
    # https://rt.cpan.org/Public/Bug/Display.html?id=102824
    # we tie STDERR in vrpipe-handler which will call this module, but since
    # we have no FILENO method in VRPipe::Persistent::InMemory, IPC::Run will
    # break. The following is a hack to work around:
    {
        no warnings 'redefine';
        *IPC::Run::_debug_fd = sub { };
    }
    
    our %ignore_keys = map { $_ => 1 } qw(study_id study_title study_accession_number sample_common_name ebi_sub_acc reference ebi_sub_md5 ebi_run_acc ebi_sub_date sample_created_date taxon_id lane);
    
    has '_irods_files_and_metadata_cache' => (
        is        => 'rw',
        isa       => 'HashRef',
        clearer   => '_clear_cache',
        predicate => '_cached',
    );
    
    has '_analysis_to_col' => (
        is      => 'rw',
        isa     => 'HashRef',
        default => sub { {} }
    );
    
    method description {
        return "Use an iRods data store to extract information from";
    }
    
    method source_description {
        return "The zone to query";
    }
    
    method method_description (Str $method) {
        if ($method eq 'all') {
            return "An element will comprise one of the files returned by imeta qu -d given the arguments you supply for the 'file_query' option (which can be the options specified directly, or the absolute path to a file containing multiple sets of imeta options, 1 set per line). The file will have all the relevant irods metadata associated with it, and a local path based on the 'local_root_dir' option (none supplied means the dataelements will have paths like irod:/abs/path/to/file.txt, and the step(s) in your pipeline will have to be able to handle such paths). To avoid spamming the irods server, the update_interval option allows you to specify the minimum number of minutes between each check for changes to files. If update_interval is not supplied it defaults to 5 seconds when testing, and 2 hours in production; the maximum is 4 hours.";
        }
        if ($method eq 'all_with_warehouse_metadata') {
            return "In addition to doing everything the all method does, it adds extra metadata found in the warehouse database to the files with the keys public_name, sample_supplier_name, sample_control, sample_cohort, taxon_id, sample_created_date and study_title (if defined). Optionally provide a comma-separated list of required keys to the required_metadata option to ignore files lacking that metadata. If any analysis has been done to a file, the associated files are stored under irods_analysis_files. If any qc files have been stored with the file, these are associated with a qc_file relationship; the require_qc_files option will skip files that do not have all the qc files specified in desired_qc_files, which is a comma-separated list of main-file-basename(excluding extension) suffixes, defaulting to _F0xB00.stats,.genotype.json,.verify_bam_id.json. The vrtrack_group option determines which group the studies your data are in are placed under - use it for grouping together studies you will analyse the same way later. The graph_filter option is a string of the form 'namespace#label#propery#value'; multiple filters can be separated by commas. The filter will look for an exact match to a property of a node that the file's node is descended from, eg. specify VRTrack#Sample#qc_failed#0 to only have files related to samples that have not been qc failed. (This method is Sanger-specific and also requires the environment variables WAREHOUSE_DATABASE, WAREHOUSE_HOST, WAREHOUSE_PORT and WAREHOUSE_USER.)";
        }
        if ($method eq 'group_by_metadata_with_warehouse_metadata') {
            return "Extension to the all_with_warehouse_metadata methods that will group files from the source according to their metadata keys. Requires the metadata_keys option which is a '|' separated list of metadata keys by which dataelements will be grouped. e.g. metadata_keys => 'sample|platform|library' will groups all files with the same sample, platform and library into one dataelement. If you use a graph_filter and the filter_after_grouping option is set (the default), grouping based on metadata will be performed first and then the filter applied with it only being necessary for one file in the group to pass the filter. If the filter_after_grouping option is not set, only files which pass the filter will be included and grouped based on their metadata.";
        }
        return '';
    }
    
    method _open_source {
        return $self->source; # the source can't really be opened
    }
    
    method _has_changed {
        my $old      = $self->_changed_marker;
        my $checksum = $self->_irods_file_metadata_checksum;
        $self->_changed_marker($checksum);
        $old || return 1;
        return ($checksum ne $old) ? 1 : 0;
    }
    
    method _update_changed_marker {
        $self->_changed_marker($self->_irods_file_metadata_checksum);
    }
    
    method _irods_file_metadata_checksum {
        # we need a quick way of looking at all the relevant metadata for all
        # our desired files, but imeta/iquest doesn't give us this, so we
        # use baton instead.
        # to avoid abusing the irods server too much we only actually calculate
        # the checksum anew every x mins; otherwise we just return the current
        # checksum
        my $im      = VRPipe::Persistent::InMemory->new;
        my $options = $self->options;
        my $update_interval;
        if ($options->{update_interval}) {
            # in production the user would supply minutes; in testing they
            # would supply seconds
            if ($VRPipe::Persistent::InMemory::deployment eq 'production') {
                $update_interval = $options->{update_interval} * 60;
            }
            else {
                $update_interval = $options->{update_interval};
            }
            
            # this used to be really slow so we had defaults or user settings
            # of like a day, so we'll set it down to at most 4 hours
            if ($update_interval > 14400) {
                $update_interval = 14400;
            }
        }
        else {
            # default to 2 hours
            $update_interval = 7200;
        }
        my $lock_key = 'irods_datasource.' . $self->_datasource_id;
        
        my $current_checksum = $self->_changed_marker;
        if ($im->is_alive($lock_key)) {
            warn "irods datasource will return without checking since it it currently being updated in another process\n" if $self->debug;
            return $current_checksum;
        }
        
        my $locked                      = $im->lock($lock_key, unlock_after => $update_interval);
        my $graph_filter                = $options->{graph_filter};
        my $filter_after_grouping       = $options->{filter_after_grouping};
        my $required_metadata           = $options->{required_metadata} || '';
        my $vrtrack_group               = $options->{vrtrack_group} || '';
        my $require_qc_files            = $options->{require_qc_files} || 0;
        my $desired_qc_files            = defined $options->{desired_qc_files} ? $options->{desired_qc_files} : '_F0xB00.stats,.genotype.json,.verify_bam_id.json';
        my $add_metadata_from_warehouse = $self->method =~ /with_warehouse_metadata$/ ? 1 : 0;
        my $local_root_dir              = $options->{local_root_dir} || '';
        
        if ($current_checksum && length($current_checksum) == 32) {
            unless ($locked) {
                warn "irods datasource will return without checking since it's been less than $update_interval seconds since the last check\n" if $self->debug;
                return $current_checksum;
            }
        }
        # else we always get the latest checksum if we have no valid checksum
        $self->debug_log("irods datasource will do a full check since it's been more than $update_interval seconds since the last check\n");
        
        #*** while neo4j suffers from too many of these running at once,
        # we'll only run 4 irods datasource updates at once, globally
        my $ids_queue_key = 'irods_datasource_queue';
        $im->enqueue($ids_queue_key, $lock_key);
        my $can_run        = 0;
        my $warned_waiting = 0;
        my $max            = 4;
        while (!$can_run) {
            my $running = 0;
            foreach my $queued ($im->queue($ids_queue_key)) {
                if ($im->is_alive($queued)) {
                    $running++;
                }
                elsif ($queued eq $lock_key) {
                    $can_run = 1 if $running < $max;
                    last;
                }
                elsif ($running < $max) {
                    next;
                }
                else {
                    last;
                }
            }
            
            if (!$can_run) {
                sleep(5);
            }
            
            if (!$can_run && !$warned_waiting) {
                $self->debug_log(" (will wait until there are not $max other irods datasources updating)\n");
                $warned_waiting = 1;
            }
        }
        
        # recheck we're not updating elsewhere, since we could have a race
        # condition following the prior is_alive() call. assert_life() avoids
        # the race condition but we didn't want to use it until after the above
        # queue was resolved.
        unless ($im->assert_life($lock_key)) {
            warn "irods datasource will return without checking since it it currently being updated in another process\n" if $self->debug;
            return $current_checksum;
        }
        
        # get the current files and their metadata and stringify it all
        my $t     = time();
        my $files = $self->_get_irods_files_and_metadata($self->_open_source(), $options->{file_query}, $local_root_dir, $add_metadata_from_warehouse, $required_metadata, $vrtrack_group, $require_qc_files, $desired_qc_files, $graph_filter, $filter_after_grouping);
        my $e     = time() - $t;
        $self->debug_log("irods _get_irods_files_and_metadata call took $e seconds\n");
        $self->_irods_files_and_metadata_cache($files);
        
        $im->dequeue($ids_queue_key, [$lock_key]);
        
        my $data = '';
        while (my ($path, $meta) = each %$files) {
            $data .= $path;
            
            foreach my $key (sort keys %$meta) {
                next if $key eq 'vrpipe_irods_order';
                next if $key eq 'irods_datasource_changed';
                next if $key eq 'irods_datasource_graph_filter';
                my $val = $meta->{$key};
                next unless defined $val;
                # limit to printable ascii so md5_hex subroutine will work
                $val =~ tr/\x20-\x7f//cd;
                $data .= "|$key:$val|";
            }
            
            if ($graph_filter) {
                my $pass = $meta->{irods_datasource_graph_filter};
                $pass ||= 0;
                $data .= "|graph_filter:$pass|";
            }
        }
        
        my $digest = md5_hex $data;
        return $digest;
    }
    
    method _get_irods_files_and_metadata (Str $zone!, Str $raw_query!, Str $local_root_dir!, Bool $add_metadata_from_warehouse!, Str $required_metadata!, Str $vrtrack_group!, Bool $require_qc_files?, Str $desired_qc_files?, Maybe[Str] $graph_filter?, Bool $filter_after_grouping?) {
        return $self->_irods_files_and_metadata_cache if $self->_cached;
        
        my $debug = $self->debug;
        
        my (undef, $gfs, $vrpipe_graph_schema, $graph) = $self->_parse_filters(undef, $graph_filter);
        
        my ($baton_metaquery, $baton_list);
        unless (($baton_metaquery = which('baton-metaquery')) && ($baton_list = which('baton-list'))) {
            $self->throw('baton-metaquery and baton-list must be in your PATH for the irods datasource to work (see http://wtsi-npg.github.io/baton/)');
        }
        
        my @required_keys;
        if ($required_metadata) {
            @required_keys = split(',', $required_metadata);
        }
        my @desired_qc_file_suffixes = split(/,/, $desired_qc_files) if $desired_qc_files;
        my $desired_qc_file_suffixes = join('|', @desired_qc_file_suffixes) if @desired_qc_file_suffixes;
        my $desired_qc_file_regex = qr/(\S+(?:$desired_qc_file_suffixes))/ if $desired_qc_file_suffixes;
        
        my ($sample_sth, $public_name, $donor_id, $supplier_name, $control, $taxon_id, $created, $gender, $internal_id);
        my ($study_sth, $study_title, $study_ac, $study_id_sth, $warehouse_study_id, $sample_donor_sth, $sample_donor_sample, $sample_donor_donor, %taxons, %genders);
        my $vrtrack;
        if ($add_metadata_from_warehouse && $ENV{WAREHOUSE_DATABASE} && $ENV{WAREHOUSE_HOST} && $ENV{WAREHOUSE_PORT} && $ENV{WAREHOUSE_USER}) {
            my $dbh = DBI->connect(
                "DBI:mysql:host=$ENV{WAREHOUSE_HOST}:port=$ENV{WAREHOUSE_PORT};database=$ENV{WAREHOUSE_DATABASE}",
                $ENV{WAREHOUSE_USER}, undef, { 'RaiseError' => 1, 'PrintError' => 0, mysql_enable_utf8 => 1 }
            );
            
            # name in warehouse corresponds to 'sample' metadata from irods, and
            # is one of the few indexed columns, so queries against it should
            # hopefully be quick. For situations where the irods metadata is
            # missing or out-of-date we want to get values from warehouse and
            # store them under the metadata key that irods uses:
            # warehouse_table,warehouse_column => irods_key
            # current_samples,supplier_name => sample_supplier_name
            # current_samples,donor_id => sample_cohort
            # current_samples,control => sample_control
            # current_samples,taxon_id => n/a (store under taxon_id)
            # current_samples,created => n/a (store under sample_created_date)
            # current_studies,name => study_title (linked to irods study_id => internal_id)
            
            # notes for /archive/GAPI/exp/infinium idat files:
            # lane name = basename of filename
            # library name = sample.beadchip.beadchip_section
            # library ssid = last 4 digits of beadchip . last 3 digits of sample_id
            # lane acc = analysis_uuid (the last one listed by imeta)
            # lane storage_path = {user dir}/analysis_uuid
            
            # notes for /archive/GAPI/gen/infinium gtc files:
            # lane name = basename of filename, or {beadchip}_{beadchip_section}
            # library name = infinium_sample
            # library ssid = last 4 digits of beadchip . last 3 digits of sample_id
            # lane acc = analysis_uuid (the last one listed by imeta)
            # lane storage_path = {user dir}/{analysis_uuid}_{basename(user dir)}
            
            # notes for /seq bam files:
            # lane name = basename of filename
            # individual name = sample_supplier_name or sample
            # individual acc = sample_accession_number
            
            my $sql = q[select public_name, donor_id, supplier_name, control, taxon_id, created, gender, internal_id from current_samples where name = ?];
            $sample_sth = $dbh->prepare($sql);
            $sample_sth->execute;
            $sample_sth->bind_columns(\($public_name, $donor_id, $supplier_name, $control, $taxon_id, $created, $gender, $internal_id));
            
            $sql       = q[select name, accession_number from current_studies where internal_id = ?];
            $study_sth = $dbh->prepare($sql);
            $study_sth->execute;
            $study_sth->bind_columns(\($study_title, $study_ac));
            
            # we don't trust irods metadata for the study id, so will query it
            $sql          = q[select study_internal_id from current_study_samples where sample_internal_id = ?];
            $study_id_sth = $dbh->prepare($sql);
            $study_id_sth->execute;
            $study_id_sth->bind_col(1, \$warehouse_study_id);
            
            # because donor_id column isn't indexed in current_samples, it's
            # quicker for us to just get the samples of all donors with this
            # query, then use the previous query to get the study details
            $sql              = q[select internal_id, donor_id from current_samples where donor_id IS NOT NULL];
            $sample_donor_sth = $dbh->prepare($sql);
            $sample_donor_sth->execute;
            $sample_donor_sth->bind_columns(\($sample_donor_sample, $sample_donor_donor));
            
            $vrtrack = VRPipe::Schema->create('VRTrack');
            $graph ||= $vrtrack->graph;
            $vrtrack_group ||= 'all_studies';
        }
        
        my %files;
        my @queries;
        if ($raw_query !~ /\s/ && -e $raw_query) {
            # it's a file containing multiple queries
            my $file = VRPipe::File->create(path => $raw_query);
            my $fh = $file->openr;
            while (<$fh>) {
                chomp;
                next if /^#/; # skip comments
                push(@queries, $_) if $_;
            }
            $file->close;
        }
        else {
            # it's a single directly-specified query
            @queries = ($raw_query);
        }
        $self->debug_log("have " . scalar(@queries) . " queries to deal with\n");
        
        # subs for working with baton-* which will be fed json
        my $json = JSON::XS->new->allow_nonref(1);
        
        my %harnesses;
        my $run_baton;
        $run_baton = sub {
            my ($cmd, $json_input, $retries) = @_;
            $retries ||= 0;
            
            # we use IPC::Run to hold open a pipe to and from baton for each
            # different baton command...
            my $cmd_str = join(' ', @$cmd);
            my ($h, $in, $out, $err);
            unless (defined $harnesses{$cmd_str}) {
                my $i = '';
                my $o = '';
                my $e = '';
                $in  = \$i;
                $out = \$o;
                $err = \$e;
                $h   = IPC::Run::harness($cmd, $in, $out, $err);
                $harnesses{$cmd_str} = [$h, $in, $out, $err];
            }
            else {
                ($h, $in, $out, $err) = @{ $harnesses{$cmd_str} };
            }
            
            # ... then send the command some json and wait for the json
            # response
            ${$err} = '';
            ${$in} .= $json_input . "\n";
            $h->pump until ${$out} =~ m{[\r\n]$}msx;
            
            if (${$err}) {
                my @errors = grep { $_ ne 'The client/server socket connection has been renewed' } split(/\n/, ${$err});
                if (@errors) {
                    $self->debug_log(" irods datasource baton error while processing [$json_input | $cmd_str]:\n  " . join("\n  ", @errors) . "\n");
                }
            }
            
            # now we decode and return the json response
            my $result;
            if (${$out}) {
                chomp(${$out});
                $result = $json->decode(${$out});
                ${$out} = '';
            }
            
            if (!$result || (ref($result) eq 'HASH' && exists $result->{error})) {
                my $err = $result ? qq/: "$result->{error}->{message}"/ : '';
                my $msg = qq/baton returned nothing from [$json_input | $cmd_str]$err/;
                if ($retries < 5) {
                    $self->debug_log(" irods datasource error: $msg; will retry\n");
                    $retries++;
                    sleep($retries);
                    
                    # we'll recurse to retry the command with a fresh connection
                    $h->finish;
                    delete $harnesses{$cmd_str};
                    return &$run_baton($cmd, $json_input, $retries);
                }
                else {
                    my $fatal = "irods datasource fatal error: $msg; giving up\n";
                    $self->debug_log($fatal);
                    die $fatal;
                }
            }
            
            return $result;
        };
        
        my $run_baton_metadata = sub {
            my ($query, $mode) = @_; # $mode is files|dirs
            my $type = $mode eq 'dirs' ? '--coll' : '--obj';
            my @cmd = ($baton_metaquery, qw(--unbuffered --zone), $zone, qw(--checksum --avu), $type);
            
            # convert query string to baton-compatible json
            # study = "Exome - Neo-antigen discovery in mouse lung cancer" and target = 1 and type = cram
            my @avus;
            foreach my $avu (split(/\s+and\s+/i, $query)) {
                my ($attribute, $o, $value) = $avu =~ /^(\S+)\s+(\S+)\s+(.+)/;
                foreach my $s (\$attribute, \$value, \$o) {
                    if (${$s} && ${$s} =~ /^['"].*['"]$/) {
                        ${$s} =~ s/^['"]//;
                        ${$s} =~ s/['"]$//;
                    }
                }
                push(@avus, { attribute => $attribute, value => $value, o => $o });
            }
            my $json_input = $json->encode({ avus => \@avus });
            
            my $baton_output = &$run_baton(\@cmd, $json_input);
            
            my @files;
            foreach my $file_hash (@$baton_output) {
                my $collection = $file_hash->{collection};
                my $basename   = $file_hash->{data_object};
                
                # convert all the avus in to a flat metadata hash
                my $meta = {};
                foreach my $avu (@{ $file_hash->{avus} }) {
                    my $attribute = $avu->{attribute};
                    next if $mode eq 'files' && $attribute =~ /^dcterms:/;
                    next if $attribute =~ /_history$/;
                    if (exists $meta->{$attribute}) {
                        my $val = $meta->{$attribute};
                        if (ref($val)) {
                            push(@$val, $avu->{value});
                        }
                        else {
                            $val = [$val, $avu->{value}];
                        }
                        $meta->{$attribute} = $val;
                    }
                    else {
                        $meta->{$attribute} = $avu->{value};
                    }
                }
                
                push(@files, { dir => $collection, $basename ? (basename => $basename) : (), metadata => $meta });
            }
            
            return \@files;
        };
        
        my $run_baton_list;
        $run_baton_list = sub {
            my ($collection, $mode, $recursive_files_only) = @_; # $mode is contents|metadata
            my $type = $mode eq 'contents' ? '--contents' : '--avu';
            my @cmd = ($baton_list, '--unbuffered', $type);
            
            # convert collection string to baton-compatible json
            # '{collection: "/seq/fluidigm/34/b9/41/1662049042/11/d7/28"}'
            my $json_input = $json->encode({ collection => $collection });
            
            my $baton_output = &$run_baton(\@cmd, $json_input);
            
            if ($mode eq 'contents') {
                my @files;
                foreach my $content (@{ $baton_output->{contents} }) {
                    my $dir      = $content->{collection};
                    my $basename = $content->{data_object};
                    next if ($basename && $basename =~ /[~\$ ]/); # just ignore annoying files we don't need anyway
                    
                    if ($recursive_files_only && !$basename) {
                        push(@files, &$run_baton_list($dir, $mode, 1));
                        next;
                    }
                    
                    push(@files, $dir . ($basename ? "/$basename" : ''));
                }
                return @files;
            }
            else {
                my $meta = {};
                foreach my $avu (@{ $baton_output->{avus} }) {
                    my $attribute = $avu->{attribute};
                    if (exists $meta->{$attribute}) {
                        my $val = $meta->{$attribute};
                        if (ref($val)) {
                            push(@$val, $avu->{value});
                        }
                        else {
                            $val = [$val, $avu->{value}];
                        }
                        $meta->{$attribute} = $val;
                    }
                    else {
                        $meta->{$attribute} = $avu->{value};
                    }
                }
                return $meta;
            }
        };
        
        # use Time::HiRes qw(gettimeofday tv_interval);
        # my @times;
        
        my $tt = time();
        my (%analysis_to_cols, %col_dates, %analysis_files, %analyses, %collections, %ils_cache, %done_donors, %donor_to_samples, %done_studies, $rel_args, $rel_args_with_props, %done_samples);
        my $ran_sample_donor_sth = 0;
        my $order                = 1;
        foreach my $query (@queries) {
            my $t                         = time();
            my $irods_files_with_metadata = &$run_baton_metadata($query, 'files');
            my $e                         = time() - $t;
            $self->debug_log(" baton query for [$query] took $e seconds to run, returning " . scalar(@$irods_files_with_metadata) . " files with their metadata\n");
            
            $t = time();
            my $f_count        = 0;
            my $f_skip_qc      = 0;
            my $f_skip_meta    = 0;
            my $f_count_actual = 0;
            FILE: foreach my $file_hash (@$irods_files_with_metadata) {
                my $collection = $file_hash->{dir};
                my $basename   = $file_hash->{basename};
                my $path       = "$collection/$basename";
                my $meta       = $file_hash->{metadata};
                $f_count++;
                print STDERR '. ' if $debug;
                $meta->{vrpipe_irods_order} = $order++;
                $meta->{expected_md5} = $file_hash->{checksum} if $file_hash->{checksum};
                
                # a Sanger-specific thing is that we can have qc-related
                # files stored in the same folder or qc sub-folder as our
                # main file, so we'll associate these files with the main
                # file if they exist, and skip if some don't exist and
                # we're in require mode
                my %qc_files;
                if ($vrtrack && $desired_qc_file_regex) {
                    my $prefix = $basename;
                    $prefix =~ s/\.[^\.]+$//;
                    my $desired_regex = qr/$prefix(?:$desired_qc_file_suffixes)/;
                    
                    my $ils_output = [];
                    if (exists $ils_cache{$collection}) {
                        $ils_output = $ils_cache{$collection};
                    }
                    else {
                        my @base_ils_output = &$run_baton_list($collection, 'contents');
                        my @qc_ils_output;
                        foreach (reverse(@base_ils_output)) {
                            if ($_ eq "$collection/qc") {
                                @qc_ils_output = &$run_baton_list("$collection/qc", 'contents');
                                last;
                            }
                        }
                        
                        foreach (@base_ils_output, @qc_ils_output) {
                            if (/^$desired_qc_file_regex$/) {
                                push(@$ils_output, $_);
                            }
                        }
                        
                        $ils_cache{$collection} = $ils_output;
                    }
                    
                    my @qc_files;
                    foreach my $path (@$ils_output) {
                        my $basename = file($path)->basename;
                        if ($basename =~ /^$desired_regex$/) {
                            push(@qc_files, $path);
                        }
                    }
                    
                    if ($require_qc_files && @qc_files != @desired_qc_file_suffixes) {
                        $f_skip_qc++;
                        $self->debug_log("  $path didn't have all required qc files\n");
                        next FILE;
                    }
                    elsif (@qc_files) {
                        foreach my $qc_file_path (@qc_files) {
                            $qc_files{$qc_file_path} = 1;
                        }
                    }
                }
                
                # grab extra info from the warehouse database if requested
                if ($sample_sth) {
                    undef $public_name;
                    undef $donor_id;
                    undef $supplier_name;
                    undef $control;
                    undef $taxon_id;
                    undef $created;
                    undef $gender;
                    undef $internal_id;
                    my $sanger_sample_id = $meta->{sample};
                    
                    if ($sanger_sample_id) {
                        $sample_sth->execute($sanger_sample_id);
                        $sample_sth->fetch;
                        if ($public_name) {
                            $meta->{public_name} = "$public_name";
                        }
                        if ($donor_id) {
                            $meta->{sample_cohort} = "$donor_id";
                        }
                        if ($supplier_name) {
                            $meta->{sample_supplier_name} = "$supplier_name";
                        }
                        if (defined $control) {
                            $meta->{sample_control} = "$control";
                        }
                        if (defined $taxon_id) {
                            $meta->{taxon_id} = "$taxon_id";
                        }
                        if (defined $created) {
                            $meta->{sample_created_date} = "$created";
                        }
                        if (defined $gender) {
                            $meta->{sample_gender} = "\u$gender";
                        }
                    }
                    
                    # override irods study_id with warehouse study_id, and
                    # set study_title
                    if (defined $internal_id) {
                        undef $warehouse_study_id;
                        $study_id_sth->execute($internal_id);
                        my (@study_ids, @study_titles, @study_acs);
                        while ($study_id_sth->fetch) {
                            push(@study_ids, "$warehouse_study_id");
                            
                            undef $study_title;
                            undef $study_ac;
                            $study_sth->execute($warehouse_study_id);
                            $study_sth->fetch;
                            if ($study_title) {
                                push(@study_titles, "$study_title");
                                $study_ac ||= ''; # ac isn't always set in warehouse
                                push(@study_acs, "$study_ac");
                            }
                        }
                        
                        # because we most likely want to populate VRTrack
                        # based on this metadata, and VRTrack can't cope
                        # with multiple study_ids per sample, we'll check
                        # the source to see if we wanted just a single study
                        # and make note of that one
                        my $preferred_i;
                        if (@study_ids > 1) {
                            $meta->{study_id} = \@study_ids;
                            
                            if ($query =~ /study_id\s+=\s+(\d+)/) {
                                my $desired = $1;
                                foreach my $i (0 .. $#study_ids) {
                                    my $id = $study_ids[$i];
                                    if ($id == $desired) {
                                        $meta->{study_id_preferred} = $desired;
                                        $preferred_i = $i;
                                        last;
                                    }
                                }
                            }
                        }
                        else {
                            $meta->{study_id} = $study_ids[0];
                        }
                        
                        if (@study_titles > 1) {
                            $meta->{study_title}            = \@study_titles;
                            $meta->{study_accession_number} = \@study_acs;
                            
                            if ($query =~ /study\s+=\s+["']([^"']+)["']/) {
                                my $desired = $1;
                                foreach my $i (0 .. $#study_titles) {
                                    my $title = $study_titles[$i];
                                    if ($title eq $desired) {
                                        $meta->{study_id_preferred} = $study_ids[$i];
                                        last;
                                    }
                                }
                            }
                        }
                        else {
                            $meta->{study_title} = $study_titles[0];
                            $meta->{study_accession_number} = $study_acs[0] if $study_acs[0];
                        }
                    }
                    
                    # another Sanger-specific thing is if we had an
                    # analysis_uuid, see if there are any collections with the
                    # same analysis_uuid and store some of the associated files
                    if (exists $meta->{analysis_uuid}) {
                        # get the most recent collection associated with an
                        # analysis done for this file
                        my @uuids = ref($meta->{analysis_uuid}) ? @{ $meta->{analysis_uuid} } : ($meta->{analysis_uuid});
                        my $analysis_to_col_hash = $self->_analysis_to_col;
                        my %collections;
                        foreach my $uuid (@uuids) {
                            if (exists $analysis_to_cols{$uuid}) {
                                foreach my $this_col (@{ $analysis_to_cols{$uuid} }) {
                                    $collections{$this_col} = 1;
                                }
                            }
                            else {
                                my $irods_dirs_with_metadata = &$run_baton_metadata("analysis_uuid = $uuid", 'dirs');
                                
                                my @these_cols;
                                foreach my $dir_hash (@$irods_dirs_with_metadata) {
                                    my $this_col = $dir_hash->{dir};
                                    $collections{$this_col} = 1;
                                    push(@these_cols, $this_col);
                                    
                                    unless (exists $col_dates{$this_col}) {
                                        my $date = $dir_hash->{metadata}->{'dcterms:created'} || '2013-01-01T12:00:00';
                                        $col_dates{$this_col} = $date;
                                        $analysis_to_col_hash->{uuids}->{$uuid}->{$this_col} = $vrtrack->date_to_epoch($date);
                                    }
                                }
                                
                                $analysis_to_cols{$uuid} = \@these_cols;
                            }
                        }
                        
                        my ($analysis_collection) = sort { $col_dates{$b} cmp $col_dates{$a} } keys %collections;
                        
                        if ($analysis_collection) {
                            my @files;
                            if (exists $analysis_files{$analysis_collection}) {
                                @files = @{ $analysis_files{$analysis_collection} };
                            }
                            else {
                                my @file_contents = &$run_baton_list($analysis_collection, 'contents', 1);
                                foreach my $file (@file_contents) {
                                    push(@files, $file);
                                    $analysis_to_col_hash->{cols}->{$analysis_collection}->{$file} = 1;
                                }
                                
                                $analysis_files{$analysis_collection} = [@files];
                            }
                            
                            if (@files) {
                                $meta->{irods_analysis_files} = @files > 1 ? \@files : $files[0];
                            }
                        }
                    }
                }
                
                foreach my $key (@required_keys) {
                    unless (defined $meta->{$key}) {
                        $f_skip_meta++;
                        next FILE;
                    }
                }
                
                # represent lab tracking metadata in the graph database
                my $graph_file;
                if ($vrtrack) {
                    # we'll only do the heavy-duty graph access if the file has
                    # changed (or is new)
                    #my $gtod = [gettimeofday];
                    my $graph_file = $vrtrack->get_file($path, 'irods:');
                    #$times[0] += tv_interval($gtod);
                    
                    if (!$graph_file || (my $changes = $self->_file_changed($path, $meta, $local_root_dir, 1, $vrtrack))) {
                        my $already_in_graph = 0;
                        if ($graph_file) {
                            my @related = $graph_file->related();
                            if (@related > 2) {
                                $already_in_graph = 1;
                            }
                        }
                        else {
                            $graph_file = $vrtrack->add_file($path, 'irods:');
                        }
                        
                        my $from_scratch = 1;
                        if ($changes && ref($changes)) { # $changes is 0 for no change, 1 for brand new, and an array ref of changes if changed
                            $meta->{vrpipe_meta_changed} = $changes;
                            
                            if ($already_in_graph) {
                                # try and quickly fix up basic changes without
                                # having to merge every single node and rel
                                my $handled   = 0;
                                my $unhandled = '';
                                my $punted    = 0;
                                foreach my $change (@$changes) {
                                    my ($key, $old, $new, $diff) = @$change;
                                    if ($key eq 'public_name') {
                                        my $sample = $vrtrack->get('Sample', { name => $meta->{sample} });
                                        if ($sample) {
                                            $sample->public_name($new);
                                            $handled++;
                                            next;
                                        }
                                    }
                                    elsif ($key eq 'sample_cohort' || $key eq 'sample_donor_id') {
                                        # this is a little bit involved; punting
                                        # for now
                                        $punted++;
                                    }
                                    
                                    $unhandled .= "$key: $diff; ";
                                }
                                
                                if ($handled == @$changes) {
                                    $from_scratch = 0;
                                }
                                elsif ($punted != @$changes) {
                                    $self->debug_log("  {temp debug; tell Sendu about this!} couldn't quickly handle all changes for $path [$unhandled]\n");
                                }
                            }
                        }
                        
                        if ($from_scratch) {
                            print STDERR '! ' if $debug;
                            
                            # relate any qc files to this file
                            foreach my $qc_file_path (keys %qc_files) {
                                my $qc_file = $vrtrack->add_file($qc_file_path, 'irods:');
                                $graph_file->relate_to($qc_file, 'qc_file');
                            }
                            
                            # if _get_irods_files_and_metadata() stored collection info, relate
                            # those now as well
                            my $analysis_to_col_hash = $self->_analysis_to_col;
                            if (exists $meta->{analysis_uuid}) {
                                foreach my $uuid (ref($meta->{analysis_uuid}) ? @{ $meta->{analysis_uuid} } : ($meta->{analysis_uuid})) {
                                    my $cols_hash = $analysis_to_col_hash->{uuids}->{$uuid} || next;
                                    my $vrtrack_analysis = $vrtrack->add('Analysis', { uuid => $uuid });
                                    $graph_file->relate_to($vrtrack_analysis, 'analysed');
                                    while (my ($this_col, $date) = each %{$cols_hash}) {
                                        my $vrtrack_col = $vrtrack->add('Collection', { path => $this_col, date => $date }, incoming => { type => 'directory', node => $vrtrack_analysis });
                                        
                                        my $file_hash = $analysis_to_col_hash->{cols}->{$this_col} || next;
                                        foreach my $file (keys %{$file_hash}) {
                                            my $vrtrack_file = $vrtrack->add_file($file, 'irods:');
                                            $vrtrack_col->relate_to($vrtrack_file, 'contains');
                                        }
                                        delete $analysis_to_col_hash->{cols}->{$this_col};
                                    }
                                    delete $analysis_to_col_hash->{uuids}->{$uuid};
                                }
                            }
                            
                            # general
                            my $group = $vrtrack->add('Group', { name => $vrtrack_group });
                            my (@studies, $preferred_study, %study_relate_to_props);
                            my @study_ids = ref $meta->{study_id} ? @{ $meta->{study_id} } : ($meta->{study_id});
                            my @study_titles = ($meta->{study_title} && ref $meta->{study_title}) ? @{ $meta->{study_title} } : ($meta->{study_title} || $meta->{study});
                            my @study_acs = ($meta->{study_accession_number} && ref $meta->{study_accession_number}) ? @{ $meta->{study_accession_number} } : ($meta->{study_accession_number});
                            my $preferred_study_id = delete $meta->{study_id_preferred};
                            foreach my $i (0 .. $#study_ids) {
                                my $study_title = $study_titles[$i] || '[no study title]';
                                my $done_studies_key = 'id:' . $study_ids[$i] . ':name:' . $study_title . ':ac:' . ($study_acs[$i] ? $study_acs[$i] : '') . ':group:' . $vrtrack_group;
                                if (exists $done_studies{$done_studies_key}) {
                                    push(@studies, $done_studies{$done_studies_key});
                                }
                                else {
                                    my $study_node = $vrtrack->add('Study', { id => $study_ids[$i], name => $study_title, $study_acs[$i] ? (accession => $study_acs[$i]) : () }, incoming => { type => 'has', node => $group });
                                    push(@studies, $study_node);
                                    $done_studies{$done_studies_key} = $study_node;
                                }
                                
                                if (defined $preferred_study_id && $preferred_study_id == $study_ids[$i]) {
                                    $preferred_study = $studies[-1];
                                    %study_relate_to_props = (properties => { preferred => 1 });
                                    
                                    # in VRPipe::File metadata we'll just store the
                                    # preferred details; this is inline with old
                                    # behaviour, is probably what the user expects,
                                    # and avoids problems if you're grouping on
                                    # study_id and then your sample gets added to a
                                    # new study
                                    $meta->{study_id}    = $study_ids[$i];
                                    $meta->{study_title} = $study_titles[$i];
                                    $meta->{study}       = $study_titles[$i] if defined $meta->{study};
                                    if ($study_acs[$i] && defined $meta->{study_accession_number}) {
                                        $meta->{study_accession_number} = $study_acs[$i];
                                    }
                                    else {
                                        delete $meta->{study_accession_number};
                                    }
                                }
                            }
                            
                            my $donor;
                            my $donor_uuid = $meta->{sample_cohort} if exists $meta->{sample_cohort};
                            if (defined $donor_uuid) {
                                unless (exists $done_donors{$donor_uuid}) {
                                    $donor = $vrtrack->add('Donor', { id => $donor_uuid });
                                    
                                    # since donor_id isn't indexed in
                                    # current_samples, we do a 1-time pull of all
                                    # sample rows in the table (that have donor_id
                                    # set)
                                    unless ($ran_sample_donor_sth) {
                                        undef $sample_donor_sample;
                                        undef $sample_donor_donor;
                                        $sample_donor_sth->execute();
                                        while ($sample_donor_sth->fetch) {
                                            $donor_to_samples{"$sample_donor_donor"}->{"$sample_donor_sample"} = 1;
                                        }
                                        $ran_sample_donor_sth = 1;
                                    }
                                    
                                    # find out all the studies that all the donor's
                                    # samples belong to
                                    my @study_details;
                                    my %good_studies;
                                    foreach my $internal_id (keys %{ $donor_to_samples{$donor_uuid} }) {
                                        undef $warehouse_study_id;
                                        $study_id_sth->execute($internal_id);
                                        while ($study_id_sth->fetch) {
                                            undef $study_title;
                                            undef $study_ac;
                                            $study_sth->execute($warehouse_study_id);
                                            $study_sth->fetch;
                                            push(@study_details, ["$warehouse_study_id", "$study_title", $study_ac ? ("$study_ac") : ()]);
                                            $good_studies{"$warehouse_study_id"} = 1;
                                        }
                                    }
                                    
                                    # a donor and all it's samples have been removed
                                    # from a study; now that we know the correct set
                                    # of studies we can remove bad links and add
                                    # missing ones
                                    my (@divorce_args, %already_linked_studies);
                                    foreach my $study ($donor->closest('VRTrack', 'Study', direction => 'incoming', all => 1, depth => 1)) {
                                        my $sid = $study->id;
                                        if (exists $good_studies{$sid}) {
                                            $already_linked_studies{$sid} = 1;
                                        }
                                        else {
                                            push(@divorce_args, [$study, $donor, 'member']);
                                        }
                                    }
                                    $graph->mass_divorce(\@divorce_args) if @divorce_args;
                                    
                                    foreach my $sd (@study_details) {
                                        next if exists $already_linked_studies{ $sd->[0] };
                                        
                                        my $study_title = $sd->[1] || '[no study title]';
                                        my $done_studies_key = 'id:' . $sd->[0] . ':name:' . $study_title . ':ac:' . ($sd->[2] ? $sd->[2] : '') . ':group:' . $vrtrack_group;
                                        my $study;
                                        if (exists $done_studies{$done_studies_key}) {
                                            $study = $done_studies{$done_studies_key};
                                        }
                                        else {
                                            $study = $vrtrack->add('Study', { id => $sd->[0], name => $study_title, $sd->[2] ? (accession => $sd->[2]) : () }, incoming => { type => 'has', node => $group });
                                            $done_studies{$done_studies_key} = $study;
                                        }
                                        
                                        my $relate_args = { from => $study, to => $donor, type => 'member' };
                                        if ($preferred_study && $study->node_id == $preferred_study->node_id) {
                                            $relate_args->{properties} = $study_relate_to_props{properties};
                                            push(@$rel_args_with_props, $relate_args);
                                        }
                                        else {
                                            push(@$rel_args, $relate_args);
                                        }
                                    }
                                    
                                    $done_donors{$donor_uuid} = $donor;
                                }
                                else {
                                    $donor = $done_donors{$donor_uuid};
                                }
                            }
                            
                            my $sample;
                            if (exists $done_samples{ $meta->{sample} }) {
                                $sample = $done_samples{ $meta->{sample} };
                            }
                            else {
                                my $sample_created_date;
                                if (defined $meta->{sample_created_date}) {
                                    # convert '2013-05-10 06:45:32' to epoch seconds
                                    $sample_created_date = $vrtrack->date_to_epoch($meta->{sample_created_date});
                                }
                                $sample = $vrtrack->add('Sample', { name => $meta->{sample}, public_name => $meta->{public_name}, id => $meta->{sample_id}, supplier_name => $meta->{sample_supplier_name}, accession => $meta->{sample_accession_number}, created_date => $sample_created_date, consent => $meta->{sample_consent}, control => $meta->{sample_control} }); #*** this is typically the slowest line in this whole graph update?!
                                
                                # update what studies this sample is linked to by
                                # removing invalid ones and adding missing ones
                                my %good_studies;
                                foreach my $study (@studies) {
                                    $good_studies{ $study->id } = 1;
                                }
                                
                                my (@divorce_args, %already_linked_studies);
                                foreach my $study ($sample->closest('VRTrack', 'Study', direction => 'incoming', all => 1, depth => 1)) {
                                    my $sid = $study->id;
                                    if (exists $good_studies{$sid}) {
                                        $already_linked_studies{$sid} = 1;
                                    }
                                    else {
                                        push(@divorce_args, [$study, $sample, 'member']);
                                    }
                                }
                                $graph->mass_divorce(\@divorce_args) if @divorce_args;
                                
                                foreach my $study (@studies) {
                                    next if exists $already_linked_studies{ $study->id };
                                    
                                    my $relate_args = { from => $study, to => $sample, type => 'member' };
                                    if ($preferred_study && $study->node_id == $preferred_study->node_id) {
                                        $relate_args->{properties} = $study_relate_to_props{properties};
                                        push(@$rel_args_with_props, $relate_args);
                                    }
                                    else {
                                        push(@$rel_args, $relate_args);
                                    }
                                }
                                
                                my $sex = 'U'; # unknown
                                if (defined $meta->{sample_gender}) {
                                    if ($meta->{sample_gender} =~ /^f/i) {
                                        $sex = 'F';
                                    }
                                    elsif ($meta->{sample_gender} =~ /^m/i) {
                                        $sex = 'M';
                                    }
                                    delete $meta->{sample_gender}; # this new thing we don't want on our VRPipe::File metadata
                                }
                                my $gender_md5 = md5_hex('sequencescape.' . $sex);
                                my $gender;
                                if (defined $genders{$gender_md5}) {
                                    $gender = $genders{$gender_md5};
                                }
                                else {
                                    $gender = $vrtrack->add('Gender', { source_gender_md5 => $gender_md5, source => 'sequencescape', gender => $sex });
                                    $genders{$gender_md5} = $gender;
                                }
                                push(@$rel_args, { from => $sample, to => $gender, type => 'gender' });
                                
                                if ($donor) {
                                    #*** we don't have a way of tracking changes to
                                    # relationships, so won't have a record if we swap
                                    # to a different donor here
                                    $donor->relate_to($sample, 'sample', selfish => 1);
                                }
                                
                                if (defined $meta->{taxon_id}) {
                                    my $taxon;
                                    my $taxon_id = $meta->{taxon_id};
                                    if (defined $taxons{$taxon_id}) {
                                        $taxon = $taxons{$taxon_id};
                                    }
                                    else {
                                        $taxon = $vrtrack->get('Taxon', { id => $meta->{taxon_id} });
                                        unless ($taxon) {
                                            $taxon = $vrtrack->add('Taxon', { id => $meta->{taxon_id}, common_name => $meta->{sample_common_name} });
                                        }
                                        $taxons{$taxon_id} = $taxon;
                                    }
                                    $taxon->relate_to($sample, 'member', selfish => 1);
                                }
                                
                                $done_samples{ $meta->{sample} } = $sample;
                            }
                            
                            $graph_file->add_properties({ manual_qc => $meta->{manual_qc}, target => $meta->{target}, md5 => $meta->{md5} });
                            my $file_connected = 0;
                            my $unique         = file($path)->basename;
                            $unique =~ s/\.gz$//;
                            $unique =~ s/\.[^\.]+$//;
                            
                            if ($local_root_dir) {
                                my $local_file = $vrtrack->add_file(file($local_root_dir, $path)->stringify);
                                $graph_file->relate_to($local_file, 'local_file');
                            }
                            
                            # bams
                            my $library;
                            if (defined $meta->{library_id}) {
                                $library = $vrtrack->add('Library', { id => $meta->{library_id}, name => $meta->{library}, tag => $meta->{tag} });
                                $sample->relate_to($library, 'prepared', selfish => 1);
                            }
                            
                            my $lane;
                            if (defined $meta->{lane}) {
                                $lane = $vrtrack->add('Lane', { unique => $unique, lane => $meta->{lane}, run => $meta->{id_run}, total_reads => $meta->{total_reads}, is_paired_read => $meta->{is_paired_read} });
                                
                                $lane->relate_to($graph_file, 'aligned', selfish => 1);
                                $file_connected = 1;
                                
                                if ($library) {
                                    $library->relate_to($lane, 'sequenced', selfish => 1);
                                }
                            }
                            
                            if ($meta->{alignment} && defined $meta->{reference}) {
                                my $aln = $vrtrack->add('Alignment', { reference => $meta->{reference} });
                                $aln->relate_to($graph_file, 'reference', selfish => 1);
                            }
                            
                            if (defined $meta->{ebi_sub_acc} && $meta->{ebi_run_acc}) {
                                my $sub = $vrtrack->add('EBI_Submission', { acc => $meta->{ebi_sub_acc}, defined $meta->{ebi_sub_date} ? (sub_date => $vrtrack->date_to_epoch($meta->{ebi_sub_date})) : () });
                                
                                my $run = $vrtrack->add('EBI_Run', { acc => $meta->{ebi_run_acc}, md5 => $meta->{ebi_sub_md5} }, incoming => { type => 'submitted', node => $graph_file }, outgoing => { type => 'submitted', node => $sub });
                            }
                            
                            # infinium idats and gtc
                            if (defined $meta->{beadchip}) {
                                my $beadchip = $vrtrack->add('Beadchip', { id => $meta->{beadchip}, design => $meta->{beadchip_design} });
                                
                                # clear any old study to beadchip rels before
                                # adding the current ones
                                foreach my $study ($beadchip->closest('VRTrack', 'Study', direction => 'incoming', all => 1, depth => 1)) {
                                    $study->divorce_from($beadchip, 'has');
                                }
                                
                                foreach my $study (@studies) {
                                    $study->relate_to($beadchip, 'has', ($preferred_study && $study->node_id == $preferred_study->node_id) ? (%study_relate_to_props) : ());
                                }
                                
                                if (defined $meta->{beadchip_section}) {
                                    my $section = $vrtrack->add('Section', { unique => $unique, section => $meta->{beadchip_section} }, incoming => { type => 'placed', node => $sample });
                                    $beadchip->relate_to($section, 'section', selfish => 1);
                                    $section->relate_to($graph_file, 'processed', selfish => 1);
                                    $file_connected = 1;
                                }
                            }
                            
                            my $isample;
                            if (defined $meta->{infinium_sample}) {
                                $isample = $vrtrack->add('Infinium_Sample', { id => $meta->{infinium_sample} });
                                $sample->relate_to($isample, 'infinium', selfish => 1);
                            }
                            
                            if (defined $meta->{infinium_plate}) {
                                my $plate = $vrtrack->add('Infinium_Plate', { id => $meta->{infinium_plate} });
                                
                                # clear any old study to plate rels before adding
                                # the current ones
                                foreach my $study ($plate->closest('VRTrack', 'Study', direction => 'incoming', all => 1, depth => 1)) {
                                    $study->divorce_from($plate, 'has');
                                }
                                
                                foreach my $study (@studies) {
                                    $study->relate_to($plate, 'has', ($preferred_study && $study->node_id == $preferred_study->node_id) ? (%study_relate_to_props) : ());
                                }
                                
                                if (defined $meta->{infinium_well}) {
                                    my $well = $vrtrack->add('Well', { unique => $meta->{infinium_plate} . '.' . $meta->{infinium_well}, well => $meta->{infinium_well} }, incoming => { type => 'placed', node => $isample || $sample });
                                    $plate->relate_to($well, 'well', selfish => 1);
                                }
                            }
                            
                            unless ($file_connected) {
                                $sample->relate_to($graph_file, 'processed');
                            }
                            
                            foreach my $these_rel_args ($rel_args, $rel_args_with_props) {
                                if (defined $these_rel_args && @$these_rel_args >= 500) {
                                    $graph->create_mass_relationships($these_rel_args);
                                    $these_rel_args = [];
                                }
                            }
                        }
                    }
                    
                    # apply graph filter
                    if ($gfs) {
                        my $pass = $self->_file_filter($graph_file, $filter_after_grouping, undef, $gfs, $vrpipe_graph_schema, $graph);
                        $meta->{irods_datasource_graph_filter} = $pass;
                    }
                }
                else {
                    my $changes = $self->_file_changed($path, $meta, $local_root_dir, 1);
                    if ($changes && ref($changes)) {
                        $meta->{vrpipe_meta_changed} = $changes;
                    }
                }
                
                $files{$path} = $meta;
                $f_count_actual++;
            }
            
            foreach my $these_rel_args ($rel_args, $rel_args_with_props) {
                if (defined $these_rel_args && @$these_rel_args) {
                    $graph->create_mass_relationships($these_rel_args);
                    $these_rel_args = [];
                }
            }
            
            $e = time() - $t;
            warn "\n" if $debug;
            my $skip_msg = '';
            $skip_msg .= "; $f_skip_meta files skipped due to not having required metadata" if $f_skip_meta;
            $skip_msg .= "; $f_skip_qc files skipped due to not having required qc files"   if $f_skip_qc;
            my $total_files = $f_count_actual < $f_count ? "$f_count_actual/$f_count passing files" : "$f_count files";
            $self->debug_log(" updating the graph db for the $total_files from [$query] took $e seconds$skip_msg\n");
        }
        
        my $total_time = time() - $tt;
        $self->debug_log("all the queries took $total_time seconds to handle, will now update the mysql db with the new set of dataelements\n");
        
        foreach my $ref (values %harnesses) {
            my $h = $ref->[0];
            $h->finish;
        }
        
        return \%files;
    }
    
    method _file_changed (Str $path, HashRef $new_metadata, Str $local_root_dir, Bool $warehouse_mode, Object $vrtrack?) {
        my ($file_abs_path, $protocol);
        if ($local_root_dir) {
            my $sub_path = $path;
            $sub_path =~ s/^\///;
            $file_abs_path = file($local_root_dir, $sub_path)->stringify;
        }
        else {
            $file_abs_path = $path;
            $protocol      = 'irods:';
        }
        
        my ($vrfile) = VRPipe::File->search({ path => $file_abs_path, $protocol ? (protocol => $protocol) : () });
        return 1 unless $vrfile;
        $vrfile = $vrfile->original;
        
        # detect metadata changes
        my $current_metadata = $vrfile->metadata;
        my @changed;
        if ($current_metadata && keys %$current_metadata) {
            while (my ($key, $val) = each %$current_metadata) {
                if ($warehouse_mode) {
                    next if exists $ignore_keys{$key};
                    next if $key =~ /_history$/;
                }
                
                next unless defined $val;
                next unless defined $new_metadata->{$key};
                if (my $diff = $self->_vals_different($val, $new_metadata->{$key})) {
                    push(@changed, [$key, $val, $new_metadata->{$key}, $diff]);
                }
            }
        }
        
        if (!@changed && $vrtrack && defined $new_metadata->{sample} && defined $new_metadata->{public_name}) {
            # double check the most important sample metadata, since $vrfile's
            # metadata could have been updated to latest irods data without
            # the graph sample node being updated to match
            my $sample = $vrtrack->get('Sample', { name => $new_metadata->{sample} });
            if ($sample) {
                my $graph_pn = $sample->public_name;
                my $new_pn   = $new_metadata->{public_name};
                if ($graph_pn ne $new_pn) {
                    push(@changed, ['public_name', $graph_pn, $new_pn, "$graph_pn => $new_pn"]);
                }
            }
            else {
                return 1;
            }
        }
        
        return @changed ? \@changed : 0;
    }
    
    method all (Defined :$handle!, Str :$file_query!, Str|Dir :$local_root_dir?, Str :$update_interval?) {
        my %args;
        $args{handle}          = $handle;
        $args{file_query}      = $file_query;
        $args{local_root_dir}  = $local_root_dir if $local_root_dir;
        $args{update_interval} = $update_interval if defined($update_interval);
        my @element_args;
        my $did = $self->_datasource_id;
        foreach my $result ($self->_all_files(%args)) {
            my $protocol   = $result->{protocol};
            my $irods_path = $result->{irods_path};
            push(@element_args, { datasource => $did, result => { paths => $result->{paths}, $protocol ? (protocol => $protocol) : (), $irods_path ? (irods_path => $irods_path) : () } });
        }
        $self->_create_elements(\@element_args);
    }
    
    method all_with_warehouse_metadata (Defined :$handle!, Str :$file_query!, Str|Dir :$local_root_dir?, Str :$update_interval?, Str :$required_metadata?, Str :$vrtrack_group?, Str :$graph_filter?, Bool :$require_qc_files = 0, Str :$desired_qc_files = '_F0xB00.stats,.genotype.json,.verify_bam_id.json') {
        my %args;
        $args{handle}           = $handle;
        $args{file_query}       = $file_query;
        $args{local_root_dir}   = $local_root_dir if $local_root_dir;
        $args{update_interval}  = $update_interval if defined($update_interval);
        $args{require_qc_files} = $require_qc_files;
        $args{desired_qc_files} = $desired_qc_files;
        
        my @element_args;
        my $did = $self->_datasource_id;
        foreach my $result ($self->_all_files(%args, add_metadata_from_warehouse => 1, $required_metadata ? (required_metadata => $required_metadata) : (), $vrtrack_group ? (vrtrack_group => $vrtrack_group) : (), $graph_filter ? (graph_filter => $graph_filter, filter_after_grouping => 0) : ())) {
            my $protocol   = $result->{protocol};
            my $irods_path = $result->{irods_path};
            push(@element_args, { datasource => $did, result => { paths => $result->{paths}, $protocol ? (protocol => $protocol) : (), $irods_path ? (irods_path => $irods_path) : () } });
        }
        warn "all_with_warehouse_metadata called _all_files and will now _create_elements for ", scalar(@element_args), " elements\n" if $self->debug;
        $self->_create_elements(\@element_args);
    }
    
    method group_by_metadata_with_warehouse_metadata (Defined :$handle!, Str :$metadata_keys!, Str :$file_query!, Str|Dir :$local_root_dir?, Str :$update_interval?, Str :$required_metadata?, Str :$vrtrack_group?, Str :$graph_filter?, Bool :$filter_after_grouping = 1, Bool :$require_qc_files = 0, Str :$desired_qc_files = '_F0xB00.stats,.genotype.json,.verify_bam_id.json') {
        my %args;
        $args{handle}           = $handle;
        $args{file_query}       = $file_query;
        $args{local_root_dir}   = $local_root_dir if $local_root_dir;
        $args{update_interval}  = $update_interval if defined($update_interval);
        $args{require_qc_files} = $require_qc_files;
        $args{desired_qc_files} = $desired_qc_files;
        
        my @meta_keys = split /\|/, $metadata_keys;
        my $group_hash;
        foreach my $result ($self->_all_files(%args, add_metadata_from_warehouse => 1, $required_metadata ? (required_metadata => $required_metadata) : (), $vrtrack_group ? (vrtrack_group => $vrtrack_group) : (), $graph_filter ? (graph_filter => $graph_filter, filter_after_grouping => $filter_after_grouping) : ())) {
            my @group_keys;
            foreach my $key (@meta_keys) {
                $self->throw("Metadata key $key not present for file " . $result->{paths}->[0]) unless (exists $result->{metadata}->{$key});
                my $val = $result->{metadata}->{$key};
                if (ref $val) {
                    $val = join(',', sort @$val);
                }
                push @group_keys, $val;
            }
            my $group_key = join '|', @group_keys;
            push(@{ $group_hash->{$group_key}->{paths} }, @{ $result->{paths} });
            $group_hash->{$group_key}->{protocol} = $result->{protocol} if defined $result->{protocol};
            if ($graph_filter && $filter_after_grouping && $result->{pass_filter}) {
                $group_hash->{$group_key}->{passes_filter} = 1;
            }
        }
        
        my $did = $self->_datasource_id;
        my @element_args;
        while (my ($group, $data) = each %$group_hash) {
            if ($graph_filter && $filter_after_grouping) {
                next unless exists $group_hash->{$group}->{passes_filter};
            }
            my $protocol = $data->{protocol};
            push(@element_args, { datasource => $did, result => { paths => $data->{paths}, group => $group, $protocol ? (protocol => $protocol) : () } });
        }
        $self->_create_elements(\@element_args);
    }
    
    method _all_files (Defined :$handle!, Str :$file_query!, Str|Dir :$local_root_dir?, Str :$update_interval?, Bool :$add_metadata_from_warehouse?, Str :$required_metadata?, Str :$vrtrack_group?, Str :$graph_filter?, Bool :$filter_after_grouping = 1, Bool :$require_qc_files = 0, Str :$desired_qc_files = '_F0xB00.stats,.genotype.json,.verify_bam_id.json') {
        # _get_irods_files_and_metadata will get called twice in row: once to
        # see if the datasource changed, and again here; _has_changed caches
        # the result, and we clear the cache after getting that data.
        $add_metadata_from_warehouse ||= 0;
        my $files = $self->_get_irods_files_and_metadata($handle, $file_query, $local_root_dir || '', $add_metadata_from_warehouse || 0, $required_metadata || '', $vrtrack_group || '', $require_qc_files, $desired_qc_files, $graph_filter, $filter_after_grouping);
        $self->_clear_cache;
        
        my $did = $self->_datasource_id;
        my @results;
        my $anti_repeat_store = {};
        foreach my $path (sort { $files->{$a}->{vrpipe_irods_order} <=> $files->{$b}->{vrpipe_irods_order} } keys %$files) {
            my $new_metadata = $files->{$path};
            delete $new_metadata->{vrpipe_irods_order};
            
            my $pass_filter = delete $new_metadata->{irods_datasource_graph_filter};
            if ($graph_filter) {
                next unless defined($pass_filter);
            }
            
            my ($file_abs_path, $protocol);
            if ($local_root_dir) {
                my $sub_path = $path;
                $sub_path =~ s/^\///;
                $file_abs_path = file($local_root_dir, $sub_path)->stringify;
            }
            else {
                $file_abs_path = $path;
                $protocol      = 'irods:';
            }
            
            # consider type to be any if not defined in the irods metadata; if
            # not a VRPipe filetype it will be treated as an any
            my $type = delete $new_metadata->{type};
            $type ||= 'any';
            
            my $vrfile = VRPipe::File->create(path => $file_abs_path, type => $type, $protocol ? (protocol => $protocol) : ())->original;
            
            # add metadata to file, detecting any changes
            my @changed_details;
            my $changes = delete $new_metadata->{vrpipe_meta_changed};
            if ($changes) {
                foreach my $change (@$changes) {
                    push(@changed_details, "$file_abs_path $change->[0]: $change->[3]");
                    last;
                }
            }
            
            if ($add_metadata_from_warehouse && exists $new_metadata->{lane}) {
                my $lane = $vrfile->basename;
                $lane =~ s/\.gz$//;
                $lane =~ s/\.[^\.]+$//;
                $new_metadata->{lane} = $lane;
            }
            
            # if there was no metadata this will add metadata to the file.
            $vrfile->add_metadata($new_metadata, replace_data => 0);
            $vrfile->add_metadata({ $protocol ? () : (irods_path => $path), defined $new_metadata->{md5} ? (expected_md5 => $new_metadata->{md5}) : () });
            
            my $result_hash = { paths => [$file_abs_path], $protocol ? (protocol => $protocol) : (irods_path => $path) };
            if (@changed_details) {
                $result_hash->{changed} = [[$vrfile, $new_metadata]];
                $self->_start_over_elements_due_to_file_metadata_change($result_hash, \@changed_details, $anti_repeat_store);
                delete $result_hash->{changed};
            }
            push(@results, { paths => [$file_abs_path], $protocol ? (protocol => $protocol) : (irods_path => $path), metadata => $new_metadata, defined($pass_filter) ? (pass_filter => $pass_filter) : () });
        }
        
        $self->_clear_file_filter_cache();
        
        return @results;
    }
}

1;
