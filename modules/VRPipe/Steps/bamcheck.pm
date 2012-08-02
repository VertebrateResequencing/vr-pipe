
=head1 NAME

VRPipe::Steps::bamcheck - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2012 Genome Research Limited.

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

class VRPipe::Steps::bamcheck with VRPipe::StepRole {
    use VRPipe::Parser;
    
    method options_definition {
        return { bamcheck_exe => VRPipe::StepOption->create(description   => 'path to your bamcheck executable',
                                                            optional      => 1,
                                                            default_value => 'bamcheck'),
                 bamcheck_options => VRPipe::StepOption->create(description => 'options to bamcheck, excluding -r and -t (which are set by reference_fasta and exome_targets_file options)',
                                                                optional    => 1),
                 reference_fasta    => VRPipe::StepOption->create(description => 'absolute path to genome reference file'),
                 exome_targets_file => VRPipe::StepOption->create(description => 'absolute path to a file describing the targets/baits used for exome pulldown (tab-delimited [chr,start,end], where start is 1-based, and end is inclusive)',
                                                                  optional    => 1) };
    }
    
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'bam', description => 'bam files', max_files => -1) };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            
            my $options      = $self->options;
            my $bamcheck_exe = $options->{bamcheck_exe};
            my $opts         = VRPipe::Steps::bamcheck->get_bamcheck_options($options);
            
            my $req = $self->new_requirements(memory => 500, time => 1);
            foreach my $bam_file (@{ $self->inputs->{bam_files} }) {
                my $ifile      = $bam_file->path;
                my $check_file = $self->output_file(output_key => 'bamcheck_files', basename => $ifile->basename . '.bamcheck', type => 'txt');
                my $ofile      = $check_file->path;
                $self->dispatch_wrapped_cmd('VRPipe::Steps::bamcheck', 'stats_from_bamcheck', ["$bamcheck_exe $opts $ifile > $ofile", $req, { output_files => [$check_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return { bamcheck_files => VRPipe::StepIODefinition->create(type        => 'txt',
                                                                    description => 'the output of bamcheck on a bam',
                                                                    max_files   => -1,
                                                                    metadata    => {
                                                                                  source_bam => 'path to the bam file this bamcheck file was created from',
                                                                                  lane       => 'lane name (a unique identifer for this sequencing run, aka read group)' }) };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Creates a bamcheck file, a file of summary stats on a bam file, and also associates some of the stats as metadata on the bam file";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method get_bamcheck_options (ClassName|Object $self: HashRef $options!) {
        my $ref = file($options->{reference_fasta});
        $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
        my $opts = '-r ' . $ref;
        
        my $targets = $options->{exome_targets_file};
        if ($targets) {
            $targets = file($targets);
            $self->throw("exome_targets_file must be an absolute path") unless $targets->is_absolute;
            $opts .= ' -t ' . $targets;
        }
        
        my $user_opts = $options->{bamcheck_options};
        if ($user_opts) {
            if ($user_opts =~ /-r|-t/) {
                $self->throw("neither -r nor -t should be supplied as bamcheck options");
            }
            
            $opts .= ' ' . $user_opts;
        }
        
        return $opts;
    }
    
    method stats_from_bamcheck (ClassName|Object $self: Str $cmd_line) {
        my ($bam_path, $check_path) = $cmd_line =~ / (\S+) > (\S+)$/;
        $bam_path || $self->throw("bad cmd line [$cmd_line]");
        my $bam_file   = VRPipe::File->get(path => $bam_path);
        my $check_file = VRPipe::File->get(path => $check_path);
        
        $bam_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $check_file->update_stats_from_disc(retries => 3);
        if ($check_file->s) {
            my $new_meta = {};
            
            # parse the bamcheck file
            my $parser = VRPipe::Parser->create('bamcheck', { file => $check_file });
            $check_file->disconnect;
            
            my $hash_key_prefix = '';
            if ($cmd_line =~ /-d/) {
                $hash_key_prefix = 'rmdup_';
            }
            elsif ($cmd_line =~ /-t/) {
                $hash_key_prefix = 'targeted_';
            }
            
            # basic stats are stored under different prefixes, because we can't
            # do something like have 'reads' metadata be only for targeted
            # regions, since that is used when processing bams to ensure no
            # truncation
            $new_meta->{ $hash_key_prefix . 'reads' }            = $parser->sequences;
            $new_meta->{ $hash_key_prefix . 'reads_mapped' }     = $parser->reads_mapped;
            $new_meta->{ $hash_key_prefix . 'bases' }            = $parser->total_length;
            $new_meta->{ $hash_key_prefix . 'bases_mapped' }     = $parser->bases_mapped;
            $new_meta->{ $hash_key_prefix . 'bases_mapped_c' }   = $parser->bases_mapped_cigar;
            $new_meta->{ $hash_key_prefix . 'bases_trimmed' }    = $parser->bases_trimmed;
            $new_meta->{ $hash_key_prefix . 'reads_paired' }     = $parser->reads_paired;
            $new_meta->{ $hash_key_prefix . 'paired' }           = $parser->is_paired;
            $new_meta->{ $hash_key_prefix . 'error_rate' }       = $parser->error_rate;
            $new_meta->{ $hash_key_prefix . 'forward_reads' }    = $parser->first_fragments;
            $new_meta->{ $hash_key_prefix . 'reverse_reads' }    = $parser->last_fragments;
            $new_meta->{ $hash_key_prefix . 'avg_read_length' }  = $parser->average_length;
            $new_meta->{ $hash_key_prefix . 'mean_insert_size' } = $parser->insert_size_average;
            $new_meta->{ $hash_key_prefix . 'sd_insert_size' }   = $parser->insert_size_standard_deviation;
            
            unless ($hash_key_prefix eq 'rmdup_') {
                # when not using -d we can still add some rmdup stats
                my $reads_duplicated = $parser->reads_duplicated;
                $new_meta->{ $hash_key_prefix . 'rmdup_reads' }        = $new_meta->{ $hash_key_prefix . 'reads' } - $reads_duplicated;
                $new_meta->{ $hash_key_prefix . 'rmdup_reads_mapped' } = $new_meta->{ $hash_key_prefix . 'reads_mapped' } - $reads_duplicated;
                my $bases_duplicated = $parser->bases_duplicated;
                $new_meta->{ $hash_key_prefix . 'rmdup_bases' }        = $new_meta->{ $hash_key_prefix . 'bases' } - $bases_duplicated;
                $new_meta->{ $hash_key_prefix . 'rmdup_bases_mapped' } = $new_meta->{ $hash_key_prefix . 'bases_mapped' } - $bases_duplicated;
                
                # coverage stats we just don't bother having if in -d mode
                $new_meta->{ $hash_key_prefix . 'mean_coverage' } = $parser->mean_coverage;
                foreach my $cov (1, 2, 5, 10, 20, 50, 100) {
                    $new_meta->{ $hash_key_prefix . "bases_of_${cov}X_coverage" } = $parser->cumulative_coverage($cov);
                }
            }
            
            # and get other metadata from bam header, but don't overwrite
            # existing info
            $parser = VRPipe::Parser->create('bam', { file => $bam_file });
            my %rg_info       = $parser->readgroup_info();
            my @rgs           = keys %rg_info;
            my $existing_meta = $bam_file->metadata;
            if (@rgs == 1) {
                my $info = $rg_info{ $rgs[0] };
                unless (defined $existing_meta->{lane}) { $new_meta->{lane} = $info->{PU} ? $info->{PU} : $rgs[0]; }
                else                                    { $new_meta->{lane} = $existing_meta->{lane} }
                unless (defined $existing_meta->{library}) { $new_meta->{library} = $info->{LB} if $info->{LB}; }
                else                                       { $new_meta->{library} = $existing_meta->{library} }
                unless (defined $existing_meta->{sample}) { $new_meta->{sample} = $info->{SM} if $info->{SM}; }
                else                                      { $new_meta->{sample} = $existing_meta->{sample} }
                unless (defined $existing_meta->{center_name}) { $new_meta->{center_name} = $info->{CN} if $info->{CN}; }
                else                                           { $new_meta->{center_name} = $existing_meta->{center_name} }
                unless (defined $existing_meta->{platform}) { $new_meta->{platform} = $info->{PL} if $info->{PL}; }
                else                                        { $new_meta->{platform} = $existing_meta->{platform} }
                unless (defined $existing_meta->{study}) { $new_meta->{study} = $info->{DS} if $info->{DS}; }
                else                                     { $new_meta->{study} = $existing_meta->{study} }
            }
            else {
                # Add metadata if same for all RGs
                my %rg_items = ( PU => 'lane', LB => 'library', SM => 'sample', CN => 'center_name', PL => 'platform', DS => 'study',);

                foreach my $rg_k (keys %rg_items ) {
                    my $rg_item = $rg_items{$rg_k};
                    unless (defined $existing_meta->{$rg_item}) {
                        foreach my $i (0..$#rgs) {
                            my $info = $rg_info{ $rgs[$i] };
                            if ($info->{$rg_k}) {
                                if ($new_meta->{$rg_item}) {
                                    if ($info->{$rg_k} ne $new_meta->{$rg_item}) {
                                        delete $new_meta->{$rg_item};
                                        last;
                                    }
                                }
                                else {
                                    $new_meta->{$rg_item} = $info->{$rg_k};
                                }
                            }
                        }
                    }
                }
            }

            if (@rgs != 1 || "$new_meta->{lane}" eq "1") {
                # call the name something we can be most sure is maximally
                # unique
                my $md5 = $bam_file->md5;
                unless ($md5) {
                    $bam_file->update_md5;
                    $md5 = $bam_file->md5;
                }
                $new_meta->{lane} = $md5;
            }
            
            $bam_file->add_metadata($new_meta);
            $check_file->add_metadata({ source_bam => $bam_file->path->stringify,
                                        lane       => $bam_file->metadata->{lane} });
        }
        else {
            $self->throw("$check_path failed to be made");
        }
    }
}

1;
