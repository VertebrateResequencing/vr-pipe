
=head1 NAME

VRPipe::StepRole - methods required of and useful for writing Steps

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

This Role must be used with all C<VRPipe::Steps::*> modules. It also provides a
host of useful methods that Steps should take advantage of in their
C<body_sub>s and elsewhere.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2015 Genome Research Limited.

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

role VRPipe::StepRole {
    use Digest::MD5;
    use VRPipe::StepStatsUtil;
    use VRPipe::Persistent::InMemory;
    use VRPipe::Schema;
    use Cwd 'abs_path';
    use File::stat;
    use Fcntl qw(:mode);
    
    our $graph;
    our $schema;
    
    method name {
        my $class = ref($self);
        my ($name) = $class =~ /.+::(.+)/;
        return $name;
    }
    requires 'options_definition';
    requires 'inputs_definition';
    requires 'body_sub';
    requires 'post_process_sub';
    requires 'outputs_definition';
    requires 'description';
    requires 'max_simultaneous';
    
    # these may be needed by body_sub and post_process_sub
    has 'step_state' => (
        is      => 'rw',
        isa     => 'VRPipe::StepState',
        trigger => \&_rebuild
    );
    
    has 'data_element' => (
        is      => 'ro',
        isa     => 'VRPipe::DataElement',
        builder => '_build_data_element',
        lazy    => 1,
        clearer => '_clear_data_element'
    );
    
    has 'output_root' => (
        is      => 'ro',
        isa     => Dir,
        coerce  => 1,
        builder => '_build_output_root',
        lazy    => 1,
        clearer => '_clear_output_root'
    );
    
    has 'options' => (
        is      => 'ro',
        isa     => 'HashRef',
        builder => '_resolve_options',
        lazy    => 1,
        clearer => '_clear_options'
    );
    
    has 'inputs' => (
        is      => 'ro',
        isa     => PersistentFileHashRef,
        builder => '_resolve_inputs',
        lazy    => 1,
        clearer => '_clear_inputs'
    );
    
    has 'outputs' => (
        is      => 'ro',
        isa     => PersistentFileHashRef,
        builder => '_build_outputs',
        lazy    => 1,
        clearer => '_clear_outputs'
    );
    has 'temps' => (
        is      => 'ro',
        isa     => PersistentArrayRef,
        builder => '_build_temps',
        lazy    => 1,
        clearer => '_clear_temps'
    );
    
    has 'previous_step_outputs' => (
        is      => 'rw',
        isa     => PreviousStepOutput,
        trigger => sub { shift->_clear_inputs }
    );
    
    has 'allow_smaller_recommended_requirements_override' => (
        is      => 'rw',
        isa     => 'Bool',
        lazy    => 1,
        builder => '_build_smaller_recommended_requirements_override'
    );
    
    has 'debug' => (
        is      => 'rw',
        isa     => 'Bool',
        default => 0
    );
    
    # to avoid memory leak we need to store our real self on our non-persistent version
    has '_persistent_step' => (
        is  => 'rw',
        isa => 'VRPipe::Step',
    );
    
    # when parse is called, we'll store our dispatched refs here
    has 'dispatched' => (
        is      => 'ro',
        traits  => ['Array'],
        isa     => 'ArrayRef',
        lazy    => 1,
        default => sub { [] },
        handles => {
            _dispatch      => 'push',
            num_dispatched => 'count'
        },
        writer  => '_set_dispatched',
        clearer => '_clear_dispatched'
    );
    
    # and we'll also store all the output files the body_sub makes
    has '_output_files' => (
        is      => 'ro',
        traits  => ['Hash'],
        isa     => 'HashRef',
        lazy    => 1,
        default => sub { {} },
        handles => { _remember_output_files => 'set' },
        writer  => '_set_output_files',
        clearer => '_clear_output_files'
    );
    has '_temp_files' => (
        is      => 'ro',
        traits  => ['Array'],
        isa     => 'ArrayRef',
        lazy    => 1,
        default => sub { [] },
        handles => { _remember_temp_file => 'push' },
        writer  => '_set_temp_files',
        clearer => '_clear_temp_files'
    );
    has '_last_output_dir' => (
        is      => 'rw',
        isa     => Dir,
        lazy    => 1,
        coerce  => 1,
        builder => '_build_last_output_dir',
        clearer => '_clear_last_output_dir'
    );
    
    method _rebuild {
        $self->_clear_data_element;
        $self->_clear_output_root;
        $self->_clear_options;
        $self->_clear_inputs;
        $self->_clear_outputs;
        $self->_clear_temps;
        $self->_clear_dispatched;
        $self->_clear_output_files;
        $self->_clear_temp_files;
        $self->_clear_last_output_dir;
    }
    
    method _build_data_element {
        my $step_state = $self->step_state || $self->throw("Cannot get data element without step state");
        return $step_state->dataelement;
    }
    
    method _build_output_root {
        my $step_state = $self->step_state || $self->throw("Cannot get output root without step state");
        my $pipeline_root = $step_state->pipelinesetup->output_root;
        
        my $de_id          = $step_state->dataelement->id;
        my $des_id         = VRPipe::DataElementState->get(dataelement => $de_id, pipelinesetup => $step_state->pipelinesetup)->id;
        my $hashing_string = 'VRPipe::DataElementState::' . $des_id;
        my @subdirs        = $self->hashed_dirs($hashing_string);
        
        my $dir = dir($pipeline_root, @subdirs, $de_id, $step_state->stepmember->step_number . '_' . $self->name);
        $self->make_path($dir);
        return $dir;
    }
    
    method _build_last_output_dir {
        return $self->output_root;
    }
    
    method _build_outputs {
        my $step_state = $self->step_state || $self->throw("Cannot get outputs without step state");
        if ($step_state->same_submissions_as) {
            $step_state = $step_state->same_submissions_as;
        }
        my $outs  = $step_state->output_files;
        my $temps = $step_state->temp_files;
        foreach my $temp (@$temps) {
            push(@{ $outs->{temp} }, $temp);
        }
        return $outs;
    }
    
    method _build_temps {
        my $step_state = $self->step_state || $self->throw("Cannot get outputs without step state");
        return $step_state->temp_files;
    }
    
    method _build_smaller_recommended_requirements_override {
        return 1;
    }
    
    method _resolve_options {
        my $step_state = $self->step_state || $self->throw("Cannot get options without step state");
        my $user_opts  = $step_state->pipelinesetup->options;
        my $hash       = $self->options_definition;
        
        my %return;
        while (my ($key, $val) = each %$hash) {
            if ($val->isa('VRPipe::StepOption')) {
                my $user_val = $user_opts->{$key};
                if (defined $user_val) {
                    my $allowed = $val->allowed_values;
                    if (@$allowed) {
                        my %allowed = map { $_ => 1 } @$allowed;
                        if (exists $allowed{$user_val}) {
                            $return{$key} = $user_val;
                        }
                        else {
                            $self->throw("'$user_val' is not an allowed option for '$key'");
                        }
                    }
                    else {
                        $return{$key} = $user_val;
                    }
                }
                else {
                    my $default  = $val->default_value;
                    my $optional = $val->optional;
                    if (defined $default && $optional ? 1 : $default ne '') {
                        $return{$key} = $default;
                    }
                    elsif (!$optional) {
                        $self->throw("the option '$key' is required");
                    }
                }
            }
            else {
                $self->throw("invalid class " . ref($val) . " supplied for option '$key' definition");
            }
        }
        
        # add in requirements overrides that pertain to us
        my $name = $self->name;
        foreach my $resource (qw(memory time cpus tmp_space local_space custom)) {
            my $user_key = $name . '_' . $resource;
            if (defined $user_opts->{$user_key}) {
                $return{ $resource . '_override' } = $user_opts->{$user_key};
            }
        }
        
        return \%return;
    }
    
    method _resolve_inputs {
        my $hash     = $self->inputs_definition;
        my $step_num = $self->step_state->stepmember->step_number;
        my ($step_adaptor) = VRPipe::StepAdaptor->search({ pipeline => $self->step_state->stepmember->pipeline->id, to_step => $step_num });
        return {} unless $step_adaptor;
        
        my $im = VRPipe::Persistent::InMemory->new();
        my %return;
        keys %$hash;      # reset the iterator in case we threw in a previous pass
        while (my ($key, $val) = each %$hash) {
            if ($val->isa('VRPipe::File')) {
                $return{$key} = [$val->e ? $val : $val->resolve];
            }
            elsif ($val->isa('VRPipe::StepIODefinition')) {
                # see if we have this $key in our previous_step_outputs or
                # via the options or data_element
                my $results;
                
                my $pso = $self->previous_step_outputs;
                if ($pso) {
                    # can our StepAdaptor adapt previous step output to this
                    # input key?
                    $results = $step_adaptor->adapt(input_key => $key, pso => $pso);
                }
                if (!$results) {
                    # can our StepAdaptor translate our dataelement into a file
                    # for this key?
                    $results = $step_adaptor->adapt(input_key => $key, data_element => $self->data_element);
                }
                if (!$results) {
                    my $opts = $self->options;
                    if ($opts && defined $opts->{$key}) {
                        $results = [VRPipe::File->get(path => $opts->{$key})];
                    }
                }
                
                if (!$results) {
                    if ($val->min_files == 0) {
                        next;
                    }
                    else {
                        $self->throw("the input file(s) for '$key' of stepstate " . $self->step_state->id . " could not be resolved");
                    }
                }
                
                # but only bother caching if we have more than 5 files
                # *** this is mainly a hack so that if we delete a file and then
                # need it again, start_over will result in the file being
                # recreated - the hack gets gatk_realign_and_variant_call.t to
                # pass, but otherwise isn't the best solution...
                my $check_note = @$results > 5;
                
                my @vrfiles;
                my @skip_reasons;
                my $wanted_type = $val->type;
                my $all_good    = 'input_files_all_good.' . $wanted_type . '.' . join(',', map { $_->id } @$results);
                my $cached      = 0;
                if ($check_note) {
                    my $note_value = $im->noted($all_good);
                    if ($note_value) {
                        if ($note_value eq 'none') {
                            next;
                        }
                        
                        my %valid_ids = map { $_->id => $_ } @$results;
                        my $ok = 1;
                        my @files;
                        foreach my $id (split(/,/, $note_value)) {
                            if (exists $valid_ids{$id}) {
                                push(@files, delete $valid_ids{$id});
                            }
                            else {
                                $ok = 0;
                                last;
                            }
                        }
                        
                        if ($ok) {
                            @vrfiles = @files;
                            $im->note($all_good, value => $note_value); # to refresh the timeout
                            $cached = 1;
                        }
                    }
                }
                
                # it is expensive to VRPipe::FileType->create, and expensive to
                # check_type(), so we create the former once and cache the final
                # result of the latter if all passed
                my $type = VRPipe::FileType->create($wanted_type, { file => 'to_be_replaced' });
                
                unless ($cached) {
                    my $skipped_not_correct_type = 0;
                    foreach my $result (@$results) {
                        unless (ref($result) && ref($result) eq 'VRPipe::File') {
                            $result = VRPipe::File->get(path => file($result)->absolute);
                        }
                        
                        unless ($wanted_type eq 'any') {
                            my $resolved = $result->resolve;
                            
                            my $has_size = 0;
                            if ($resolved->s) {
                                $has_size = 1;
                                $type->file($resolved->path);
                                unless ($type->check_type) {
                                    # the type check will fail if the file doesn't
                                    # really exist, so make sure
                                    $resolved->update_stats_from_disc;
                                    if ($resolved->s) {
                                        push(@skip_reasons, "file " . $result->path . " was not the correct type ($wanted_type)");
                                        $skipped_not_correct_type++;
                                        next;
                                    }
                                    else {
                                        $result->update_stats_from_disc;
                                        $has_size = 0;
                                    }
                                }
                            }
                            unless ($has_size) {
                                my $db_type = $resolved->type;
                                
                                if ($db_type) {
                                    # if that's an auto-generated type then it
                                    # should be treated as an 'any' for this
                                    # purpose
                                    eval "require VRPipe::FileType::$db_type;";
                                    if ($@) {
                                        $db_type = 'any';
                                    }
                                }
                                
                                if ($db_type && $wanted_type ne $db_type && $db_type ne 'any') {
                                    push(@skip_reasons, "file " . $result->path . " does not exist so its type can't be checked properly, but it doesn't seem to be correct: expected type $wanted_type and got type $db_type");
                                    next;
                                }
                            }
                        }
                        elsif (!$result->e) {
                            # missing_input_files() assumes we've already done
                            # the confirmation of file existence
                            my $resolved = $result->resolve;
                            unless ($resolved->e) {
                                $resolved->update_stats_from_disc;
                            }
                        }
                        
                        push(@vrfiles, $result);
                    }
                    if ($check_note && (@vrfiles || $val->min_files == 0) && ((@vrfiles == @$results) || (@vrfiles + $skipped_not_correct_type == @$results))) {
                        $im->note($all_good, value => @vrfiles ? join(',', map { $_->id } @vrfiles) : 'none');
                    }
                }
                
                my $max_allowed = $val->max_files;
                my $min_allowed = $val->min_files;
                if (!@vrfiles && @skip_reasons) {
                    $self->throw("for '$key' of stepstate " . $self->step_state->id . " none of the input files had a suitable type:\n" . join("\n", @skip_reasons)) if ($min_allowed > 0);
                    next;
                }
                my $num_results = @vrfiles;
                if ($max_allowed == -1) {
                    $max_allowed = $num_results;
                }
                if ($min_allowed == -1) {
                    $min_allowed = $num_results;
                }
                unless ($num_results >= $min_allowed && $num_results <= $max_allowed) {
                    $self->throw("there were $num_results input file(s) for '$key' of stepstate " . $self->step_state->id . ", which does not fit the allowed range $min_allowed..$max_allowed");
                }
                
                $return{$key} = [map { $_->e ? $_ : $_->resolve } @vrfiles];
            }
            else {
                $self->throw("invalid class " . ref($val) . " supplied for input '$key' value definition");
            }
        }
        
        return \%return;
    }
    
    method _missing (PersistentFileHashRef $hash, PersistentHashRef $defs, Bool $inputs_mode = 0) {
        # (in inputs_mode we do not have to check file existence and type since
        # inputs() already did that)
        my $debug = $self->debug;
        warn "       in _missing()\n" if $debug;
        my $t = time();
        
        my %file_type_objs;
        
        my (@missing, @messages, %missing_metadata);
        # check the files we actually need/output are as expected
        while (my ($key, $val) = each %$hash) {
            my $def     = $defs->{$key};
            my $check_s = 1;
            my @required_metadata_keys;
            if ($def && $def->isa('VRPipe::StepIODefinition')) {
                $check_s                = $def->check_existence;
                @required_metadata_keys = $def->required_metadata_keys;
            }
            
            foreach my $file (@$val) {
                my $resolved = $inputs_mode ? $file : $file->resolve;
                
                my $ignore_s = 0;
                if ($check_s && $file->protocol ne 'file:/') {
                    # we can't check the existence of non-local files
                    $ignore_s = 1;
                }
                
                # we may be in a situation where $file has been moved elsewhere
                # in the past, but now we've recreated a possibly different
                # $file, so $resolved is out-of-date (or possibly doesn't exist
                # any more)
                if (!$ignore_s && $check_s && !$inputs_mode) {
                    # double-check incase the step did not update_stats_from_disc
                    my $actual_s   = $resolved->check_file_size_on_disc();
                    my $resolved_s = $resolved->s;
                    if ((defined $actual_s && defined $resolved_s && $actual_s != $resolved_s) || (!defined $actual_s && $resolved_s) || (defined $actual_s && !defined $resolved_s)) {
                        $resolved->update_stats_from_disc(keep_stat_cache => 1);
                    }
                    $file->update_stats_from_disc(keep_stat_cache => 1) unless $resolved->id == $file->id;
                }
                if ($file->e && (!$resolved->s || $file->mtime > $resolved->mtime)) {
                    $file->moved_to(undef);
                    $file->update;
                    $resolved = $file;
                }
                
                if (!$ignore_s && $check_s && !$resolved->s) {
                    push(@missing, $file->path);
                    push(@messages, $file->path . ($resolved->e ? " is an empty file." : " does not exist."));
                }
                else {
                    my $bad = 0;
                    
                    # check the filetype is correct (except for temp files,
                    # since the check is expensive)
                    if (!$ignore_s && $check_s && !$inputs_mode && $key ne 'temp') {
                        my $type_str = $resolved->type;
                        unless (exists $file_type_objs{$type_str}) {
                            $file_type_objs{$type_str} = VRPipe::FileType->create($type_str, { file => 'to_be_replaced' });
                        }
                        my $type = $file_type_objs{$type_str};
                        $type->file($resolved->path);
                        unless ($type->check_type) {
                            push(@messages, $resolved->path . " exists, but is the wrong type (not a " . $resolved->type . " file)!");
                            $bad = 1;
                        }
                    }
                    
                    # check the expected metadata keys exist
                    if (@required_metadata_keys) {
                        my $meta = $file->metadata;
                        foreach my $key (@required_metadata_keys) {
                            unless (exists $meta->{$key}) {
                                push(@messages, $file->path . " exists, but lacks required metadata key $key!");
                                $bad = 1;
                                $missing_metadata{ $file->path } = 1;
                            }
                        }
                    }
                    
                    if ($bad) {
                        push(@missing, $file->path);
                    }
                }
            }
        }
        
        my $e = time() - $t;
        warn "       _missing() returning after $e seconds\n" if $debug;
        
        return (\@missing, \@messages, \%missing_metadata);
    }
    
    method missing_input_files {
        my $debug = $self->debug;
        warn "      missing_input_files() will get inputs() and definitions\n" if $debug;
        my $t    = time();
        my $hash = $self->inputs;
        my $e1   = time() - $t;
        $t = time();
        my $defs = $self->inputs_definition;
        my $e2   = time() - $t;
        warn "      inputs() took $e1 seconds and inputs_definition() took $e2 seconds; will call _missing()\n" if $debug;
        return $self->_missing($hash, $defs, 1);
    }
    
    method missing_output_files {
        my $debug = $self->debug;
        warn "      missing_output_files() will get outputs() and definitions\n" if $debug;
        my $hash = $self->outputs;
        my $defs = $self->outputs_definition;
        
        # check that we don't have any outputs defined in the definition that
        # no files were made for (we must not use while each for the defs hash,
        # since $defs is the same for multiple calls of this method for
        # different dataelementstates in piplinesetup trigger(), and if one
        # fails with a throw, the others would not enter the loop and look like
        # they had no missing output files and instantly complete, even though
        # they failed horribly)
        foreach my $key (keys %$defs) {
            next if exists $hash->{$key};
            my $val = $defs->{$key};
            next if $val->min_files == 0;
            
            # we used to just throw here, as it would indicate a permanent error
            # in the step body_sub code. However it is also possible for this
            # to happen due to a bad database update where the StepOutputFile
            # row does not get created even though the Submissions do. For that
            # case we risk an infinite restart situation and restart this step
            return ([1], ["'$key' was defined as an output of step " . $self->step_state->stepmember->step_number . " (" . $self->step_state->stepmember->step->name . "), yet no output file was made with that output_key (either the Step is not correctly written, or the database failed to update correctly when the Step was parsed)"]);
        }
        
        warn "      , will call _missing()\n" if $debug;
        return $self->_missing($hash, $defs);
    }
    
    method output_file (Str :$output_key?, File|Str :$basename!, FileType :$type!, Dir|Str :$output_dir?, Dir|Str :$sub_dir?, HashRef :$metadata?, Str :$protocol?, Bool :$temporary = 0) {
        if (!$temporary && !$output_key) {
            $self->throw("output_key is required");
        }
        if (!$temporary && $output_key eq 'temp') {
            $self->throw("'temp' is a reserved output_key");
        }
        
        $output_dir ||= $self->output_root;
        $output_dir = dir($output_dir);
        if ($sub_dir) {
            $output_dir = dir($output_dir, $sub_dir);
        }
        $self->throw("output_dir must be absolute ($output_dir)") unless $output_dir->is_absolute;
        if (!$protocol || ($protocol eq 'file:/')) {
            my $lod = $self->_last_output_dir;
            if (!defined $lod || ($lod ne $output_dir)) {
                $self->make_path($output_dir);
                $self->_last_output_dir($output_dir);
            }
        }
        
        my $vrfile = VRPipe::File->create(path => file($output_dir, $basename), type => $type, $protocol ? (protocol => $protocol) : ());
        $vrfile->add_metadata($metadata) if $metadata;
        
        if ($temporary) {
            my $root = $self->output_root;
            $self->throw("temporary files must be placed within the output_root") unless $output_dir =~ /^$root/;
            $self->_remember_temp_file($vrfile);
        }
        else {
            my $hash = $self->_output_files;
            my $files = $hash->{$output_key} || [];
            push(@$files, $vrfile);
            $self->_remember_output_files($output_key => $files);
        }
        
        return $vrfile;
    }
    
    method set_cmd_summary (VRPipe::StepCmdSummary $cmd_summary) {
        my $step_state = $self->step_state;
        $step_state->cmd_summary($cmd_summary);
        $step_state->update;
    }
    
    method _run_coderef (Str $method_name) {
        # $self is a VRPipe::Step, even when *_sub was defined in a
        # VRPipe::Steps::subclass; to regain full benefits of inheritance in
        # those *_sub subs, we'll load the real module and use that instead
        # of $self
        my $non_persistent = $self->_from_non_persistent;
        my $return;
        if ($non_persistent) {
            #*** this is very ugly - is there a better way?
            $non_persistent->step_state($self->step_state);
            $non_persistent->previous_step_outputs($self->previous_step_outputs);
            $non_persistent->_persistent_step($self);
            $non_persistent->debug($self->debug);
            
            my $ref = $non_persistent->$method_name();
            $return = &$ref($non_persistent);
            
            if ($method_name eq 'body_sub') {
                $self->_set_output_files($non_persistent->_output_files);
                $self->_set_temp_files($non_persistent->_temp_files);
                $self->_last_output_dir($non_persistent->_last_output_dir);
                $self->_set_dispatched($non_persistent->dispatched);
            }
        }
        else {
            my $ref = $self->$method_name();
            $return = &$ref($self);
        }
        
        return $return;
    }
    
    method parse {
        # if we have missing input files, check to see if some other step
        # created them, and start those steps over in the hopes the files will
        # be recreated; otherwise throw
        my $debug = $self->debug;
        warn "    parse will check missing input files\n" if $debug;
        my ($missing, $messages, $missing_metadata) = $self->missing_input_files;
        if (@$missing) {
            my $with_recourse = 0;
            my %states_to_restart;
            foreach my $path (@$missing) {
                # there's no recourse if the file was actually just missing some
                # metadata, not physically missing
                next if exists $missing_metadata->{$path};
                
                my $file = VRPipe::File->get(path => $path);
                my $count = 0;
                my $state;
                foreach my $sof (VRPipe::StepOutputFile->search({ file => $file->id }, { prefetch => 'stepstate' })) {
                    $count++;
                    $state = $sof->stepstate;
                }
                
                if ($count == 1) {
                    $with_recourse++;
                    push(@{ $states_to_restart{ $state->id } }, $file->resolve->path);
                }
            }
            
            if ($with_recourse == @$missing) {
                while (my ($state_id, $files) = each %states_to_restart) {
                    #*** would have to make inactive pipeline active as well...
                    #    how to make sure Manager parses this step soon?...
                    my $state = VRPipe::StepState->get(id => $state_id);
                    if ($state->complete) {
                        $state->pipelinesetup->log_event("Calling StepState->start_over to regenerate needed input files (@$files)", stepstate => $state->id, dataelement => $state->dataelement->id);
                        $state->start_over;
                    }
                }
                return "Called start_over to regenerate needed input files";
            }
            else {
                my $step_state = $self->step_state;
                return "There is a problem with the input files for step " . $self->name . " (for data element " . $self->data_element->id . ", pipelinesetup " . $step_state->pipelinesetup->id . ", stepstate " . $step_state->id . "):\n" . join("\n", @$messages);
            }
        }
        else {
            warn "    - no missing files\n" if $debug;
        }
        
        # actually run the body_sub
        warn "    will run body_sub\n" if $debug;
        my $t = time();
        $self->_run_coderef('body_sub');
        my $e = time() - $t;
        warn "    ran the body_sub in $e seconds\n" if $debug;
        
        # store output and temp files on the StepState
        my $output_files = $self->_output_files;
        
        if (keys %$output_files) {
            my $step_state = $self->step_state;
            $step_state->output_files($output_files);
            $step_state->update;
            my @ofiles;
            while (my ($key, $files) = each %$output_files) {
                push(@ofiles, map { $_->path } @$files);
            }
            
            if (@ofiles < 100) {
                my %ofiles = map { $_ => 1 } @ofiles;
                # we could have so many files that attempting to store all their
                # paths in a single pipelinesetuplog would overflow the maximum
                # size for the message column in the database, which means our
                # parse would fail, so we create multiple log lines, one for each
                # output file
                foreach my $ofile (keys %ofiles) {
                    $step_state->pipelinesetup->log_event("StepRole->parse() ran the body_sub and created new StepOutputFile $ofile", stepstate => $step_state->id, dataelement => $step_state->dataelement->id);
                }
            }
            warn "    - stored ", scalar(keys %$output_files), " output files on the stepstate\n" if $debug;
        }
        else {
            warn "    - no output files\n" if $debug;
        }
        
        my $temp_files = $self->_temp_files;
        if (@$temp_files) {
            my $step_state = $self->step_state;
            $step_state->temp_files($temp_files);
            $step_state->update;
            warn "    - stored ", scalar(@$temp_files), " temp files on the stepstate\n" if $debug;
        }
        
        my $dispatched = $self->num_dispatched;
        if ($dispatched) {
            return 0;
        }
        else {
            # nothing dispatched, run the post_process and return any errors
            # from that
            warn "    will run post_process\n" if $debug;
            return $self->post_process;
        }
    }
    
    method post_process {
        my $debug = $self->debug;
        warn "     post_process() will run the post_process_sub\n" if $debug;
        my $ok        = $self->_run_coderef('post_process_sub');
        my $stepstate = $self->step_state;
        my $our_ss_id = $stepstate->id;
        
        my $error;
        if ($ok) {
            warn "     post_process_sub returned ok, will check missing_output_files\n" if $debug;
            my ($missing, $messages) = $self->missing_output_files;
            warn "     - will unlink_temp_files\n" if $debug;
            $stepstate->unlink_temp_files;
            if (@$missing) {
                $stepstate->pipelinesetup->log_event("Calling StepState->start_over because post_process had a problem with the output files: " . join("\n", @$messages), stepstate => $our_ss_id, dataelement => $stepstate->dataelement->id);
                $stepstate->start_over;
                $error = "There was a problem with the output files, so the stepstate was started over:\n" . join("\n", @$messages);
                warn "     - $error\n" if $debug;
            }
            else {
                warn "     no missing output files, will handle file ownership\n" if $debug;
                
                # change group ownership on files and enable world access to
                # all parent dirs if user set group on the setup
                my $ps    = $stepstate->pipelinesetup;
                my $group = $ps->unix_group;
                if ($group) {
                    my $hash = $self->outputs;
                    my @paths;
                    while (my ($key, $val) = each %$hash) {
                        foreach my $file (@$val) {
                            # but only do this if we were the first to create
                            # the files
                            my $oss = $file->output_by(1);
                            $oss || next;
                            next unless $oss->id == $our_ss_id;
                            
                            # later on we don't need to do anything to files
                            # that don't exist, and trying to stat a
                            # non-existant file causes us to die and possibly
                            # result in an infinite-restart loop!
                            $file->reselect_values_from_db;
                            $file->e || next;
                            
                            my $path = $file->path;
                            if (-l $path) {
                                my $real = abs_path($path);
                                if ($real ne $path) {
                                    my ($real_file) = VRPipe::File->search({ path => $real });
                                    if ($real_file) {
                                        $oss = $real_file->output_by(1);
                                        next unless ($oss && ($oss->id == $our_ss_id));
                                    }
                                }
                            }
                            
                            push(@paths, $path);
                        }
                    }
                    
                    my (undef, undef, $gid) = getgrnam($group);
                    if ($gid) {
                        warn "     will change permissions for ", scalar(@paths), " files\n" if $debug;
                        # change the group on the files
                        chown $<, $gid, @paths;
                        
                        # in case we are running as root we'll try and change
                        # the user of the files as well
                        my (undef, undef, $uid) = getpwnam($ps->user);
                        if ($uid && $uid != $<) {
                            chown $uid, $gid, @paths;
                            
                            # make sure we still have write access to the files,
                            # but don't chmod if files already are world
                            # readable
                            my $all_can_read = 1;
                            foreach my $path (@paths) {
                                my $mode = stat($path)->mode;
                                my $readable = ($mode & S_IRUSR) && ($mode & S_IRGRP) && ($mode & S_IROTH);
                                unless ($readable) {
                                    $all_can_read = 0;
                                    last;
                                }
                            }
                            
                            unless ($all_can_read) {
                                chmod 0660, @paths; # -rw-rw----
                            }
                        }
                        
                        # try and make parent dirs accessible
                        my %dirs;
                        foreach my $path (@paths) {
                            my $dir = file($path)->dir;
                            $dirs{$dir} = 1;
                            my $num_parents = $dir->dir_list;
                            for (1 .. $num_parents) {
                                $dir = $dir->parent;
                                $dirs{$dir} = 1;
                            }
                        }
                        chmod 0755, keys %dirs; # -rwxrwxr-x
                    }
                }
                
                return;
            }
        }
        else {
            $stepstate->unlink_temp_files;
            $stepstate->pipelinesetup->log_event("Calling StepState->start_over because post_process did not return true", stepstate => $our_ss_id, dataelement => $stepstate->dataelement->id);
            $stepstate->start_over;
            $error = "The post_process_sub did not return true, so did a start_over on the stepstate.";
            warn "     $error\n" if $debug;
        }
        
        return $error;
    }
    
    method new_requirements (Int :$memory!, Int :$time!, Int :$cpus?, Int :$tmp_space?, Int :$local_space?, HashRef :$custom?) {
        # time used to be specified in hours, but now in seconds
        if ($time < 60) {
            $time *= 60 * 60;
        }
        
        # get the current mean+2sd memory and time of past runs of this step
        my $ssu = VRPipe::StepStatsUtil->new(step => $self->isa('VRPipe::Step') ? $self : $self->_persistent_step);
        my ($rec_mem, $rec_time);
        my $ss = $self->step_state;
        if ($ss) {
            my $setup = $ss->pipelinesetup;
            $rec_mem = $ssu->recommended_memory(pipelinesetup => $setup);
            $rec_time = $ssu->recommended_time(pipelinesetup => $setup);
        }
        $rec_mem  ||= $ssu->recommended_memory;
        $rec_time ||= $ssu->recommended_time;
        
        # if we have recommendations, override the settings passed in from the
        # step body_sub
        my $allow_override = $self->allow_smaller_recommended_requirements_override;
        if ($rec_mem && ($allow_override ? 1 : $rec_mem > $memory)) {
            $memory = $rec_mem;
        }
        if ($rec_time && ($allow_override ? 1 : $rec_time > $time)) {
            $time = $rec_time;
        }
        
        # user can override by providing pipelinesetup options
        my $options = $self->options;
        if ($options->{memory_override} && $options->{memory_override} > $memory) {
            $memory = $options->{memory_override};
        }
        if ($options->{time_override}) {
            my $override_time = $options->{time_override};
            if ($override_time < 60) {
                $override_time *= 60 * 60;
            }
            if ($override_time > $time) {
                $time = $override_time;
            }
        }
        #*** and the other resources?...
        
        # round time up to the nearest xmins, depending on how much time is
        # needed
        my $rounder;
        if ($time < 900) {
            $rounder = 300; # 5mins
        }
        elsif ($time < 3600) {
            $rounder = 900; # 15mins
        }
        elsif ($time < 43200) {
            $rounder = 1800; # 30mins
        }
        else {
            $rounder = 3600; # 60mins
        }
        my $lower_bound = $time - ($time % $rounder);
        $time = $lower_bound + $rounder;
        
        # due to overheads of running things via vrpipe-handler, we have a
        # minimum memory req of 500MB and will increase memory by 250MB
        $memory += 250;
        if ($memory < 500) {
            $memory = 500;
        }
        
        my $req = VRPipe::Requirements->create(
            memory => $memory,
            time   => $time,
            $cpus        ? (cpus        => $cpus)        : (),
            $tmp_space   ? (tmp_space   => $tmp_space)   : (),
            $local_space ? (local_space => $local_space) : (),
            $custom      ? (custom      => $custom)      : ()
        );
        
        return $req;
    }
    
    method dispatch (ArrayRef $aref) {
        my $extra_args = $aref->[2] || {};
        $extra_args->{dir} ||= $self->_last_output_dir;
        $aref->[2] = $extra_args;
        $self->_dispatch($aref);
    }
    
    method dispatch_vrpipecode (Str $code, VRPipe::Requirements $req, HashRef $extra_args?) {
        my $deployment = VRPipe::Persistent::SchemaBase->database_deployment;
        
        # use lib for anything that has been added to INC, but only if we're
        # testing, since production cannot (must not) work with temp altered INC
        my $use_lib = '';
        my $deploy  = '';
        if ($deployment eq 'testing') {
            use lib;
            my %orig_inc = map { $_ => 1 } @lib::ORIG_INC;
            my %new_lib;
            foreach my $inc (@INC) {
                unless (exists $orig_inc{$inc}) {
                    $new_lib{$inc} = 1;
                }
            }
            $use_lib = join(' ', map { '-I' . dir($_)->absolute } ('modules', 't', keys %new_lib)) . ' ';
            $deploy = '=testing';
        }
        
        my $cmd = qq[perl $use_lib-MVRPipe$deploy -e "$code"];
        $self->dispatch([$cmd, $req, $extra_args]);
    }
    
    method dispatch_md5sum (VRPipe::File $vrfile, Maybe[Str] $expected_md5) {
        my $path = $vrfile->path;
        my $req = $self->new_requirements(memory => 500, time => 1);
        
        if ($expected_md5) {
            return $self->dispatch_vrpipecode(qq[use Digest::MD5; VRPipe::Manager->get->disconnect; open(FILE, q[$path]) or die q[Could not open file $path]; binmode(FILE); if (Digest::MD5->new->addfile(*FILE)->hexdigest eq q[$expected_md5]) { VRPipe::File->get(path => q[$path], md5 => q[$expected_md5]); } else { die q[md5sum of $path does not match expected value] }], $req);
        }
        else {
            return $self->dispatch_vrpipecode(qq[use Digest::MD5; VRPipe::Manager->get->disconnect; open(FILE, q[$path]) or die q[Could not open file $path]; binmode(FILE); VRPipe::File->get(path => q[$path], md5 => Digest::MD5->new->addfile(*FILE)->hexdigest);], $req);
        }
    }
    
    # using this lets you run a bunch of perl code wrapping a command line exe,
    # yet still keep the command line visible so you know the most important
    # bit of what was run
    method dispatch_wrapped_cmd (Str $class, Str $method, ArrayRef $dispatch_args) {
        my ($cmd, $req, $extra_args) = @$dispatch_args;
        
        eval "require $class;";
        $self->throw("'$class' is not a valid class name: $@") if $@;
        
        $class->can($method) || $self->throw("$method is not a valid method of $class");
        my $code = $class->isa('VRPipe::Persistent') ? '' : "use $class; ";
        $code .= "$class->$method(q[$cmd]);";
        
        my @args = ($code, $req);
        push(@args, $extra_args) if $extra_args;
        
        return $self->dispatch_vrpipecode(@args);
    }
    
    method common_metadata (ArrayRef['VRPipe::File'] $files!) {
        my %meta;
        foreach my $file (@$files) {
            my $file_meta = $file->metadata;
            while (my ($key, $val) = each %$file_meta) {
                if (ref($val)) {
                    $val = join(',!,', @$val);
                }
                $meta{$key}->{$val}++;
            }
        }
        
        # only keep metadata common to all files
        my $common_meta = {};
        foreach my $key (keys %meta) {
            my @vals = keys %{ $meta{$key} };
            next unless (@vals == 1 && $meta{$key}->{ $vals[0] } == @$files);
            my @val = split(/,!,/, $vals[0]);
            $common_meta->{$key} = @val > 1 ? \@val : $val[0];
        }
        
        return $common_meta;
    }
    
    method combined_metadata (ArrayRef['VRPipe::File'] $files!) {
        my %meta;
        foreach my $file (@$files) {
            my $file_meta = $file->metadata;
            while (my ($key, $val) = each %$file_meta) {
                if (exists $meta{$key}) {
                    my $current = $meta{$key};
                    if (ref($current) && ref($current) eq 'ARRAY') {
                        push(@{ $meta{$key} }, ref($val) && ref($val) eq 'ARRAY' ? @{$val} : $val);
                    }
                    else {
                        $meta{$key} = [$current, ref($val) && ref($val) eq 'ARRAY' ? @{$val} : $val];
                    }
                }
                else {
                    $meta{$key} = $val;
                }
            }
        }
        
        return \%meta;
    }
    
    method element_meta {
        my %element_meta = %{ $self->step_state->dataelement->result };
        delete @element_meta{qw(paths lane group)};
        return %element_meta;
    }
    
    method handle_override_options (HashRef $meta!) {
        my $options = $self->options;
        
        return $options unless (exists $meta->{chunk_override_file} && exists $meta->{chrom} && exists $meta->{from} && exists $meta->{to});
        
        my $override_options = $options;
        
        my $chrom  = $meta->{chrom};
        my $from   = $meta->{from};
        my $to     = $meta->{to};
        my $region = "${chrom}_${from}-${to}";
        
        my $override_file = $meta->{chunk_override_file};
        my $override      = do $override_file;
        
        foreach my $option (keys %$override_options) {
            if (exists $override->{"$region"}->{"$option"}) {
                $override_options->{"$option"} = $override->{"$region"}->{"$option"};
            }
        }
        
        return $override_options;
    }
    
    # a step can use normal VRPipe::Schema methods to create whatever graph
    # nodes it likes to store some results, but these will not be associated
    # with the stepstate (the pipelinesetup, datalement and step) that created
    # them; supplying the result nodes to this sub (either in the body_sub where
    # we created the Step as an instance with step_state passed to new(), or
    # while running as a Job from code made by dispatch_vrpipecode/wrapped_cmd
    # where we can find out the stepstate id from an env var) will give the
    # desired association
    method result_nodes (ClassName|Object $self: ArrayRef $nodes, Object $ss?) {
        $schema ||= VRPipe::Schema->create('VRPipe');
        
        unless ($ss) {
            if (defined $ENV{VRPIPE_STEPSTATE}) {
                $ss = VRPipe::StepState->get(id => $ENV{VRPIPE_STEPSTATE});
            }
            elsif (ref($self)) {
                $ss = $self->step_state;
            }
            else {
                return;
            }
        }
        
        # if and when we stop using MySQL and switch completely to the graph db,
        # there will already be a stepstate node in the graph connected to
        # pipelinesetup etc, but for now we'll create the appropriate hierarchy
        # right now
        my $schema_ss = $schema->ensure_state_hierarchy($ss);
        
        foreach my $result_node (@$nodes) {
            $schema_ss->relate_to($result_node, 'result');
        }
    }
    
    # you can specify multiple inputs here, but only 1 output
    #*** if this is an issue, this could be reimplemented to work with multiple
    #     outputs as well...
    method relate_input_to_output (ClassName|Object $self: Str|Object|ArrayRef[Str|Object] $input, Str $type, Str|Object $output, HashRef $output_file_meta?) {
        $schema ||= VRPipe::Schema->create('VRPipe');
        
        my $output_file = ref($output) ? $output : $schema->add('File', { path => $output });
        if ($output_file_meta) {
            $output_file->add_properties($output_file_meta);
        }
        
        my (@input_hash_specs, @inputs);
        foreach my $input ((ref($input) && ref($input) eq 'ARRAY') ? @$input : ($input)) {
            if (ref($input)) {
                push(@inputs, $input);
            }
            else {
                push(@input_hash_specs, { path => $input });
            }
        }
        if (@input_hash_specs) {
            push(@inputs, $schema->add('File', \@input_hash_specs));
        }
        
        if (@inputs == 1) {
            $inputs[0]->relate_to($output_file, $type);
        }
        else {
            my @spec_list;
            foreach my $input_file (@inputs) {
                push(@spec_list, { from => $input_file, to => $output_file, type => $type });
            }
            $graph ||= VRPipe::Persistent::Graph->new();
            $graph->create_mass_relationships(\@spec_list);
        }
        
        $self->result_nodes([$output_file]);
    }
}

1;
