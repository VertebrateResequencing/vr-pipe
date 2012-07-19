
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

Copyright (c) 2011-2012 Genome Research Limited.

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
    has 'step_state' => (is      => 'rw',
                         isa     => 'VRPipe::StepState',
                         trigger => \&_rebuild);
    
    has 'data_element' => (is      => 'ro',
                           isa     => 'VRPipe::DataElement',
                           builder => '_build_data_element',
                           lazy    => 1,
                           clearer => '_clear_data_element');
    
    has 'output_root' => (is      => 'ro',
                          isa     => Dir,
                          coerce  => 1,
                          builder => '_build_output_root',
                          lazy    => 1,
                          clearer => '_clear_output_root');
    
    has 'options' => (is      => 'ro',
                      isa     => 'HashRef',
                      builder => '_resolve_options',
                      lazy    => 1,
                      clearer => '_clear_options');
    
    has 'inputs' => (is      => 'ro',
                     isa     => PersistentFileHashRef,
                     builder => '_resolve_inputs',
                     lazy    => 1,
                     clearer => '_clear_inputs');
    
    has 'outputs' => (is      => 'ro',
                      isa     => PersistentFileHashRef,
                      builder => '_build_outputs',
                      lazy    => 1,
                      clearer => '_clear_outputs');
    has 'temps' => (is      => 'ro',
                    isa     => ArrayRefOfPersistent,
                    builder => '_build_temps',
                    lazy    => 1,
                    clearer => '_clear_temps');
    
    has 'previous_step_outputs' => (is      => 'rw',
                                    isa     => PreviousStepOutput,
                                    trigger => sub { shift->_clear_inputs });
    
    has 'allow_smaller_recommended_requirements_override' => (is      => 'rw',
                                                              isa     => 'Bool',
                                                              lazy    => 1,
                                                              builder => '_build_smaller_recommended_requirements_override');
    
    # when parse is called, we'll store our dispatched refs here
    has 'dispatched' => (is      => 'ro',
                         traits  => ['Array'],
                         isa     => 'ArrayRef',
                         lazy    => 1,
                         default => sub { [] },
                         handles => { _dispatch      => 'push',
                                      num_dispatched => 'count' },
                         writer  => '_set_dispatched',
                         clearer => '_clear_dispatched');
    
    # and we'll also store all the output files the body_sub makes
    has '_output_files' => (is      => 'ro',
                            traits  => ['Hash'],
                            isa     => 'HashRef',
                            lazy    => 1,
                            default => sub { {} },
                            handles => { _remember_output_files => 'set' },
                            writer  => '_set_output_files',
                            clearer => '_clear_output_files');
    has '_temp_files' => (is      => 'ro',
                          traits  => ['Array'],
                          isa     => 'ArrayRef',
                          lazy    => 1,
                          default => sub { [] },
                          handles => { _remember_temp_file => 'push' },
                          writer  => '_set_temp_files',
                          clearer => '_clear_temp_files');
    has '_last_output_dir' => (is      => 'rw',
                               isa     => Dir,
                               lazy    => 1,
                               coerce  => 1,
                               builder => '_build_last_output_dir',
                               clearer => '_clear_last_output_dir');
    
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
        
        return dir($pipeline_root, @subdirs, $de_id, $step_state->stepmember->step_number . '_' . $self->name);
    }
    
    method _build_last_output_dir {
        return $self->output_root;
    }
    
    method _build_outputs {
        my $step_state = $self->step_state || $self->throw("Cannot get outputs without step state");
        my $outs       = $step_state->output_files;
        my $temps      = $step_state->temp_files;
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
                elsif (!$val->optional) {
                    $self->throw("the option '$key' is required");
                }
                else {
                    my $default = $val->default_value;
                    if (defined $default) {
                        $return{$key} = $default;
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
        
        my %return;
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
                
                my $num_results = @$results;
                my $max_allowed = $val->max_files;
                my $min_allowed = $val->min_files;
                if ($max_allowed == -1) {
                    $max_allowed = $num_results;
                }
                if ($min_allowed == -1) {
                    $min_allowed = $num_results;
                }
                unless ($num_results >= $min_allowed && $num_results <= $max_allowed) {
                    $self->throw("there were $num_results input file(s) for '$key' of stepstate " . $self->step_state->id . ", which does not fit the allowed range $min_allowed..$max_allowed");
                }
                
                my @vrfiles;
                my @skip_reasons;
                foreach my $result (@$results) {
                    unless (ref($result) && ref($result) eq 'VRPipe::File') {
                        $result = VRPipe::File->get(path => file($result)->absolute);
                    }
                    
                    my $wanted_type = $val->type;
                    unless ($wanted_type eq 'any') {
                        my $resolved = $result->resolve;
                        my $has_size = 0;
                        if ($resolved->s) {
                            $has_size = 1;
                            my $type = VRPipe::FileType->create($wanted_type, { file => $resolved->path });
                            unless ($type->check_type) {
                                # the type check will fail if the file doesn't
                                # really exist, so make sure
                                $resolved->update_stats_from_disc;
                                if ($resolved->s) {
                                    push(@skip_reasons, "file " . $result->path . " was not the correct type ($wanted_type)");
                                    next;
                                }
                                else {
                                    $result->update_stats_from_disc;
                                    $has_size = 0;
                                }
                            }
                        }
                        unless ($has_size) {
                            my $db_type = $result->type;
                            if ($db_type && $wanted_type ne $db_type) {
                                push(@skip_reasons, "file " . $result->path . " was not the correct type, expected type $wanted_type and got type $db_type");
                                next;
                            }
                        }
                    }
                    
                    push(@vrfiles, $result);
                }
                if (!@vrfiles && @skip_reasons) {
                    $self->throw("none of the input files had a suitable type:\n" . join("\n", @skip_reasons));
                }
                
                $return{$key} = [map { $_->e ? $_ : $_->resolve } @vrfiles];
            }
            else {
                $self->throw("invalid class " . ref($val) . " supplied for input '$key' value definition");
            }
        }
        
        return \%return;
    }
    
    method _missing (PersistentFileHashRef $hash, PersistentHashRef $defs) {
        my @missing;
        # check the files we actually need/output are as expected
        while (my ($key, $val) = each %$hash) {
            my $def     = $defs->{$key};
            my $check_s = 1;
            if ($def && $def->isa('VRPipe::StepIODefinition')) {
                $check_s = $def->check_existence;
            }
            
            foreach my $file (@$val) {
                $file->reselect_values_from_db;
                my $resolved = $file->resolve;
                if ($check_s && !$resolved->s) {
                    # double-check incase the step did not update_stats_from_disc
                    $resolved->update_stats_from_disc(retries => 3);
                    unless ($resolved->s) {
                        push(@missing, $file->path);
                    }
                    elsif ($file->id != $resolved->id) {
                        $file->update_stats_from_disc;
                    }
                }
                else {
                    my $bad = 0;
                    
                    # check the filetype is correct
                    if ($check_s) {
                        my $type = VRPipe::FileType->create($resolved->type, { file => $resolved->path });
                        unless ($type->check_type) {
                            $self->warn($resolved->path . " exists, but is the wrong type!");
                            $bad = 1;
                        }
                    }
                    
                    # check the expected metadata keys exist
                    if ($def && $def->isa('VRPipe::StepIODefinition')) {
                        my @needed = $def->required_metadata_keys;
                        if (@needed) {
                            my $meta = $file->metadata;
                            foreach my $key (@needed) {
                                unless (exists $meta->{$key}) {
                                    $self->warn($file->path . " exists, but lacks metadata key $key, needed by def " . $def->id . "!");
                                    $bad = 1;
                                }
                            }
                        }
                    }
                    
                    if ($bad) {
                        push(@missing, $file->path);
                    }
                }
            }
        }
        
        return @missing;
    }
    
    method missing_input_files {
        return $self->_missing($self->inputs, $self->inputs_definition);
    }
    
    method missing_output_files {
        my $hash = $self->outputs;
        my $defs = $self->outputs_definition;
        
        # check that we don't have any outputs defined in the definition that
        # no files were made for
        while (my ($key, $val) = each %$defs) {
            next if exists $hash->{$key};
            next if $val->min_files == 0;
            $self->throw("'$key' was defined as an output, yet no output file was made with that output_key");
        }
        
        return $self->_missing($hash, $defs);
    }
    
    method _run_coderef (Str $method_name) {
        my $ref = $self->$method_name();
        return &$ref($self);
    }
    
    method output_file (Str :$output_key?, File|Str :$basename!, FileType :$type!, Dir|Str :$output_dir?, Dir|Str :$sub_dir?, HashRef :$metadata?, Bool :$temporary = 0) {
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
        $self->make_path($output_dir); #*** repeated, potentially unecessary filesystem access...
        $self->_last_output_dir($output_dir);
        
        my $vrfile = VRPipe::File->create(path => file($output_dir, $basename), type => $type);
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
    
    method parse {
        # if we have missing input files, check to see if some other step
        # created them, and start those steps over in the hopes the files will
        # be recreated; otherwise throw
        my @missing = $self->missing_input_files;
        if (@missing) {
            my $with_recourse = 0;
            my %states_to_restart;
            foreach my $path (@missing) {
                my $file = VRPipe::File->get(path => $path);
                my $resolved = $file->resolve;
                next if $resolved->s; # there's no recourse if the file was actually just missing some metadata, not physically missing
                my $count = 0;
                my $state;
                foreach my $sof (VRPipe::StepOutputFile->search({ file => $file->id }, { prefetch => 'stepstate' })) {
                    $count++;
                    $state = $sof->stepstate;
                }
                
                if ($count == 1) {
                    $with_recourse++;
                    push(@{ $states_to_restart{ $state->id } }, $resolved->path);
                }
            }
            
            if ($with_recourse == @missing) {
                while (my ($state_id, $files) = each %states_to_restart) {
                    #*** would have to make inactive pipeline active as well...
                    #    how to make sure Manager parses this step soon?...
                    my $state = VRPipe::StepState->get(id => $state_id);
                    if ($state->complete) {
                        $self->debug("to regenerate needed input files (@$files) for stepstate " . $self->step_state->id . ", will start stepstate $state_id over again");
                        $state->start_over;
                    }
                }
                return 0;
            }
            else {
                my $step_state = $self->step_state;
                $self->throw("Required input files are missing for step " . $self->name . " (for data element " . $self->data_element->id . ", pipelinesetup " . $step_state->pipelinesetup->id . ", stepstate " . $step_state->id . "): [@missing]");
            }
        }
        
        # $self is a VRPipe::Step, even when body_sub was defined in a
        # VRPipe::Steps::subclass; to regain full benefits of inheritence in
        # those body_sub subs, we'll load the real module and use that instead
        # of $self
        my $non_persistent = $self->_from_non_persistent;
        if ($non_persistent) {
            #*** this is pretty ugly - is there a better way?
            $non_persistent->step_state($self->step_state);
            $non_persistent->previous_step_outputs($self->previous_step_outputs);
            
            $non_persistent->_run_coderef('body_sub');
            
            $self->_set_output_files($non_persistent->_output_files);
            $self->_set_temp_files($non_persistent->_temp_files);
            $self->_last_output_dir($non_persistent->_last_output_dir);
            $self->_set_dispatched($non_persistent->dispatched);
        }
        else {
            $self->_run_coderef('body_sub');
        }
        
        # store output and temp files on the StepState
        my $output_files = $self->_output_files;
        if (keys %$output_files) {
            my $step_state = $self->step_state;
            $step_state->output_files($output_files);
            $step_state->update;
        }
        my $temp_files = $self->_temp_files;
        if (@$temp_files) {
            my $step_state = $self->step_state;
            $step_state->temp_files($temp_files);
            $step_state->update;
        }
        
        # return true if we've finished the step
        my $dispatched = $self->num_dispatched;
        if ($dispatched) {
            return 0;
        }
        else {
            return $self->post_process;
        }
    }
    
    method post_process {
        my $ok         = $self->_run_coderef('post_process_sub');
        my $stepstate  = $self->step_state;
        my $debug_desc = "step " . $self->name . " failed for data element " . $self->data_element->id . " and pipelinesetup " . $self->step_state->pipelinesetup->id . " (stepstate " . $stepstate->id . ")";
        
        if ($ok) {
            my @missing = $self->missing_output_files;
            $stepstate->unlink_temp_files;
            if (@missing) {
                $self->warn("Some output files are missing (@missing) for $debug_desc");
                $stepstate->start_over;
                return 0;
            }
            else {
                return 1;
            }
        }
        else {
            $stepstate->unlink_temp_files;
            $self->warn("The post-processing part of $debug_desc");
            $stepstate->start_over;
            return 0;
        }
    }
    
    method new_requirements (Int :$memory!, Int :$time!, Int :$cpus?, Int :$tmp_space?, Int :$local_space?, HashRef :$custom?) {
        # get the current mean+2sd memory and time of past runs of this step
        my $ssu      = VRPipe::StepStatsUtil->new(step => $self->isa('VRPipe::Step') ? $self : VRPipe::Step->get(name => $self->name));
        my $rec_mem  = $ssu->recommended_memory;
        my $rec_time = $ssu->recommended_time;
        
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
        if (defined $options->{memory_override} && $options->{memory_override} > $memory) {
            $memory = $options->{memory_override};
        }
        if (defined $options->{time_override} && $options->{time_override} > $time) {
            $time = $options->{time_override};
        }
        #*** and the other resources?...
        
        return
          VRPipe::Requirements->create(memory => $memory,
                                       time   => $time,
                                       $cpus        ? (cpus        => $cpus)        : (),
                                       $tmp_space   ? (tmp_space   => $tmp_space)   : (),
                                       $local_space ? (local_space => $local_space) : (),
                                       $custom      ? (custom      => $custom)      : ());
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
                    $new_lib{ file($inc)->absolute } = 1;
                }
            }
            my @new_lib = keys %new_lib;
            if (@new_lib) {
                $use_lib = "use lib(qw(@new_lib)); ";
            }
            
            $deploy = 'VRPipe::Persistent::SchemaBase->database_deployment(q[testing]); ';
        }
        
        my $cmd = qq[perl -MVRPipe::Persistent::Schema -e "$use_lib$deploy$code"];
        $self->dispatch([$cmd, $req, $extra_args]);
    }
    
    method dispatch_md5sum (VRPipe::File $vrfile, Maybe[Str] $expected_md5) {
        my $path = $vrfile->path;
        my $req = $self->new_requirements(memory => 500, time => 1);
        
        if ($expected_md5) {
            return $self->dispatch_vrpipecode(qq[use Digest::MD5; open(FILE, q[$path]) or die q[Could not open file $path]; binmode(FILE); if (Digest::MD5->new->addfile(*FILE)->hexdigest eq q[$expected_md5]) { VRPipe::File->get(path => q[$path], md5 => q[$expected_md5]); } else { die q[md5sum of $path does not match expected value] }], $req);
        }
        else {
            return $self->dispatch_vrpipecode(qq[use Digest::MD5; open(FILE, q[$path]) or die q[Could not open file $path]; binmode(FILE); VRPipe::File->get(path => q[$path], md5 => Digest::MD5->new->addfile(*FILE)->hexdigest);], $req);
        }
    }
    
    # using this lets you run a bunch of perl code wrapping a command line exe,
    # yet still keep the command line visible so you know the most important
    # bit of what was run
    method dispatch_wrapped_cmd (ClassName $class, Str $method, ArrayRef $dispatch_args) {
        my ($cmd, $req, $extra_args) = @$dispatch_args;
        
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
            foreach my $key (keys %$file_meta) {
                $meta{$key}->{ $$file_meta{$key} } += 1;
            }
        }
        # Only keep metadata common to all files
        my $common_meta = {};
        foreach my $key (keys %meta) {
            my @vals = keys %{ $meta{$key} };
            next unless (@vals == 1 && $meta{$key}->{ $vals[0] } == @$files);
            $common_meta->{$key} = $vals[0];
        }
        return $common_meta;
    }
}

1;
