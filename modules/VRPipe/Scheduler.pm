
=head1 NAME

VRPipe::Scheduler - a generic interface to job scheduling systems

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

In order to manage the execution of command lines across the many nodes of a
compute cluster, some job scheduling system should be in place. Examples
include LSF and Grid Engine. B<VRPipe> submits the work it wants done to a
scheduler; this Scheduler provides a single consistent interface to all the
possible schedulers.

Currently only LSF is supported. Support for other shedulers can be added by
creating a C<VRPipe::Schedulers::[name]> class that implements the
C<VRPipe::SchedulerMethodsRole> role.

For doing work on the local machine when there isn't a cluster available (or
for testing purposes), there is also a 'Local' scheduler supplied with
B<VRPipe>.

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

class VRPipe::Scheduler extends VRPipe::Persistent {
    use VRPipe::Config;
    my $vrp_config = VRPipe::Config->new();
    use VRPipe::Persistent::SchemaBase;
    use VRPipe::SchedulerMethodsFactory;
    use VRPipe::Interface::CmdLine;
    use DateTime;
    
    has 'type' => (
        is                   => 'rw',
        isa                  => Varchar [64],
        builder              => 'default_type',
        traits               => ['VRPipe::Persistent::Attributes'],
        is_key               => 1,
        allow_key_to_default => 1
    );
    
    has 'output_root' => (
        is      => 'ro',
        isa     => Dir,
        coerce  => 1,
        builder => 'default_output_root',
        lazy    => 1
    );
    
    method default_type (ClassName|Object $self:) {
        my $method_name = VRPipe::Persistent::SchemaBase->database_deployment . '_scheduler';
        my $type        = $vrp_config->$method_name();
        return lc($type);
    }
    
    method default_output_root (ClassName|Object $self:) {
        my $method_name = VRPipe::Persistent::SchemaBase->database_deployment . '_logging_directory';
        my $root        = $vrp_config->$method_name();
        return "$root";       # stringify what could be a VRPipe::Base::Configuration::Env
    }
    
    # VRPipe::Schedulers::[type] classes will provide scheduler-specific
    # methods
    has 'scheduler_instance' => (
        is      => 'ro',
        isa     => 'Object',
        builder => '_instantiate_method_class',
        lazy    => 1,
        handles => 'VRPipe::SchedulerMethodsRole'
    );
    
    method _instantiate_method_class (ClassName|Object $self:) {
        return VRPipe::SchedulerMethodsFactory->create(lc($self->type), {});
    }
    
    method start_scheduler {
        my $cmd = $self->start_command;
        system("$cmd > /dev/null 2> /dev/null");
    }
    
    method stop_scheduler {
        my $cmd = $self->stop_command;
        system("$cmd > /dev/null 2> /dev/null");
    }
    
    method ensure_running (Str :$cmd!, VRPipe::Requirements :$requirements!, PositiveInt :$count = 1, Str :$cwd?) {
        $cmd =~ s/^\S+perl/perl/;   # different nodes may have different perls installed at different locations
        
        # get the details of everything already in the scheduler for this cmd,
        # removing from the queue anything not currently running when we're over
        # the desired count
        my ($scheduled_cmd_count, $running_sidaids) = $self->command_status(cmd => $cmd, max => $count);
        return @$running_sidaids if $scheduled_cmd_count >= $count;
        
        # submit $still_needed new jobs that all run $cmd
        my $still_needed       = $count - $scheduled_cmd_count;
        my $scheduler_cmd_line = join(
            ' ',
            $self->submit_command,
            $self->submit_args(
                requirements => $requirements,
                stdo_file    => '/dev/null',
                stde_file    => '/dev/null',
                $cwd ? (cwd => $cwd) : (),
                cmd   => $cmd,
                count => $still_needed
            )
        );
        system($scheduler_cmd_line);
        
        return @$running_sidaids;
    }
    
    method output_dir (PersistentObject $for) {
        my $root_dir = $self->output_root;
        
        my $hashing_string = ref($for) . '::' . $for->id;
        my @subdirs        = $self->hashed_dirs($hashing_string);
        
        $hashing_string =~ s/\:/_/g;
        my $output_dir = dir($root_dir, @subdirs, $hashing_string);
        
        $output_dir->mkpath;
        return $output_dir;
    }
    
    __PACKAGE__->make_persistent();
}

1;
