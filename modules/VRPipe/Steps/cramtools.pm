=head1 NAME

VRPipe::Steps::cramtools - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::cramtools extends VRPipe::Steps::java {
    has 'cramtools_path' => (is => 'rw',
                             isa => Dir,
                             coerce => 1);
    
    around _build_standard_options {
        return [@{$self->$orig}, 'cramtools_path'];
    }
    
    our %CRAMTOOLS_VERSIONS;
    has 'cramtools_version' => (is => 'ro',
                                isa => 'Str',
                                lazy => 1,
                                builder => 'determine_cramtools_version');
    
    method jar (ClassName|Object $self:) {
        return file($self->cramtools_path, 'cramtools.jar');
    }
    
    method determine_cramtools_version (ClassName|Object $self:) {
        my $cramtools_jar = $self->jar->stringify;
        unless (defined $CRAMTOOLS_VERSIONS{$cramtools_jar}) {
            my $jvm_args = $self->jvm_args(50);
            my $java_exe = $self->java_exe;
            $CRAMTOOLS_VERSIONS{$cramtools_jar} = VRPipe::StepCmdSummary->determine_version(qq[$java_exe $jvm_args -jar $cramtools_jar -h], 'v([\d\.\-]+[a-z\d]*)');
        }
        return $CRAMTOOLS_VERSIONS{$cramtools_jar};
    }
    
    around options_definition {
        return { %{$self->$orig},
                 reference_fasta => VRPipe::StepOption->get(description => 'absolute path to genome reference file used to do the mapping'),
                 cramtools_path => VRPipe::StepOption->get(description => 'path to cramtools jar file', optional => 1, default_value => "$ENV{CRAMTOOLS}"),
                };
    }
    method inputs_definition {
        return { };
    }
    method body_sub {
        return sub { return 1; };
    }
    method outputs_definition {
        return { };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Generic step for using the cramtools";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
}

1;
