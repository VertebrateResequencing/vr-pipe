
=head1 NAME

VRPipe::Steps::picard - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>. Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011,2012,2015 Genome Research Limited.

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

class VRPipe::Steps::picard extends VRPipe::Steps::java {
    has 'picard_path' => (
        is  => 'rw',
        isa => 'Str'
    );
    
    has '+memory_multiplier' => (default => 0.7);
    
    around _build_standard_options {
        return [@{ $self->$orig }, 'picard_path'];
    }
    
    our %PICARD_VERSIONS;
    has 'picard_version' => (
        is      => 'ro',
        isa     => 'Str',
        lazy    => 1,
        builder => 'determine_picard_version'
    );
    
    method determine_picard_version (ClassName|Object $self:) {
        my $picard_path = $self->picard_path;
        unless (defined $PICARD_VERSIONS{$picard_path}) {
            my $version = 0;
            
            if (-d $picard_path) {
                opendir(my $dh, $picard_path) || $self->throw("Could not open picard directory $picard_path");
                foreach (readdir $dh) {
                    if (/^picard-([\d\.]+)\.jar/) {
                        $version = $1;
                        last;
                    }
                }
                closedir($dh);
            }
            else {
                if ($picard_path =~ /picard-(\d+\.\d+)\.jar/) {
                    $version = $1;
                }
                else {
                    $version = `java -jar $picard_path ViewSam --version 2>&1`;
                    ($version) = $version =~ /^(\d+\.\d+)/;
                }
            }
            
            $PICARD_VERSIONS{$picard_path} = $version;
        }
        return $PICARD_VERSIONS{$picard_path};
    }
    
    method jar (ClassName|Object $self: Str $basename!) {
        my $path = $self->picard_path;
        if (-d $path) {
            return file($path, $basename);
        }
        else {
            $basename =~ s/\.jar$//;
            return "$path $basename";
        }
    }
    
    around options_definition {
        return { %{ $self->$orig }, picard_path => VRPipe::StepOption->create(description => 'path to the directory containing Picard jar files (old Picard), or the path to the picard.jar file (new Picard)', $ENV{PICARD} ? (default_value => $ENV{PICARD}) : ()) };
    }
    
    method inputs_definition {
        return {};
    }
    
    method body_sub {
        return sub { return 1; };
    }
    
    method outputs_definition {
        return {};
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Generic step for using Picard";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
