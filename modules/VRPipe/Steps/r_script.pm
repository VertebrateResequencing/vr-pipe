
=head1 NAME

VRPipe::Steps::r_script - a step

=head1 DESCRIPTION

Generic step for steps using Rscript command, providing params for the Rscript
executable and R libaries pathset, and generating a standard Rscript command
prefix

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

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

class VRPipe::Steps::r_script with VRPipe::StepRole {
    use POSIX qw(ceil);
    
    has 'rscript_cmd' => (
        is  => 'rw',
        isa => 'Str'
    );
    
    has 'r_libs' => (
        is  => 'rw',
        isa => 'Str'
    );
    
    has 'standard_options' => (
        is      => 'ro',
        isa     => 'ArrayRef',
        lazy    => 1,
        builder => '_build_standard_options'
    );
    
    method _build_standard_options {
        return ['rscript_cmd', 'r_libs'];
    }
    
    method rscript_cmd_prefix {
        my $options = $self->options;
        my $r_libs  = $self->r_libs;
        my $return;
        if ($r_libs) {
            $return = "export R_LIBS=" . $self->r_libs . ";";
        }
        $return .= $self->rscript_cmd;
        return $return;
    }
    
    method handle_standard_options (HashRef $options) {
        foreach my $method (@{ $self->standard_options }) {
            next unless defined $options->{$method};
            $self->$method($options->{$method});
        }
    }
    
    method options_definition {
        return {
            rscript_cmd => VRPipe::StepOption->create(description => 'path to your Rscript executable and default arguments', optional => 1, default_value => 'Rscript --vanilla'),
            r_libs => VRPipe::StepOption->create(description => 'R libraries path set', optional => 1, $ENV{R_LIBS} ? (default_value => $ENV{R_LIBS}) : ()),
        };
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
        return "Generic step for steps using R in batch mode";
    }
    
    method max_simultaneous {
        return 0;                 # meaning unlimited
    }
}

1;
