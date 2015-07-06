
=head1 NAME

VRPipe::Steps::gatk - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::gatk extends VRPipe::Steps::java {
    has 'gatk_path' => (
        is     => 'rw',
        isa    => Dir,
        coerce => 1
    );
    
    has 'gatk_key' => (
        is  => 'rw',
        isa => 'Str',
    );
    
    around _build_standard_options {
        return [@{ $self->$orig }, 'gatk_path', 'gatk_key'];
    }
    
    our %GATK_VERSIONS;
    has 'gatk_version' => (
        is      => 'ro',
        isa     => 'Str',
        lazy    => 1,
        builder => 'determine_gatk_version'
    );
    
    method jar (ClassName|Object $self:) {
        return file($self->gatk_path, 'GenomeAnalysisTK.jar');
    }
    
    method determine_gatk_version (ClassName|Object $self:) {
        my $gatk_jar = $self->jar->stringify;
        unless (defined $GATK_VERSIONS{$gatk_jar}) {
            my $jvm_args = $self->jvm_args(50);
            my $java_exe = $self->java_exe;
            $GATK_VERSIONS{$gatk_jar} = VRPipe::StepCmdSummary->determine_version(qq[$java_exe $jvm_args -jar $gatk_jar -h], 'v([\d\.\-]+[a-z\d]*)');
        }
        return $GATK_VERSIONS{$gatk_jar};
    }
    
    method gatk_prefix (ClassName|Object $self: Int $memory, Str|Dir $dir?) {
        my $jvm_args = $self->jvm_args($memory, $dir ? $dir : ());
        my $gatk_key = $self->gatk_key;
        my $prefix   = $self->java_exe . qq[ $jvm_args -jar ] . $self->jar;
        
        # we want to always run with --phone_home NO_ET, but can only do that
        # in v1.5+ if we have a gatk_key
        my $version = $self->gatk_version;
        $version =~ s/^(\d+\.\d+).+/$1/;
        if ($version < 1.5) {
            $prefix .= ' --phone_home NO_ET';
        }
        elsif ($gatk_key) {
            $prefix .= qq[ --phone_home NO_ET --gatk_key $gatk_key];
        }
        
        return $prefix;
    }
    
    around options_definition {
        return {
            %{ $self->$orig },
            reference_fasta => VRPipe::StepOption->create(description => 'absolute path to genome reference file used to do the mapping'),
            gatk_path       => VRPipe::StepOption->create(description => 'path to directory containing GATK jar files', $ENV{GATK} ? (default_value => $ENV{GATK}) : ()),
            gatk_key        => VRPipe::StepOption->create(description => 'path to key file to avoid GATK phoning home', optional => 1, $ENV{GATK_KEY} ? (default_value => $ENV{GATK_KEY}) : ()),
        };
    }
    
    around options {
        my $options         = $self->$orig;
        my $reference_fasta = file($options->{reference_fasta});
        $self->throw("reference_fasta must be an absolute path") unless $reference_fasta->is_absolute;
        return $options;
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
        return "Generic step for using the GenomeAnalysisToolkit (GATK)";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}

1;
