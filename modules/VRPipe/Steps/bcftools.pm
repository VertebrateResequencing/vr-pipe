
=head1 NAME

VRPipe::Steps::bcftools - a step

=head1 DESCRIPTION

Generic non-functional step for bcftools-using steps to inherit from

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2014 Genome Research Limited.

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

class VRPipe::Steps::bcftools with VRPipe::StepRole {
    has 'bcftools_version_string' => (
        is      => 'ro',
        isa     => 'Str',
        lazy    => 1,
        builder => '_build_bcftools_version_string'
    );
    
    has 'bcftools_v1' => (
        is      => 'ro',
        isa     => 'Bool',
        lazy    => 1,
        builder => '_build_bcftools_v1'
    );
    
    method _build_bcftools_version_string {
        return VRPipe::StepCmdSummary->determine_version($self->options->{bcftools_exe}, '^Version: (.+)$');
    }
    
    method _build_bcftools_v1 {
        # the dev release was v 0.[01].*, and then the move to using htslib
        # resulted in usage changes, starting with dev release 0.2+, which is
        # planned to become version 1.0; so we say if this is v1 usage or not
        my $version_string = $self->bcftools_version_string;
        my ($major, $minor) = $version_string =~ /(\d+)\.(\d+)/; #0.1.19 vs 0.2.0
        return $major > 0 || $minor > 1;
    }
    
    method options_definition {
        return {
            bcftools_exe => VRPipe::StepOption->create(
                description   => 'path to bcftools executable',
                optional      => 1,
                default_value => 'bcftools'
            ),
        };
    }
    
    method inputs_definition {
        return {};
    }
    
    method body_sub {
        return sub { };
    }
    
    method outputs_definition {
        return {};
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Generic step for steps using bcftools - don't use directly";
    }
    
    method max_simultaneous {
        return 0;                                                # meaning unlimited
    }
    
    method _bcftools_calling_command {
        return $self->bcftools_v1 ? 'call' : 'view';
    }
    
    method _bcftools_calling_defaults {
        return $self->bcftools_v1 ? '-m' : '-p 0.99 -vcgN';
    }
    
    method _bcftools_calling_for_old_style_snps_only {
        return $self->bcftools_v1 ? '-cvV indels' : '-cgvI';
    }
    
    method _bcftools_compressed_vcf_output (Maybe[Str] $filter?) {
        if ($filter) {
            return "| $filter | bgzip -c";
        }
        return $self->bcftools_v1 ? '-O z' : '| bgzip -c ';
    }
    
    method _bcftools_samples_option {
        return $self->bcftools_v1 ? '-S' : '-s';
    }
    
    method _bcftools_site_files_option (File $path) {
        my $arg = $self->bcftools_v1 ? 'T' : 'l';
        return " -$arg $path";
    }
}

1;
