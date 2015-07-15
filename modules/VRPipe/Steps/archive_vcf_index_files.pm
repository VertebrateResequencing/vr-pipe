
=head1 NAME

VRPipe::Steps::archive_vcf_index_files - a step

=head1 DESCRIPTION

Archive TBI or CSI index files to another disk

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2015 Genome Research Limited.

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

class VRPipe::Steps::archive_vcf_index_files extends VRPipe::Steps::archive_index_files {
    method inputs_definition {
        return {
            vcf_file       => VRPipe::StepIODefinition->create(type => 'aln', description => 'a VCF or BCF file that should be archived'),
            vcf_index_file => VRPipe::StepIODefinition->create(type => 'idx', description => 'a TBI or CSI index file that should be archived'),
        };
    }
    
    around inputs {
        my $inputs = $self->$orig();
        $inputs->{file}  = $inputs->{vcf_file};
        $inputs->{index} = $inputs->{vcf_index_file};
        return $inputs;
    }
    
    method description {
        return "Move a TBI or CSI index file to the same archive directory as it's associated VCF or BCF file";
    }

}

1;
