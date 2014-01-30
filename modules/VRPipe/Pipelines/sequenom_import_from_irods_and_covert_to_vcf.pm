
=head1 NAME

VRPipe::Pipelines::sequenom_import_from_irods_and_covert_to_vcf - a pipeline

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2013 Genome Research Limited.

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

class VRPipe::Pipelines::sequenom_import_from_irods_and_covert_to_vcf with VRPipe::PipelineRole {
    method name {
        return 'sequenom_import_from_irods_and_covert_to_vcf';
    }
    
    method description {
        return 'Import sequenom CSV files from iRODS, then convert the CSV files to VCFs, adding sequenom_gender metadata.';
    }
    
    method step_names {
        (
            'irods_get_files_by_basename', #1
            'sequenom_csv_to_vcf',         #2
        );
    }
    
    method adaptor_definitions {
        (
            { from_step => 0, to_step => 1, to_key   => 'basenames' },
            { from_step => 1, to_step => 2, from_key => 'local_files', to_key => 'csv_files' },
        );
    }
    
    method behaviour_definitions {
        (
            { after_step => 2, behaviour => 'delete_outputs', act_on_steps => [1], regulated_by => 'cleanup', default_regulation => 1 },
        );
    }

}

1;
