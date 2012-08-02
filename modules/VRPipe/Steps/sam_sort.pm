
=head1 NAME

VRPipe::Steps::sam_sort - a step

=head1 DESCRIPTION

Sort a sam file.

Currently this step sorts a sam file by co-oridinate.

*** more documentation to come





=head1 AUTHOR

NJWalker <nw11@sanger.ac.uk>.

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

class VRPipe::Steps::sam_sort with VRPipe::StepRole {
    use File::Basename;
    use Data::Dumper;
    
    method options_definition {
        return {};
    }
    
    method inputs_definition {
        return { sam_file => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'A sam file') };
    }
    
    method body_sub {
        return sub {
            my $self = shift;
            #my $options = $self->options;
            my ($infile) = @{ $self->inputs->{sam_file} };
            my $name = $infile->basename;
            my $outfile = $self->output_file(output_key => 'sorted_sam',
                                             basename   => $name . '.sort' . '.sam',
                                             type       => 'txt',
                                             metadata   => $infile->metadata);
            my $infile_path  = $infile->path;
            my $outfile_path = $outfile->path;
            my $req          = $self->new_requirements(memory => 1500, time => 1);
            my $cmd          = qq{ grep '^\@' $infile_path > $outfile_path; grep -v '^\@' $infile_path | sort -k 3,3V -k 4,4n >> $outfile_path; };
            $self->dispatch([qq[$cmd], $req, { output_files => [$outfile] }]);
          }
    }
    
    method outputs_definition {
        return { sorted_sam => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'Sorted sam file by co-ordinates.'), };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Sort a sam file using - grep '^\@' file_to_sort.sam > sorted_file.sam; grep-v '^\@' file_to_sort.sam | sort -k 3,3 -k4,4n >> sorted_file.sam";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }

}
