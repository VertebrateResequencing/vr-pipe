
=head1 NAME

VRPipe::Steps::bam_metadata_with_sex - a step

=head1 DESCRIPTION

Extends bam_metadata to add Sample Gender meta-data to bam files

=head1 AUTHOR

Chris Joyce <cj5@sanger.ac.uk>.

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

class VRPipe::Steps::bam_metadata_with_sex extends VRPipe::Steps::bam_metadata {
    around options_definition {
        my $options = $self->$orig;
        return { %{ $self->$orig },
                 sample_sex_file => VRPipe::StepOption->create(description => 'File listing the sex (eg M or F) of samples, assumed_sex is used if no file provided',                      optional => 1),
                 assumed_sex     => VRPipe::StepOption->create(description => 'If sex is not present for a sample in the sample sex file (or no file provided), then this sex is assumed', optional => 1), };
    }
    
    method post_process_sub {
        return sub {
            my $self = shift;
            
            my $options         = $self->options;
            my $assumed_sex     = $options->{assumed_sex};
            my $sample_sex_file = $options->{sample_sex_file};
            
            my %sample_sex;
            
            if ($sample_sex_file) {
                my $sample_file = Path::Class::File->new($sample_sex_file);
                $self->throw("sample_sex_file must be an absolute path") unless $sample_file->is_absolute;
                
                my $fh = $sample_file->openr;
                while (<$fh>) {
                    chomp;
                    my ($sample, $sex) = split;
                    $sample_sex{$sample} = $sex;
                }
                $fh->close;
            }
            
            foreach my $bamfile (@{ $self->inputs->{bam_files} }) {
                my $meta = $bamfile->metadata;
                unless (defined $meta->{sex}) {
                    my $bam_sample = $meta->{sample};
                    
                    my $bam_sex;
                    if ($sample_sex{$bam_sample}) {
                        $bam_sex = $sample_sex{$bam_sample};
                    }
                    else {
                        $self->throw("Cannot assign a sex for sample $bam_sample") unless $assumed_sex;
                        $bam_sex = $assumed_sex;
                    }
                    $bamfile->add_metadata({ sex => $bam_sex });
                }
            }
            return 1;
        };
    }
    
    around _build_meta_to_check {
        my $meta_to_check = $self->$orig;
        return [@{$meta_to_check}, 'sample'];
    }
    
    method description {
        return "Extends bam_metadata to Add Sample Gender meta-data to bam files";
    }

}

1;
