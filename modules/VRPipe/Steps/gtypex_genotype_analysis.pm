=head1 NAME

VRPipe::Steps::gtypex_genotype_analysis - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

John Maslen <jm23@sanger.ac.uk>.

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

class VRPipe::Steps::gtypex_genotype_analysis with VRPipe::StepRole {
    method options_definition {
        return { gtype_confidence => VRPipe::StepOption->create(description => 'confidence level to be used as a cutoff to determine whether genotyping has passed or not',
                                                             optional => 1,
                                                             default_value => 1.05) };
    }
    method inputs_definition {
        return { gtypex_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1,
        		                                       description => 'gtypex files containing likelihood scores for genotyping',
                                                               metadata => {expected_sample => 'name of expected sample', source_bam => 'input bam path'}) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $min_ratio = $options->{gtype_confidence};
            my $req = $self->new_requirements(memory => 3900, time => 1);
	    
            foreach my $gtypex (@{$self->inputs->{gtypex_files}}) {
		my $source_bam = $gtypex->metadata->{source_bam};
                my $gtypex_path = $gtypex->path;
                my $cmd = "use VRPipe::Steps::gtypex_genotype_analysis; VRPipe::Steps::gtypex_genotype_analysis->analyse_gtypex_output(gtypex => q[$gtypex_path], confidence => q[$min_ratio], source_bam => q[$source_bam]);";
                $self->dispatch_vrpipecode($cmd, $req);
            }
        };
    }
    method outputs_definition {
        return { };
    }
    method post_process_sub {
        return sub { return 1; }
    }
    method description {
        return "The genotype status is determined by analysing the gtypex file of likelihood scores of sample genotypes and comparing the score ratio (over the 'next best' genotype) with the confidence level cutoff provided; results are stored as metadata (gtype_analysis) on the bam file";
    }
    method max_simultaneous {
        return 0;
    }
    
    method analyse_gtypex_output (ClassName|Object $self: Str|File :$gtypex, Str|File :$source_bam, Num :$confidence) {
        my $gtypex_file = VRPipe::File->get(path => $gtypex);
        my $meta = $gtypex_file->metadata;
        my $expected = $meta->{expected_sample};
        my $bam_file = VRPipe::File->get(path => $source_bam);
        my $fh = $gtypex_file->openr;
        $bam_file->disconnect;
        my $has_data = 0;
    	my ($hit1,$hit2,$gtype1,$lhood1,$gtype2,$lhood2);
	my $entrp = <$fh>;
    	while (my $line = <$fh>) {
    	    if ( $has_data && defined($hit1) && defined($hit2) ) { last; }
            if ( !($line =~ m{sample\s+(?:.*/)?(\S+?)(?:\.snp)?\s+likelihood \d+ over \d+ sites, score (\S+),} ) ) { 
                $self->throw("Could not parse $gtypex: $hit1") 
            }
	    if ( $expected && $1 eq $expected ) { $has_data = 1; } 
            if ( !defined $hit1 ) {
                $hit1   = $line;
                $gtype1 = $1;
                $lhood1 = $2;
            }
            elsif ( !defined $hit2 ) {
                $hit2   = $line; 
                $gtype2 = $1;
                $lhood2 = $2;
            }
        }
        close $fh;
        
        my $expected_gtype2 = ($expected eq $gtype2 && $lhood1 == $lhood2) ? 1 : 0;
        if ( $expected && !$has_data ) { $expected = 0; }
        my $ratio = $lhood1!=0 ? $lhood2/$lhood1 : $lhood2/1e-6;
	$ratio = sprintf("%0.3f", $ratio);
        my $gt_status;
        if ( $expected_gtype2 ) {
            $gt_status = "status=confirmed expected=$expected found=$gtype2 ratio=$ratio"; 
        }
        elsif ( $ratio<$confidence ) {
            if ( $expected ) {
                $gt_status = "status=unconfirmed expected=$expected found=$gtype1 ratio=$ratio"; 
            }
            else {
               	$gt_status = "status=unknown expected=none found=$gtype1 ratio=$ratio";
            }	
        }
        elsif ( !$expected ) {
            $gt_status = "status=candidate expected=none found=$gtype1 ratio=$ratio"; 
        }
        elsif ( $expected eq $gtype1 ) {
            $gt_status = "status=confirmed expected=$expected found=$gtype1 ratio=$ratio"; 
        }
        else {
            $gt_status = "status=wrong expected=$expected found=$gtype1 ratio=$ratio";
        }
        
	my $new_meta = {gtype_analysis => $gt_status};
	$bam_file->add_metadata($new_meta, replace_data => 1);
	$gtypex_file->add_metadata($new_meta, replace_data => 1);
    }
}

1;
