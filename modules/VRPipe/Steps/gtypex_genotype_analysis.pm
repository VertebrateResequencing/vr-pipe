use VRPipe::Base;

class VRPipe::Steps::gtypex_genotype_analysis with VRPipe::StepRole {
    method options_definition {
        return { gtype_confidence => VRPipe::StepOption->get(description => 'confidence level to be used as a cutoff to determine whether genotyping has passed or not',
                                                             optional => 1,
                                                             default_value => 1.05) };
    }
    method inputs_definition {
        return { gtypex_files => VRPipe::StepIODefinition->get(type => 'txt', max_files => -1,
        		                                       description => 'gtypex files containing likelihood scores for genotyping',
                                                               metadata => {sample => 'name of expected sample', source_bam => 'input bam path'}}) };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $min_ratio = $options->{gtype_confidence};
            my $req = $self->new_requirements(memory => 3900, time => 1);
	    
            foreach my $gtypex (@{$self->inputs->{gtypex_files}}) {
		my $meta = $gtypex->metadata;
		my $expected = $meta->{sample};
		my $source_bam = $meta->{source_bam};
            	my $bam_basename = file($source_bam)->basename;
                my $gtypex_path = $gtypex->path;          
                my $gtype_file = $self->output_file(output_key => 'gtype_files',
                                                    basename => $bam_basename.'.gtype',
                                                    type => 'txt',
                                                    metadata => {sample => $expected, source_bam => $source_bam});
                my $gtype_path = $gtype_file->path;
                my $cmd = "use VRPipe::Steps::gtypex_genotype_analysis; VRPipe::Steps::gtypex_genotype_analysis->analyse_gtypex_output(gtypex => q[$gtypex_path], gtype => q[$gtype_path], confidence => q[$min_ratio]);";
                $self->dispatch_vrpipecode($cmd, $req, {output_files => [$gtype_file]});
            }
        };
    }
    method outputs_definition {
        return { gtype_files => VRPipe::StepIODefinition->get(type => 'txt',
							      max_files => -1,
							      description => 'file of genotype status report for the sample',
							      metadata => {sample => 'sample name'}) };
    }
    method post_process_sub {
        return sub { return 1; }
    }
    method description {
        return "The genotype status is determined by analysing the gtypex file of likelihood scores of sample genotypes and comparing the score ratio (over the 'next best' genotype) with the confidence level cutoff provided";
    }
    method max_simultaneous {
        return 0;
    }
    
    method analyse_gtypex_output (ClassName|Object $self: Str|File :$gtypex, Str|File :$gtype, Num :$confidence) {
        my $gtypex_file = VRPipe::File->get(path => $gtypex);
        my $meta = $gtypex_file->metadata;
        my $expected = $meta->{sample};
        my $gtype_file = VRPipe::File->get(path => $gtype);
        my $fh = $gtypex_file->openr;
        my $ofh = $gtype_file->openw;
        $gtype_file->disconnect;
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
        my $gt_status;
        if ( $expected_gtype2 ) {
            $gt_status = "status=confirmed expected=$expected found=$gtype2 ratio=$ratio\n"; 
        }
        elsif ( $ratio<$confidence ) {
            if ( $expected ) {
                $gt_status = "status=unconfirmed expected=$expected found=$gtype1 ratio=$ratio\n"; 
            }
            else {
               	$gt_status = "status=unknown expected=none found=$gtype1 ratio=$ratio\n";
            }	
        }
        elsif ( !$expected ) {
            $gt_status = "status=candidate expected=none found=$gtype1 ratio=$ratio\n"; 
        }
        elsif ( $expected eq $gtype1 ) {
            $gt_status = "status=confirmed expected=$expected found=$gtype1 ratio=$ratio\n"; 
        }
        else {
            $gt_status = "status=wrong expected=$expected found=$gtype1 ratio=$ratio\n";
        }
        print $ofh $gt_status;
        close $ofh;
    }
}

1;
