use VRPipe::Base;

class VRPipe::Steps::genotype_analysis with VRPipe::StepRole {

	method options_definition {
        return { gtype_confidence => VRPipe::StepOption->get(description => 'confidence level to be used as a cutoff to determine whether genotyping has passed or not',
                                                            optional => 1,
                                                            default_value => 1.05),
               };
    }
    
    method inputs_definition {
        return { 
        		 gtypex_files => VRPipe::StepIODefinition->get(type => 'txt', max_files => -1,
        		                                               description => 'gtypex files containing likelihood scores for genotyping',
                                                               metadata => {sample => 'sample name'}
                                                               )
        		 };
	}
	
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $min_ratio = $options->{gtype_confidence};
            foreach my $gtypex (@{$self->inputs->{gtypex_files}}) {
            	my $meta = $gtypex->metadata;
            	$gtypex->basename =~ /\S+\.bam/;
            	my $bam_basename = $&;
                $self->throw("gtypex file has no sample metadata - cannot analyse") unless $meta->{sample};
                my $expected = $meta->{sample};
                my $req = $self->new_requirements(memory => 3900, time => 1);
                my $gtype_file = $self->output_file(output_key => 'gtype_files',
                                              basename => $bam_basename.'.gtype',
                                              type => 'txt',
                                              metadata => {sample => $expected});
                my $gtype_path = $gtype_file->path;
                my $gtypex_path = $gtypex->path;
                open my $fh,'<',$gtypex_path or $self->throw("$gtypex_path; $!");
                my $has_data = 0;
    			my ($hit1,$hit2,$gtype1,$lhood1,$gtype2,$lhood2);
				my $entrp = <$fh>;
    			while (my $line=<$fh>) {
        			if ( $has_data && defined($hit1) && defined($hit2) ) { last; }
                    if ( !($line =~ m{sample\s+(?:.*/)?(\S+?)(?:\.snp)?\s+likelihood \d+ over \d+ sites, score (\S+),} ) ) { 
            			$self->throw("Could not parse $gtypex_path: $hit1") 
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
                elsif ( $ratio<$min_ratio ) { 
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
                my $cmd = qq[echo "$gt_status" > $gtype_path];
                $self->dispatch([$cmd, $req]);
			}
		};
	}

    method outputs_definition {
        return { gtype_files => VRPipe::StepIODefinition->get(type => 'txt',
                                                                  max_files => -1,
                                                                  description => 'file of genotype status report for the sample',
                                                                  metadata => {sample => 'sample name'}
                                                                  )
         };
    }
    
    method post_process_sub {
        return sub { return 1; }
    }
    
    method description {
        return "This step takes a gtypex file with associated sample name metadata and compares the sample name with the likelihood scores to determine genotype status";
    }

    method max_simultaneous {
        return 0;
    }
}

1;
