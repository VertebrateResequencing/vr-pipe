=head1 NAME

VRPipe::Steps::bam_substitution_rate - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

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

class VRPipe::Steps::bam_substitution_rate with VRPipe::StepRole {
    use VRPipe::Parser;
    
    method options_definition {
        return { reference_fasta => VRPipe::StepOption->create(description => 'absolute path to genome reference file'),
                 samtools_exe => VRPipe::StepOption->create(description => 'path to samtools executable', 
							 optional => 1, 
							 default_value => 'samtools') };
    }
    method inputs_definition {
        return { bam_files => VRPipe::StepIODefinition->create(type => 'bam',
                                                            max_files => -1,
                                                            description => '1 or more bam files',
                                                            metadata => {mapped_fastqs => 'comma separated list of the fastq file(s) that were mapped',
                                                                         optional => ['mapped_fastqs']}),
                 dict_file => VRPipe::StepIODefinition->create(type => 'txt',
                                                            description => 'a sequence dictionary file for your reference fasta') };
    }
    method body_sub {
        return sub {
            my $self = shift;
            my $options = $self->options;
            my $ref = Path::Class::File->new($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path") unless $ref->is_absolute;
            
            my $samtools_exe = $options->{samtools_exe};
            my ($dict_file) = @{$self->inputs->{dict_file}};
            my $dict_path = $dict_file->path->stringify;
            my $req = $self->new_requirements(memory => 3000, time => 1);
            
            foreach my $bam_file (@{$self->inputs->{bam_files}}) {
                my $bam_path = $bam_file->path->stringify;
                my $this_cmd = qq[$samtools_exe mpileup -q 20 -uf $ref $bam_path | bcftools view -cgvI - | grep -vc \\#];
                $self->dispatch_vrpipecode("use VRPipe::Steps::bam_substitution_rate; VRPipe::Steps::bam_substitution_rate->calculate_substitution_rate(snp_count_from => q[$this_cmd], dict_file => q[$dict_path], bam_path => q[$bam_path ]);", $req);
            }
        };
    }
    method outputs_definition {
        return { };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Calculates the substitution rate given a mapped bam, storing the result as metadata on the bam (and and source bam/fastqs it was made from).";
    }
    method max_simultaneous {
        return 0; # meaning unlimited
    }
    
    method calculate_substitution_rate (ClassName|Object $self: Str :$snp_count_from!, Str|File :$dict_file!, Str|File :$bam_path!) {
        my $bam_file = VRPipe::File->get(path => $bam_path);
        
        $bam_file->disconnect;
        open(my $cmd_fh, "$snp_count_from |") || $self->throw("failed to run [$snp_count_from]");
        my $snps;
        while (<$cmd_fh>) {
            chomp;
            $snps = $_;
        }
        close($cmd_fh);
        
        my $pars = VRPipe::Parser->create('dict', {file => $dict_file});
        my $ref_bases = $pars->total_length;
        
        my $substitution_rate = sprintf("%0.3f", $snps / $ref_bases);
        
        $bam_file->add_metadata({substitution_rate => $substitution_rate});
        my $bam_meta = $bam_file->metadata;
        foreach my $key (qw(source_bam source_fastq mapped_fastqs)) {
            my @source_paths = split(',', $bam_meta->{$key} || next);
            foreach my $path (@source_paths) {
                VRPipe::File->get(path => $path)->add_metadata({substitution_rate => $substitution_rate});
            }
        }
    }
}

1;
