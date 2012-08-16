
=head1 NAME

VRPipe::Steps::cufflinks - a step

=head1 DESCRIPTION

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

class VRPipe::Steps::cufflinks with VRPipe::StepRole {
    use File::Basename;
    use Data::Dumper;
    
    method options_definition {
        return {
            cufflinks_exe    => VRPipe::StepOption->create(description => 'path to your cufflinks executable', optional => 1, default_value => 'cufflinks'),
            reference_fasta  => VRPipe::StepOption->create(description => 'path to genome file e.g. mm9.fa'),
            known_genes_path => VRPipe::StepOption->create(description => 'path to known genes file file e.g. knownGenesMm9.gtf'),
            gene_mask_path   => VRPipe::StepOption->create(description => 'path to gene mask file e.g. GeneMaskMm9.gtf')
        };
    }
    
    method inputs_definition {
        return {
            # sequence file - fastq for now
            sam_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => '1 or more sam files')
        };
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $cufflinks_exe = $options->{cufflinks_exe};
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'cufflinks', version => VRPipe::StepCmdSummary->determine_version($cufflinks_exe, 'cufflinks v(.+)\n'), summary => 'cufflinks'));
            my $req              = $self->new_requirements(memory => 500, time => 1); #16GB RAM? Could be 8GB?
            my @input_file       = @{ $self->inputs->{sam_files} };
            my $reference_fasta  = $options->{reference_fasta};
            my $known_genes_path = $options->{known_genes_path};
            my $gene_mask_path   = $options->{gene_mask_path};
            
            $self->throw("One input file expected") unless (@input_file == 1);
            my $output_file_1 = $self->output_file(
                output_key => 'transcripts_gtf',
                basename   => "transcripts.gtf",
                type       => 'txt',
                metadata   => $input_file[0]->metadata
            );
            my $output_file_2 = $self->output_file(
                output_key => 'isoforms_fpkm_tracking',
                basename   => "isoforms.fpkm_tracking",
                type       => 'txt',
                metadata   => $input_file[0]->metadata
            );
            my $output_file_3 = $self->output_file(
                output_key => 'genes_fpkm_tracking',
                basename   => "genes.fpkm_tracking",
                type       => 'txt',
                metadata   => $input_file[0]->metadata
            );
            
            my $output_file_dir = $output_file_1->dir->stringify;
            my $input_file_path = $input_file[0]->path;
            
            #TO DO:  Build the command up  - allow more options.
            
            my $cmd = "cufflinks --num-threads 32 --frag-bias-correct $reference_fasta --GTF-guide $known_genes_path --upper-quartile-norm --max-mle-iterations 10000 --total-hits-norm --max-bundle-frags 100000000 -M $gene_mask_path -o $output_file_dir --multi-read-correct $input_file_path";
            
            my $out = $self->dispatch([qq[$cmd], $req, { output_files => [$output_file_1, $output_file_2, $output_file_3] }]);
          }
    }
    
    method outputs_definition {
        return {
            transcripts_gtf        => VRPipe::StepIODefinition->create(type => 'txt', description => 'cufflinks assembled isoforms in gtf format'),
            isoforms_fpkm_tracking => VRPipe::StepIODefinition->create(type => 'txt', description => 'cufflinks estimated isoform-level expression values in FPKM Tracking Format'),
            genes_fpkm_tracking    => VRPipe::StepIODefinition->create(type => 'txt', description => 'cufflinks estimated gene-level expression values in the generic FPKM Tracking Format')
        };
    }
    
    method description {
        return "Step for cufflinks: transcript assembly and quantification for rna-seq";
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
}
