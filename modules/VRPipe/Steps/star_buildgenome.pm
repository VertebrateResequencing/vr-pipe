
=head1 NAME

VRPipe::Steps::star_buildgenome - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Yasin Memari <ym3@sanger.ac.uk>.

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

class VRPipe::Steps::star_buildgenome with VRPipe::StepRole {
    method options_definition {
        return {
            star_exe => VRPipe::StepOption->create(
                description   => 'path to your STAR executable',
                optional      => 1,
                default_value => 'STAR'
            ),
            star_genomeGenerate_options => VRPipe::StepOption->create(
                description   => 'options to STAR genomeGenerate, excluding genome directory option (default specifies parameters required to estimate the memory size)',
                optional      => 1,
                default_value => '--genomeSAindexNbases 14 --runThreadN 16 --limitIObufferSize 150000000'
            ),
            reference_fasta => VRPipe::StepOption->create(description => 'absolute path to genome reference file'),
        };
    }
    
    method inputs_definition {
        return {};
    }
    
    method body_sub {
        return sub {
            my $self          = shift;
            my $options       = $self->options;
            my $star_exe      = $options->{star_exe};
            my $star_gen_opts = $options->{star_genomeGenerate_options};
            if ($star_gen_opts =~ /runMode|genomeDir/) {
                $self->throw("star_genomeGenerate_options should not include --runMode and --genomeDir");
            }
            if ($star_gen_opts !~ /genomeSAindexNbases/ || $star_gen_opts !~ /runThreadN/ || $star_gen_opts !~ /limitIObufferSize/) {
                $self->throw("star_genomeGenerate_options should include --genomeSAindexNbases, --limitIObufferSize and --runThreadN");
            }
            
            my $ref = file($options->{reference_fasta});
            $self->throw("reference_fasta must be an absolute path, $ref") unless $ref->is_absolute;
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe => 'STAR', version => 0, summary => 'STAR --runMode genomeGenerate --genomeDir ' . $ref->dir . ' --genomeFastaFiles $reference_fasta ' . $star_gen_opts));
            
            $self->output_file(output_key => 'params_file', output_dir => $ref->dir, basename => 'genomeParameters.txt', type => 'txt');
            $self->output_file(output_key => 'log_file', output_dir => $ref->dir, basename => 'Log.out', type => 'txt');
            my @chr_affix = qw(Start Name NameLength Length);
            foreach (@chr_affix) {
                $self->output_file(output_key => 'text_files', output_dir => $ref->dir, basename => "chr$_.txt", type => 'txt');
            }
            foreach (qw(Genome SA SAindex)) {
                $self->output_file(output_key => 'binary_files', output_dir => $ref->dir, basename => $_, type => 'bin');
            }
            $self->output_file(output_key => 'sjdb_file', output_dir => $ref->dir, basename => "sjdbList.out.tab", type => 'txt');
            
            my ($cpus)         = $star_gen_opts =~ m/--runThreadN\s*(\d+)/;
            my ($SAindexBases) = $star_gen_opts =~ m/--genomeSAindexNbases\s*(\d+)/;
            my ($buffer)       = $star_gen_opts =~ m/--limitIObufferSize\s*(\d+)/;
            my $size           = -s $ref;
            my $memory         = int((10 * $size + 6 * 4**$SAindexBases + $buffer * $cpus) / 1024 / 1024); #MB
            my $req = $self->new_requirements(memory => $memory, time => 4, $cpus ? (cpus => $cpus) : ());
            
            $self->dispatch(["cd " . $ref->dir . "; $star_exe --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles $ref $star_gen_opts", $req, { block_and_skip_if_ok => 1 }]);
        };
    }
    
    method outputs_definition {
        return {
            params_file  => VRPipe::StepIODefinition->create(type => 'txt', description => 'file containing genome generation parameters'),
            log_file     => VRPipe::StepIODefinition->create(type => 'txt', description => 'log file'),
            text_files   => VRPipe::StepIODefinition->create(type => 'txt', description => 'other text files containing chromosome names and lengths'),
            binary_files => VRPipe::StepIODefinition->create(type => 'bin', description => 'binary genome sequence and suffix array files'),
            sjdb_file    => VRPipe::StepIODefinition->create(type => 'txt', description => 'splice junction database files', check_existence => 0),
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Builds a STAR genome sequence and suffix array from a fasta file.";
    }
    
    method max_simultaneous {
        return 0;                                                                                          # meaning unlimited
    }
}

1;
