
=head1 NAME

VRPipe::Steps::bam_reheader - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

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

class VRPipe::Steps::bam_reheader with VRPipe::StepRole {
    method options_definition {
        return {
            samtools_exe => VRPipe::StepOption->create(
                description   => 'path to your samtools executable',
                optional      => 1,
                default_value => 'samtools'
            ),
            header_comment_file => VRPipe::StepOption->create(description => 'path to your file containing SAM comment lines to include in the header', optional => 1)
        };
    }
    
    method inputs_definition {
        return {
            bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => '1 or more bam files',
                metadata    => {
                    lane             => 'lane name (a unique identifer for this sequencing run, aka read group)',
                    library          => 'library name',
                    sample           => 'sample name',
                    center_name      => 'center name',
                    platform         => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                    study            => 'name of the study',
                    insert_size      => 'expected library insert size if paired',
                    mean_insert_size => 'calculated mean insert size if paired',
                    bases            => 'total number of base pairs',
                    reads            => 'total number of reads (sequences)',
                    paired           => '0=unpaired reads were mapped; 1=paired reads were mapped',
                    optional         => ['bases', 'lane', 'library', 'insert_size', 'mean_insert_size', 'sample', 'center_name', 'platform', 'study']
                }
            ),
            dict_file => VRPipe::StepIODefinition->create(
                type        => 'txt',
                description => 'a sequence dictionary file for your reference fasta'
            )
        };
    }
    
    method body_sub {
        return sub {
            my $self      = shift;
            my $options   = $self->options;
            my $samtools  = $options->{samtools_exe};
            my $dict_path = $self->inputs->{dict_file}->[0]->path;
            my $comment   = '';
            if ($options->{header_comment_file}) {
                my $comment_path = file($options->{header_comment_file});
                $self->throw("header_comment_file must be an absolute path if it is supplied") unless $comment_path->is_absolute;
                $comment .= ", comment => q[$comment_path]";
            }
            
            my $req = $self->new_requirements(memory => 1000, time => 1);
            my $step_state = $self->step_state->id;
            foreach my $bam (@{ $self->inputs->{bam_files} }) {
                my $bam_path        = $bam->path;
                my $bam_meta        = $bam->metadata;
                my $basename        = $bam->basename;
                my $headed_bam_file = $self->output_file(
                    output_key => 'headed_bam_files',
                    basename   => $basename,
                    type       => 'bam',
                    metadata   => $bam_meta
                );
                
                my $headed_bam_path = $headed_bam_file->path;
                
                my $header_file = $self->output_file(
                    basename  => $basename . '.header',
                    type      => 'txt',
                    temporary => 1
                );
                
                my $this_cmd = "use VRPipe::Steps::bam_reheader; VRPipe::Steps::bam_reheader->reheader_and_check(samtools => q[$samtools], dict => q[$dict_path], output => q[$headed_bam_path], step_state => $step_state, bam => q[$bam_path]$comment);";
                $self->dispatch_vrpipecode($this_cmd, $req, { output_files => [$headed_bam_file, $header_file] });
            }
        };
    }
    
    method outputs_definition {
        return {
            headed_bam_files => VRPipe::StepIODefinition->create(
                type        => 'bam',
                max_files   => -1,
                description => 'a bam file with good header',
                metadata    => {
                    lane             => 'lane name (a unique identifer for this sequencing run, aka read group)',
                    library          => 'library name',
                    sample           => 'sample name',
                    center_name      => 'center name',
                    platform         => 'sequencing platform, eg. ILLUMINA|LS454|ABI_SOLID',
                    study            => 'name of the study, put in the DS field of the RG header line',
                    insert_size      => 'expected library insert size if paired',
                    mean_insert_size => 'calculated mean insert size if paired',
                    bases            => 'total number of base pairs',
                    reads            => 'total number of reads (sequences)',
                    paired           => '0=unpaired reads were mapped; 1=paired reads were mapped; 2=mixture of paired and unpaired reads were mapped',
                    optional         => ['bases', 'lane', 'library', 'insert_size', 'mean_insert_size', 'sample', 'center_name', 'platform', 'study']
                }
            )
        };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Replaces a bam header so that it has complete sequence information, a good RG line, and chained PG lines";
    }
    
    method max_simultaneous {
        return 0;          # meaning unlimited
    }
    
    method reheader_and_check (ClassName|Object $self: Str|File :$samtools!, Str|File :$dict!, Str|File :$output!, Persistent :$step_state!, Str|File :$bam!, Str|File :$comment?) {
        # make a nice sam header
        my $header_file = VRPipe::File->get(path => $output . '.header');
        my $header_path = $header_file->path;
        my $hfh         = $header_file->openw;
        
        my $header_lines = 0;
        print $hfh "\@HD\tVN:1.0\tSO:coordinate\n";
        $header_lines++;
        
        # copy over the SQ lines from the dict file
        my $dict_file = VRPipe::File->get(path => $dict);
        my $dfh = $dict_file->openr;
        while (<$dfh>) {
            next unless /^\@SQ/;
            print $hfh $_;
            $header_lines++;
        }
        $dict_file->close;
        
        my $bam_file        = VRPipe::File->get(path => $bam);
        my $headed_bam_file = VRPipe::File->get(path => $output);
        
        my %tech_to_platform = (
            SLX       => 'ILLUMINA',
            '454'     => 'LS454',
            SOLID     => 'ABI_SOLID',
            ILLUMINA  => 'ILLUMINA',
            'LS454'   => 'LS454',
            ABI_SOLID => 'ABI_SOLID',
        );
        
        my $existing_meta = $bam_file->metadata;
        if (exists $tech_to_platform{ $existing_meta->{platform} }) {
            $existing_meta->{platform} = $tech_to_platform{ $existing_meta->{platform} };
        }
        
        # construct the RG lines from the bam metadata if the lane metadata is present
        # otherwise copy the RG lines from the existing header
        if (defined $existing_meta->{lane} && !($existing_meta->{lane} =~ /,/)) {
            print $hfh "\@RG\tID:" . $existing_meta->{lane};
            if (defined $existing_meta->{library}) {
                print $hfh "\tLB:" . $existing_meta->{library};
            }
            if (defined $existing_meta->{sample}) {
                print $hfh "\tSM:" . $existing_meta->{sample};
            }
            if (defined $existing_meta->{insert_size}) {
                print $hfh "\tPI:" . sprintf("%0.0f", $existing_meta->{insert_size});
            }
            elsif (defined $existing_meta->{mean_insert_size}) {
                print $hfh "\tPI:" . sprintf("%0.0f", $existing_meta->{mean_insert_size});
            }
            if (defined $existing_meta->{center_name}) {
                print $hfh "\tCN:" . $existing_meta->{center_name};
            }
            if (defined $existing_meta->{platform}) {
                print $hfh "\tPL:" . $existing_meta->{platform};
            }
            if (defined $existing_meta->{study}) {
                print $hfh "\tDS:" . $existing_meta->{study};
            }
            print $hfh "\n";
            $header_lines++;
        }
        else {
            my $pars = VRPipe::Parser->create('bam', { file => $bam_file });
            my %readgroup_info = $pars->readgroup_info();
            $self->throw("No readgroup info found and no single lane metadata for bam $bam. Cannot reheader.") unless (scalar keys %readgroup_info > 0);
            foreach my $rg (sort keys %readgroup_info) {
                print $hfh "\@RG\tID:" . $rg;
                if (defined $readgroup_info{$rg}->{LB}) {
                    print $hfh "\tLB:" . $readgroup_info{$rg}->{LB};
                }
                if (defined $readgroup_info{$rg}->{SM}) {
                    print $hfh "\tSM:" . $readgroup_info{$rg}->{SM};
                }
                if (defined $readgroup_info{$rg}->{PI}) {
                    print $hfh "\tPI:" . sprintf("%0.0f", $readgroup_info{$rg}->{PI});
                }
                if (defined $readgroup_info{$rg}->{CN}) {
                    print $hfh "\tCN:" . $readgroup_info{$rg}->{CN};
                }
                if (defined $readgroup_info{$rg}->{PL}) {
                    print $hfh "\tPL:" . $readgroup_info{$rg}->{PL};
                }
                if (defined $readgroup_info{$rg}->{DS}) {
                    print $hfh "\tDS:" . $readgroup_info{$rg}->{DS};
                }
                if (defined $readgroup_info{$rg}->{PU}) {
                    print $hfh "\tPU:" . $readgroup_info{$rg}->{PU};
                }
                if (defined $readgroup_info{$rg}->{DT}) {
                    print $hfh "\tDT:" . $readgroup_info{$rg}->{DT};
                }
                print $hfh "\n";
                $header_lines++;
            }
        }
        
        # Construct a chain of PG lines for the header by looking at previous
        # steps in our pipeline. If this the bam was produced from a vrpipe
        # datasource, find it's parent elements and include in the program chain
        # If this bam wasn't produced by vrpipe, include any original PG lines
        # that were intitially already present
        my $this_step_state = VRPipe::StepState->get(id => $step_state);
        my $pipelinesetup   = $this_step_state->pipelinesetup;
        my $dataelement     = $this_step_state->dataelement;
        my $stepmember      = $this_step_state->stepmember;
        my $this_stepm_id   = $stepmember->id;
        my $pipeline        = $stepmember->pipeline;
        my @program_chain   = $self->program_chain(pipelinesetup => $pipelinesetup, stepmember => $stepmember, dataelement => $dataelement);
        if (defined $existing_meta->{original_pg_chain}) {
            my %pg_lines;
            foreach my $link (@program_chain) {
                $pg_lines{ $link->{program_id} }->{ $link->{previous_program} }->{ $link->{program_name} }->{ $link->{program_version} }->{ $link->{command_line} } = 1;
            }
            
            my @orig_chain;
            foreach my $pg (split("\n", $existing_meta->{original_pg_chain})) {
                my ($pid) = $pg =~ /\tID:([^\t]+)/;
                my ($pn)  = $pg =~ /\tPN:([^\t]+)/;
                my ($pv)  = $pg =~ /\tVN:([^\t]+)/;
                my ($cl)  = $pg =~ /\tCL:([^\t]+)/;
                my ($pp)  = $pg =~ /\tPP:([^\t]+)/;
                $pp ||= 'null';
                push(@orig_chain, { program_id => $pid, program_name => $pn, program_version => $pv, command_line => $cl, previous_program => $pp }) unless exists $pg_lines{$pid}->{$pp}->{$pn}->{$pv}->{$cl};
                $pg_lines{$pid}->{$pp}->{$pn}->{$pv}->{$cl} = 1;
            }
            
            if (@orig_chain) {
                if (@program_chain) {
                    $program_chain[0]->{previous_program} = $orig_chain[-1]->{program_id};
                }
                
                @program_chain = (@orig_chain, @program_chain);
            }
        }
        foreach my $pg_line (@program_chain) {
            print $hfh "\@PG\tID:" . $pg_line->{program_id};
            if ($pg_line->{program_name}) {
                print $hfh "\tPN:" . $pg_line->{program_name};
            }
            if ($pg_line->{previous_program} ne 'null') {
                print $hfh "\tPP:" . $pg_line->{previous_program};
            }
            if ($pg_line->{program_version}) {
                print $hfh "\tVN:" . $pg_line->{program_version};
            }
            if ($pg_line->{command_line}) {
                print $hfh "\tCL:", $pg_line->{command_line};
            }
            print $hfh "\n";
            $header_lines++;
        }
        if ($comment) {
            my $comment_file = VRPipe::File->create(path => $comment);
            my $cfh = $comment_file->openr;
            while (<$cfh>) {
                next unless /^\@CO/;
                $self->throw("Bad line in header comment file $comment: $_") unless /^\@CO\t/;
                print $hfh $_;
                $header_lines++;
            }
            $comment_file->close;
        }
        $header_file->close;
        
        my $cmd_line = qq[$samtools reheader $header_path $bam > $output];
        $headed_bam_file->disconnect;
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        my $expected_lines = $existing_meta->{reads} + $header_lines;
        $headed_bam_file->update_stats_from_disc(retries => 3);
        my $actual_lines = $headed_bam_file->lines;
        
        if ($actual_lines == $expected_lines) {
            return 1;
        }
        else {
            $headed_bam_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_lines lines were generated in the reheaded bam file, yet there were $expected_lines records in the input bam files and header");
        }
    }
    
    method program_chain (ClassName|Object $self: VRPipe::PipelineSetup :$pipelinesetup!, VRPipe::DataElement :$dataelement!, VRPipe::StepMember :$stepmember!) {
        my %chains;
        foreach my $program ($self->command_history(pipelinesetup => $pipelinesetup, stepmember => $stepmember, dataelement => $dataelement)) {
            foreach my $readgroup (@{ $program->{element_readgroups} }) {
                push @{ $chains{$readgroup} }, $program;
            }
        }
        my @readgroups = sort keys %chains;
        
        my @program_chain;
        my %pg_lines;
        my %pp;
        my $branch = 1;
        my $dup_id = 1;
        while (@{ $chains{ $readgroups[-1] } }) {
            foreach my $readgroup (@readgroups) {
                my $program = shift @{ $chains{$readgroup} };
                
                my $pid = $program->{program_id};
                my $pn  = $program->{program_name};
                my $pv  = $program->{program_version};
                my $cl  = $program->{command_line};
                my $pp  = $pp{$readgroup} || 'null';
                
                if ($pp =~ /[^\.]+([\.\d]+)$/) {
                    $pid .= $1;
                }
                
                if (exists $pg_lines{"$pid"} && !(exists $pg_lines{"$pid"}->{"$pp"}->{"$pn"}->{"$pv"}->{"$cl"})) {
                    if (exists $pg_lines{"$pid"}->{"$pp"}) {
                        while (exists $pg_lines{"$pid.$dup_id"}) {
                            last if (exists $pg_lines{"$pid.$dup_id"}->{"$pp"}->{"$pn"}->{"$pv"}->{"$cl"});
                            $dup_id++;
                        }
                        $pid .= ".$dup_id";
                    }
                    else {
                        $pid .= ".$branch";
                        $branch++;
                    }
                }
                
                my $pg_hash = { program_id => $pid, program_name => $pn, program_version => $pv, command_line => $cl, previous_program => $pp };
                push @program_chain, $pg_hash unless (exists $pg_lines{"$pid"}->{"$pp"}->{"$pn"}->{"$pv"}->{"$cl"});
                $pg_lines{"$pid"}->{"$pp"}->{"$pn"}->{"$pv"}->{"$cl"} = 1;
                $pp{$readgroup} = $pid;
            }
        }
        return @program_chain;
    }
    
    method command_history (ClassName|Object $self: VRPipe::PipelineSetup :$pipelinesetup!, VRPipe::DataElement :$dataelement!, VRPipe::StepMember :$stepmember!) {
        my $this_stepm_id = $stepmember->id;
        my $pipeline      = $stepmember->pipeline;
        
        my @history;
        if ($pipelinesetup->datasource->type eq 'vrpipe') {
            my $vrpipe_sources = $pipelinesetup->datasource->_source_instance->vrpipe_sources;
            foreach my $link (VRPipe::DataElementLink->search({ child => $dataelement->id }, { prefetch => [qw(pipelinesetup parent)] })) {
                my $this_pipelinesetup = $link->pipelinesetup;
                my $setup_id           = $link->pipelinesetup->id;
                next unless exists $vrpipe_sources->{$setup_id};
                my $this_stepmember = VRPipe::StepMember->get(id => $vrpipe_sources->{$setup_id}->{final_smid});
                push @history, $self->command_history(pipelinesetup => $this_pipelinesetup, dataelement => $link->parent, stepmember => $this_stepmember);
            }
        }
        
        my @readgroups = $self->element_readgroups($dataelement);
        return () unless @readgroups;
        foreach my $stepm ($pipeline->step_members) {
            my ($step_state) = VRPipe::StepState->search({ pipelinesetup => $pipelinesetup, stepmember => $stepm, dataelement => $dataelement });
            next unless $step_state;
            my $cmd_summary = $step_state->cmd_summary;
            next unless $cmd_summary;
            my $step_name = $stepm->step->name;
            my $pg_hash = { program_id => $step_name, program_name => $cmd_summary->exe, program_version => $cmd_summary->version, command_line => $cmd_summary->summary, element_readgroups => \@readgroups };
            push @history, $pg_hash;
            last if $stepm->id == $this_stepm_id;
        }
        return @history;
    }
    
    method element_readgroups (ClassName|Object $self: VRPipe::DataElement $dataelement!) {
        my %readgroups;
        my $files = $dataelement->files || $self->throw("data element " . $dataelement->id . " gave a result with no paths");
        foreach my $bam (@$files) {
            my $metadata = $bam->metadata;
            next unless defined $metadata->{lane};
            foreach my $lane (split ',', $metadata->{lane}) {
                $readgroups{$lane} = 1;
            }
        }
        foreach my $link (VRPipe::DataElementLink->search({ child => $dataelement->id }, { prefetch => 'parent' })) {
            foreach my $rg ($self->element_readgroups($link->parent)) {
                $readgroups{$rg} = 1;
            }
        }
        return sort keys %readgroups;
    }
}

1;
