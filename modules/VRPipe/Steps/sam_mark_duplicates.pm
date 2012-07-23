
=head1 NAME

VRPipe::Steps::sam_mark_duplicates - a step

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Shane McCarthy <sm15@sanger.ac.uk>.

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

class VRPipe::Steps::sam_mark_duplicates extends VRPipe::Steps::picard {
    around options_definition {
        return { %{ $self->$orig }, markdup_options => VRPipe::StepOption->create(description => 'command line options for Picard MarkDuplicates', optional => 1, default_value => 'ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT'), };
    }
    
    method inputs_definition {
        return { sam_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => '1 or more coordinate-sorted sam files') };
    }
    
    method body_sub {
        return sub {
            my $self    = shift;
            my $options = $self->options;
            $self->handle_standard_options($options);
            my $markdup_jar  = $self->jar('MarkDuplicates.jar');
            my $markdup_opts = $options->{markdup_options};
            
            $self->set_cmd_summary(VRPipe::StepCmdSummary->create(exe     => 'picard',
                                                                  version => $self->picard_version(),
                                                                  summary => 'java $jvm_args -jar MarkDuplicates.jar INPUT=$sam_file OUTPUT=$markdup_sam_file ' . $markdup_opts));
            
            my $req = $self->new_requirements(memory => 5800, time => 2);
            foreach my $sam (@{ $self->inputs->{sam_files} }) {
                my $sam_base     = $sam->basename;
                my $sam_meta     = $sam->metadata;
                my $markdup_base = $sam_base;
                $markdup_base =~ s/sam$/markdup.sam/;
                my $markdup_sam_file = $self->output_file(output_key => 'markdup_sam_files',
                                                          basename   => $markdup_base,
                                                          type       => 'txt',
                                                          metadata   => $sam_meta);
                
                my $temp_dir = $options->{tmp_dir} || $markdup_sam_file->dir;
                my $jvm_args = $self->jvm_args($req->memory, $temp_dir);
                
                my $this_cmd = $self->java_exe . " $jvm_args -jar $markdup_jar INPUT=" . $sam->path . " OUTPUT=" . $markdup_sam_file->path . " $markdup_opts";
                #   Doing a normal dispatch for now as issue with below
                #   $self->dispatch_wrapped_cmd('VRPipe::Steps::sam_mark_duplicates', 'markdup_and_check', [$this_cmd, $req, { output_files => [$markdup_sam_file] }]);
                $self->dispatch([qq[$this_cmd], $req, { output_files => [$markdup_sam_file] }]);
            }
        };
    }
    
    method outputs_definition {
        return { markdup_sam_files => VRPipe::StepIODefinition->create(type => 'txt', max_files => -1, description => 'a sam file with duplicates marked') };
    }
    
    method post_process_sub {
        return sub { return 1; };
    }
    
    method description {
        return "Mark duplicates in a sam files using picard";
    }
    
    method max_simultaneous {
        return 0;            # meaning unlimited
    }
    
    method markdup_and_check (ClassName|Object $self: Str $cmd_line) {
        my ($in_path, $out_path) = $cmd_line =~ /INPUT=(\S+) OUTPUT=(\S+)/;
        warn $cmd_line;
        $in_path  || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        $out_path || $self->throw("cmd_line [$cmd_line] was not constructed as expected");
        
        my $in_file  = VRPipe::File->get(path => $in_path);
        my $out_file = VRPipe::File->get(path => $out_path);
        
        $in_file->disconnect;
        warn $cmd_line;      ########
        system($cmd_line) && $self->throw("failed to run [$cmd_line]");
        
        $out_file->update_stats_from_disc(retries => 3);
        my $expected_reads = $in_file->metadata->{reads} || $in_file->num_records;
        my $actual_reads = $out_file->num_records;
        
        if ($actual_reads == $expected_reads) {
            return 1;
        }
        else {
            $out_file->unlink;
            $self->throw("cmd [$cmd_line] failed because $actual_reads reads were generated in the output sam file, yet there were $expected_reads reads in the original sam file");
        }
    }
}
1;
