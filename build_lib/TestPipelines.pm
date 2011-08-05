package TestPipelines;
use strict;
use warnings;
use TestPersistentReal;
use Exporter 'import';
use Path::Class;
use lib "t";

our @EXPORT = qw(get_output_dir handle_pipeline output_subdirs);

our $manager = VRPipe::Manager->get();

sub get_output_dir {
    my $sub_dir = shift;
    my $scheduler = VRPipe::Scheduler->get();
    my $output_dir = dir($scheduler->output_root, 'pipelines_test_output', $sub_dir);
    $scheduler->remove_tree($output_dir);
    $scheduler->make_path($output_dir);
    return $output_dir;
}

sub handle_pipeline {
    my $give_up = 100;
    while (! $manager->trigger) {
        last if $give_up-- <= 0;
        $manager->handle_submissions;
        sleep(1);
    }
    my $all_created = 1;
    foreach my $ofile (@_) {
        unless (-s $ofile) {
            warn "$ofile is missing\n";
            $all_created = 0;
        }
    }
    return $all_created;
}

sub output_subdirs {
    my $element_id = shift;
    return ($manager->hashed_dirs('VRPipe::DataElement::'.$element_id), $element_id);
}

1;