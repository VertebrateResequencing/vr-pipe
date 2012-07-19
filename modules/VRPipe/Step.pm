
=head1 NAME

VRPipe::Step - defines what work should be done, as a component of a pipeline

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

A L<VRPipe::Pipeline> is comprised of one or more Steps. A Step defines what
work should be done. It also defines what inputs it accepts, and what it
outputs, along with what options are available. Typically the code that is run
(in the C<body_sub()>) when a Step is processed by B<VRPipe> looks at the
inputs and options and comes up with a command line that needs executing. It is
important to note that Steps do not actually run the command lines themselves;
they only define what is supposed to run. Truly trivial steps may do the work
themselves only if that work will always take less than a second, regardless of
input and options.

Steps should also be written generically, with the idea that they could become
a component of many different pipelines. They should not care about where their
input files are, what their filenames are, nor where their output files will
go. The B<VRPipe> system will take care of all of this for the Step
automatically. The system also handles other things, like making sure the
command line ran successfully (exited with 0), and that the input files were
available, and that the output files were actually created.

If a Step (once it's command line is exectuted - something it doesn't have to
worry about) will result in any output files being generated (even tempory
ones), it must always tell the system by using the C<output_file()> method.
That method, along with other critical methods that should be used, can be
found in L<VRPipe::StepRole>.

Like Pipelines, A Step is stored in the database and can be created on-the-fly.
They can also be written into Perl module files. Creating a
C<VRPipe::Steps::[step_name]> class (which must implement the StepRole role)
will result in that Step automatically being entered into the database, and
updating the code in the C<body_sub()> will also update the database
automatically. Steps should be carefully named so that their utility is clear,
specific and likely to be unique without collision in the future. For example,
do not name a Step that "calls SNPs" the 'call' step, because in the future you
might want a Step that "calls Indels" and they can't both be named 'call'. Do
not even name it 'call_snps', because in the future you may want to call SNPs
in a completey different way, and you can't update the Step without affecting
previously created (possibly still running) PipelineSetups. A better name might
be 'call_snps_with_methodX'. Remember that users have to look at a list of Step
names when deciding what Steps they want in a custom Pipeline they're putting
together, so a good name is important; try to be consistent with other existing
Step names as well.

In addition to (usually) defining a command line that should be executed (at
some point), a Step also defines the L<VRPipe::Requirements> needed to run that
Job using the C<new_requirements()> method. This is just a sensible default; it
can be overridden by the user or by the system when appropriate.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

class VRPipe::Step extends VRPipe::Persistent with VRPipe::StepRole {
    use VRPipe::StepNonPersistentFactory;
    
    has 'name' => (is     => 'rw',
                   isa    => Varchar [64],
                   traits => ['VRPipe::Persistent::Attributes'],
                   is_key => 1);
    
    has 'options_definition' => (is      => 'rw',
                                 isa     => PersistentHashRef,
                                 traits  => ['VRPipe::Persistent::Attributes'],
                                 default => sub { {} });
    
    has 'inputs_definition' => (is     => 'rw',
                                isa    => PersistentHashRef,
                                traits => ['VRPipe::Persistent::Attributes']);
    
    has 'body_sub' => (is     => 'rw',
                       isa    => 'CodeRef',
                       traits => ['VRPipe::Persistent::Attributes']);
    
    has 'post_process_sub' => (is     => 'rw',
                               isa    => 'CodeRef',
                               traits => ['VRPipe::Persistent::Attributes']);
    
    has 'outputs_definition' => (is     => 'rw',
                                 isa    => PersistentHashRef,
                                 traits => ['VRPipe::Persistent::Attributes']);
    
    has 'max_simultaneous' => (is      => 'rw',
                               isa     => IntSQL [4],
                               traits  => ['VRPipe::Persistent::Attributes'],
                               default => 0); # meaning unlimited
    
    has 'description' => (is          => 'rw',
                          isa         => Text,
                          traits      => ['VRPipe::Persistent::Attributes'],
                          is_nullable => 1);
    
    __PACKAGE__->make_persistent();
}

1;
