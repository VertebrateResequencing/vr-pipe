#!/usr/bin/env perl
use strict;
use warnings;
use POE qw(Filter::Reference);

BEGIN {
    use Test::Most tests => 6;
    
    use_ok('t::VRPipe::SimpleSpawner');
    use_ok('t::VRPipe::LivingSpawner');
    use_ok('VRPipe::Cmd');
    use_ok('MooseX::Workers::Job');
}

ok my $ss = t::VRPipe::SimpleSpawner->new(), 'made a new SimpleSpawner using default processes';
$ss->run();
is_deeply $ss->output, ['simple'], 'output from a single process was ok';

package Manager;
use Moose;
with qw(MooseX::Workers);
use POE qw(Filter::Reference Filter::Line);

sub run {
    warn "will spawn and run\n";
    $_[0]->spawn(
        sub {
            sleep 3;
            print POE::Filter::Reference->new->put([ {msg => "Hello World\n"} ]);
            print STDERR "Hey look, an error message\n";
        }
    );
    POE::Kernel->run();
}

sub stdout_filter  { new POE::Filter::Reference }
sub stderr_filter  { new POE::Filter::Line }

sub worker_stdout  {  
    my ( $self, $result ) = @_;
    warn "in worker_stdout\n";
    print $result->{msg};
}
sub worker_stderr  {
    my ( $self, $stderr_msg ) = @_;
    warn $stderr_msg;
}

sub worker_manager_start {  }
sub worker_manager_stop  {  }
sub max_workers_reached  {  }
sub worker_error   {  }
sub worker_done    {  }
sub worker_started {  }
sub sig_child      {  }
sub sig_TERM       {  }
no Moose;

package main;
Manager->new->run();

exit;

my @jobs;
foreach my $name ('foo', 'bar', 'henry') {
    my $job = MooseX::Workers::Job->new(
        name    => $name,
        command => sub { for (1..3) { sleep 1; print POE::Filter::Reference->new->put([ {msg => $name} ]); } }
    );
    push(@jobs, $job);
}

ok $ss = t::VRPipe::SimpleSpawner->new(processes => \@jobs), 'made a new SimpleSpawner';
$ss->run();
is_deeply $ss->output, ['foo', 'bar', 'henry', 'foo', 'bar', 'henry', 'foo', 'bar', 'henry'], 'output from 3 processes interleaved as expected';

exit;