#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 4;
    use VRPipeTest;
    $ENV{EMAIL_SENDER_TRANSPORT} = 'Test';
    use_ok('VRPipe::Interface::BackEnd');
}

ok my $backend = VRPipe::Interface::BackEnd->new(deployment => 'testing'), 'made a BackEnd instance';

warning_like { $backend->log('log_test', email_admin => 1, force_when_testing => 1) } qr/^\d\d\d\d\/\d\d\/\d\d \d\d:\d\d:\d\d \| pid \d+ \| log_test$/, 'log() generates a warning';
my @deliveries = Email::Sender::Simple->default_transport->deliveries;
my $email      = $deliveries[0]->{email};
my $email_body = $email->get_body;
$email_body =~ s/\s+$//;
is_deeply [scalar(@deliveries), $email->get_header("Subject"), $email->get_header('From'), $email_body], [1, 'VRPipe Server message', '"VRPipe Server" <vrpipe@do.not.reply>', 'log_test'], 'log() also sent an email';

#*** we need lots more testing of all the BackEnd methods...

done_testing;
exit;
