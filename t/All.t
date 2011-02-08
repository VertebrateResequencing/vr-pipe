#!/usr/bin/env perl
# test basics of all perl modules and scripts in the distribution;
use strict;
use warnings;

use Test::Strict;

all_perl_files_ok('scripts');

# we use MooseX::Declare, so none of our modules have "use strict"
$Test::Strict::TEST_STRICT = 0;
all_perl_files_ok('modules');

exit;
