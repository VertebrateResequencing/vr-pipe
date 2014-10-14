VRPipe is a generic pipeline system designed for use by the Vertebrate
Resequencing team at the Sanger Institute, but also for broad use for anyone.

See the [wiki](https://github.com/VertebrateResequencing/vr-pipe/wiki) for
further details, including comprehensive installation and usage guides. There
is also a complete guide for using VRPipe in Amazon's cloud, suitable even for
those who have never used the cloud before.

The rest of this document covers the essentials, and includes some installation
advice (how to install problematic CPAN dependencies) not present on the wiki.

NB: when upgrading, always read the IMPORTANT_NOTES file and follow the guidance
for all versions since your currently installed version.


# Install Pre-Requisites

You will need the source version of samtools compiled with -fPIC and -m64 in the
CFLAGS, and the environment variable SAMTOOLS pointing to that source directory
(which should now contain bam.h and libbam.a, and MUST contain the samtools
executable). Samtools can be downloaded from here:
http://sourceforge.net/projects/samtools/files/samtools/
(Note that tests use the samtools in your PATH, not the samtools in your
SAMTOOLS dir; production pipelines will use the samtools at the absolute path
you configure them with during setup - the default will be the first one in your
PATH.)

You need Redis 2.6 or greater installed. Compile it as normal ('make') and
simply include redis-server in your PATH. You do not have to configure it in any
way. But note that its maximum number of connections is capped by the max file
descriptors limit imposed on the user that will run vrpipe-server, so make sure
that limit (ulimit -n) is higher than the number of CPUs in your cluster.
Redis can be downloaded from here: http://redis.io/download

It is recommended that you set PERL_INLINE_DIRECTORY to ~/.Inline and create
that directory.

You require Module::Build in order for the Build.PL script to work. It is
recommended you 'install Bundle::CPAN' using 'cpan' prior to attempting setup.


First, do the initial setup and dependency installation:  
`perl Build.PL`

If this says you have "ERRORS/WARNINGS FOUND IN PREREQUISITES" try:  
`./Build installdeps`  
to install missing prerequisites from CPAN.

installdeps will install things automatically in some random order, but ideally
you should manually install certain modules in a specific order:

1. EV
2. AnyEvent

Possible problems you may encounter in installing your missing prerequisites:

* You may find that you have to unset your LANG environment variable temporarily
  to install all CPAN prerequisites.
* You will likely find that Inline::Filters fails its preprocess.t test; this is
  ok and you can manually use `cpan` to:  
  `cpan> force install Inline::Filters`
* If you have trouble with Authen::PAM (required by Authen::Simple::PAM), you
  may have to install it manually instead of via a non-interactive system like
  cpan minus. Enter the build directory where it tried to build Authen::PAM,
  then do:
  `perl Makefile.PL`
  `make`
  `make test`
  `make install`  
* You may find that MooseX::Types::Parameterizable fails its tests. In that case
  you can manually use 'cpan' to install an earlier version that should work:  
  `cpan> install JJNAPIORK/MooseX-Types-Parameterizable-0.07.tar.gz`
* If you find that AnyEvent::ForkManager fails tests and emits errors
  mentioning WNOHANG, you will have to cd to the cpan build directory and edit
  line 11 of blib/lib/AnyEvent/ForkManager.pm to read "use POSIX qw(WNOHANG);"
  (without the quotes). Then you can do `make test` and `make install`.
* Twiggy may fail its disconnect.t test. This is safe to ignore; just:  
  `cpan> force install Twiggy`
* You may encounter difficulties installing Proc::ProcessTable, which is one of
  Proc::Killfam's dependencies. In that case, try manually installing an earlier
  version of it before retrying the Proc::Killfam install normally:  
  `cpan> install DURIST/Proc-ProcessTable-0.44.tar.gz`
* Some installs of Perl 5.8.8 can have trouble with Class::Accessor::Grouped (a
  prerequisite of DBIx::Class), where some tests fail due to desired error
  messages not looking quite right. This should be harmless in reality, so just
  force install it:  
  `cpan> force install Class::Accessor::Grouped`  
  and then try the DBIx::Class installation again normally.
* Time::Format may fail some tests with errors mentioning problems with
  Date::Manip. As long as you actually have Date::Manip installed, these
  failures are probably fine and you can force the install:  
  `cpan> force install Time::Format`
* Redis may fail one its tests due to an error message not being in the
  expected format. This is safe to ignore:  
  `cpan> force install Redis`
* If all the prerequisites seem to be installed but running 'perl Build.PL'
  gives you an error message mentioning problems with any kind of Moose module,
  one possibility is that you have a mixture of some new Moose modules and old
  ones. You can solve this by updating all your Moose-related modules to their
  latest versions. Eg. if you have cpanm installed, try:  
  `cpanm Class::MOP Fey Fey::ORM MooseX::Aliases MooseX::ClassAttribute
    MooseX::LazyRequire MooseX::Method::Signatures MooseX::NonMoose MooseX::POE
    MooseX::Role::Cmd MooseX::Role::Parameterized MooseX::SemiAffordanceAccessor
    MooseX::Singleton MooseX::StrictConstructor`
* If you get as far as running `./Build test` (see below), and find some tests
  failing with an error message that mentions IO::Uncompress::Gunzip, try
  installing the latest version of that:  
  `cpan> install IO::Uncompress::Gunzip`


# Configure

Once the prerequisites are installed you will be asked a number of setup
questions such as the details of your production and testing databases. It is
recommended to use a stand-alone relational database such as MySQL. If you're
only using code such as a parser or the bas function, or only running pipelines
on very small datasets, you can use SQLite, supplying a file path for the
database name. However SQLite is not recommended for production use.

Note that when using MySQL, the server should be configured to allow at least
1/3 as many connections as you have CPUs in your cluster. It must also allow
the transaction isolation level to be changed from a default of REPEATABLE-READ
to READ-COMMITTED, and gap-locks must be possible
(innodb_locks_unsafe_for_binlog must be off). In practice this means that if
you have database replication set up, you may find that you have to reconfigure
the server to binlog-mode=MIXED and restart/redo your replication server.
To avoid possible issues we recommend MySQL 5.5 with the following settings in
my.cnf, where the idea is to maximise transaction concurrency and performance:

* max_connections = 1000
* binlog_cache_size=32M
* innodb_log_file_size=512M
* innodb_log_files_in_group = 2
* innodb_file_per_table = 1
* innodb_fast_shutdown=0
* innodb_buffer_pool_size = 1G
* innodb_change_buffering=all
* innodb_flush_log_at_trx_commit=0
* innodb_log_buffer_size = 32M
* innodb_thread_concurrency=16
* innodb_concurrency_tickets=5000
* innodb_commit_concurrency=50
* innodb_autoinc_lock_mode=2
* binlog_format=mixed
* transaction_isolation = REPEATABLE-READ
* innodb_locks_unsafe_for_binlog = 0

The SGE (Sun|Oracle GridEngine) scheduler has only been tested with SGE 8.1.3
from: https://arc.liv.ac.uk/trac/SGE. If tests get stuck and nothing happens
with no obvious error messages, try running:  
`perl -Imodules -It scripts/vrpipe-server --deployment testing --farm testing_farm -f --debug start`  
which should provide an error message. You may then have to alter or subclass
VRPipe::Schedulers::sge; bug reports or patches welcome.

Once setup is complete modules/VRPipe/SiteConfig.pm will have been created. It
is this file that holds your site-wide configuration of VRPipe.

# Test

To test the code prior to using it:  
`./Build test`

(Note that there is currently an issue where the first time you ever run tests
they will fail because the test database wasn't created yet; start the tests
running until you get past the DataSource.t test, then ctrl-c to kill the tests,
then start them again: they should now work)

A number of environment variables affect which tests run and how. In most cases
it is fine to not worry about these extra variables. They are described here for
the sake of completeness.

* VRPIPE_TEST_PIPELINES, when true, will fully test all the pipelines and steps
                       instead of just the core functionality of the system.
                       These tests take a very long time - potentially hours -
                       so this is off by default. If you're interested in
                       running a particular pipeline, you can test just that
                       pipeline by setting this variable to true and then
                       running:
               `./Build test --test_files t/VRPipe/Pipelines/[name].t --verbose`
* GATK, pointing to a directory containing GATK jar files, will enable tests that
      require GATK
* GATK2, pointing to a directory containing GATK version 2+ jar files, will enable
      tests that require GATK version 2 or higher
* PICARD, pointing to a directory containing Picard jar files, will enable tests
      that require Picard
* CRAMTOOLS, pointing to a directory containing Cramtools jar files, will enable
      tests that require Cramtools
* TRIMMOMATIC_JAR_PATH, pointing to the timmomatic .jar file, will enable tests
      that require Trimmomatic
* CONVEX_R_LIB, pointing to a directory containing your R lib files, will enable
      tests of Convex
* VRPIPE_VRTRACK_TESTDB will enable testing of the VRTrack DataSource, using the
                      value supplied as the database name (other database
                      details will come from the standard VRTrack environment
                      variables)
                      
Additionally, tests that require certain executables will not run unless you
have those executables in your PATH.


# Install

`./Build install`

(or just include the modules subdirectory in your PERL5LIB, and include the
scripts subdirectory in your PATH)

To create your production database (testing database is created automatically):  
`vrpipe-db_deploy`

If in the future the VRPipe code is updated and there is a change to the schema,
you will need to stop VRPipe, update your code, then run:  
`vrpipe-db_upgrade`


# Usage

Ensure that vrpipe-server is always running on a machine that can submit jobs
to your job scheduler (it must also be able to ssh to all other machines without
a password), and that has at least 2GB of free memory (which probably means the
machine must have at least 3GB physical memory), especially if using the local
scheduler:  
`vrpipe-server --farm my_farm start`

Now you can use `vrpipe-create_step` and `vrpipe-create_pipeline` to create
pipelines if necessary. Then end-users can set up a pipeline using `vrpipe-setup`.
vrpipe-setup and vrpipe-status will be the most commonly used scripts. The web
front-end is currently limited to providing status, but it does so a bit nicer
and a bit faster than vrpipe-status cmd line tool. You can alter what it
displays by manually appending vrpipe-status arguments to the url, eg.
http://[serverurl]/status?brief=1&incomplete=1&user=me

Access to output files can be achieved with `vrpipe-output`. When you have an
output file but can't remember how it was made, use `vrpipe-fileinfo`.

Each vrpipe-* script has a --help option, eg:  
`vrpipe-status --help`

As a general point, any Perl script you write that wants to use some VRPipe
code should typically "use VRPipe;". After that most things should work without
further "use" statements. To choose the deployment you can "import" the one you
want. eg. to work with your testing database:  
`use VRPipe 'testing';`  
or on the command-line:  
`perl -MVRPipe=testing -e '...'`


# When Things Go Wrong

For normal errors caused by failing jobs in a pipeline, and for most other
issues, you should receive an email that either specifies the problem, or asks
you to investigate (and fix) with vrpipe-submissions. Using vrpipe-status will
also point out issues in case you're not receiving emails.

In the rare event that vrpipe-status and vrpipe-submissions do not reveal the
problem, you can investigate further by looking at the logs. vrpipe-server keeps
a log file in the logging directory you configured during installation (when you
ran 'perl Build.PL'). For problems with a particular setup, use vrpipe-logs to
see what has been happening with it.

Sometimes it can happen that a step is failing because the output of a previous
step is corrupt. VRPipe will only ever automatically keeps retrying the current
step, so it will keep failing. vrpipe-submissions only gives you an easy way to
force the reattempt of a failed submission (ie. which will be for the current
step), so it isn't helpful either. In this case you should use vrpipe-elements,
which provides an easy way to 'go back' and reset steps that VRPipe thought
completed successfully. Just be careful not to redo any step that outputs files
that might be used by other running jobs in your setup, or by other setups.


# External Software

Some software need environment variables setup (using setenv in csh or export in
bash). The following list shows the name of the software, the environment
variable you need to set, and the value you should set it to, separated by
commas.

* samtools,SAMTOOLS,/path/to/samtools/source_directory
* cramtools,CRAMTOOLS,/path/to/cramtools_jar_files

eg. to have cramtools work properly in the pipelines you might do:  
`setenv CRAMTOOLS /path/to/cramtools_jar_files`  
or  
`export CRAMTOOLS=/path/to/cramtools_jar_files`  
depending on what shell you are using.

These values are used internally and are independent of values you might set
when asked for eg. the path to your cramtools jars when running vrpipe-setup.

A number of other environment variables are used to setup default values for
when running vrpipe-setup, but are not required to use the software (ie. you can
just manually enter the required path each time you setup the pipeline):

* GATK,GATK,/path/to/GATK_jar_files
* GATK2,GATK2,/path/to/GATK2_jar_files
* picard,PICARD,/path/to/picard_jar_files
* R,R_LIB,/path/to/R_library_files
* bismark,BISMARK_GENOME_FOLDER,/path/to/bismark_genome_files
* trimmomatic,TRIMMOMATIC_JAR_PATH,/path/to/trimmomatic_jar_file


# COPYRIGHT & LICENSE

Copyright (c) 2011-2013 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.

The usage of a range of years within a copyright statement contained within this
distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads 'Copyright (c) 2005, 2007-
2009, 2011-2012' should be interpreted as being identical to a statement that
reads 'Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012' and a copyright
statement that reads "Copyright (c) 2005-2012' should be interpreted as being
identical to a statement that reads 'Copyright (c) 2005, 2006, 2007, 2008, 2009,
2010, 2011, 2012'."