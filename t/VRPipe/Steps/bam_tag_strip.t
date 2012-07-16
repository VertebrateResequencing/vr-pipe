#!/usr/bin/env perl
use strict;
use warnings;
use Path::Class;

BEGIN {
    use Test::Most tests => 6;
    use VRPipeTest (required_env => [qw(VRPIPE_TEST_PIPELINES SAMTOOLS)]);
    use TestPipelines;
    
    use_ok('VRPipe::Steps::bam_strip_tags');
}

my ($output_dir, $pipeline, $step) = create_single_step_pipeline('bam_strip_tags', 'bam_files');
is_deeply [$step->id, $step->description], [1, 'Strips tags from bam files'], 'bam_strip_tags step created and has correct description';

my $tags_to_strip = 'OQ XM XG XO';
my @tags = split(/\s+/, $tags_to_strip);
my $tag_bam = VRPipe::File->create(path => file(qw(t data 2822_6.pe.bam))->absolute)->path;
my $strip_bam = VRPipe::File->create(path => file($output_dir, 'strip.bam'))->path;
my $ok = VRPipe::Steps::bam_strip_tags->tag_strip($tag_bam, $strip_bam, tags_to_strip => [@tags]);
is $ok, 1,  'tag_strip() ran okay';

my @expected_reads = ("IL3_2822:6:1:20:1902\t121\tfake_chr1\t36591\t37\t5M1I55M\t=\t36591\t0\tTGCCTAACATTGAAAATGTCAAGTAAAGACAGATGTGTCAACGTTACTTGACAGACCTGTA\t?>?<?>/445224;97556294/359543443595:79998649;7<=999?>?=B@==>>\tRG:Z:2822_6\tXT:A:U\tNM:i:3\tSM:i:37\tAM:i:0\tX0:i:1\tX1:i:0\tXM:i:2\tXO:i:1\tXG:i:1\tMD:Z:1T0T57",
"IL3_2822:6:1:20:1902\t181\tfake_chr1\t36591\t0\t53M1S\t=\t36591\t0\tTGGTGGGTCTAAATACACCCCCCGCCGGGTCTGAATATCAATACCATAGAAACA\t)-3;*0**%&.,.&%&,&(2;:32%*==656)2.==;<>\@9=0=@@@@@><9==\tRG:Z:2822_6\tXC:i:53",
"IL3_2822:6:1:24:477\t69\tfake_chr1\t70725\t0\t*\t=\t70725\t0\tTCAACATCTGTTTCTTCATCAAAACTGTATTGATTCCAAATGAAATTACGGATATGACCAA\t=99::==:=<:?:79;89:777;7169344434778356;/43476875351-7;@@=<:9\tRG:Z:2822_6",
"IL3_2822:6:1:24:477\t137\tfake_chr1\t70725\t37\t40M14S\t=\t70725\t0\tTTTCAGCCCGCGTTAAGAAATAGCCAGGGGAGGAATGAATCATGCGCCAGAAAA\t=>><@><;=2&246//</<8<==0&2=<<(7;8277092&2,00&*=0&,2(.&\tRG:Z:2822_6\tXC:i:40\tXT:A:U\tNM:i:2\tSM:i:37\tAM:i:0\tX0:i:1\tX1:i:0\tXM:i:2\tXO:i:0\tXG:i:0\tMD:Z:10A13A15",
"IL3_2822:6:1:21:402\t117\tfake_chr1\t104074\t0\t*\t=\t104074\t0\tTTCAAGTCCCCTTCGCACTAATCCCATTACCCTTTGTTGTATCAAGTGATTTTTCAACATC\t;:???\@81)'\$\$\$\$\$\$\$\$\$\$++'\$\$\$'\$\$\$'155338:9;89<=<<;;=>>>>??><@>>?\tRG:Z:2822_6",
"IL3_2822:6:1:21:402\t185\tfake_chr1\t104074\t37\t11S43M\t=\t104074\t0\tGCGGGCTTGATTGAGAAATTAGATGCTTGGGATGCTAAATATTCTGAAACATTA\t)%%%''+'.'0'295;3+2++<-7>::5)=66/=:35=>3;9>=:>>:;>=*/;\tRG:Z:2822_6\tXC:i:43\tXT:A:U\tNM:i:1\tSM:i:37\tAM:i:0\tX0:i:1\tX1:i:0\tXM:i:1\tXO:i:0\tXG:i:0\tMD:Z:17T25",
"IL3_2822:6:1:21:1773\t69\tfake_chr1\t170008\t0\t*\t=\t170008\t0\tTGCTGGCATGATGGCCCCCATTTCTTCTAGTTTTTATGGCAAGAAGACCCTGCTCAGATCA\t>>>\@A?@@@<<:A;:6=?6)6@>9==?'\$'<;??<<<?>:;64<<;9759<===;9>>@;;\tRG:Z:2822_6",
"IL3_2822:6:1:21:1773\t137\tfake_chr1\t170008\t37\t53M1S\t=\t170008\t0\tGACACCACCTCAGTTCCTGTACGGATGTCCACGCCATTTTCCAACATCTTCATC\t)9=1<=3:1*191&1;;/64=;=<0806==9467;1\%1,7<;1/:)0)/4=6,,\tRG:Z:2822_6\tXC:i:53\tXT:A:U\tNM:i:1\tSM:i:37\tAM:i:0\tX0:i:1\tX1:i:0\tXM:i:1\tXO:i:0\tXG:i:0\tMD:Z:47T5",
"IL3_2822:6:1:24:1620\t69\tfake_chr1\t179410\t0\t*\t=\t179410\t0\tGCGCTCCTTGGGTACCGAACACCAGACCAGATACGNGCAAAATGGCTAGACGAATGACCAA\t====:;;;98::/:998689::878686387575'\$'4459''\$21+'\$\$\$\$\$%=@>>>=<\tRG:Z:2822_6",
"IL3_2822:6:1:24:1620\t137\tfake_chr1\t179410\t37\t54M\t=\t179410\t0\tTGTGCAACTGTTAGCAGAGCATTTTGTACCAAGTGCTTCTTACTGGCTGGCAGG\t>>>==========:5====5======4:85=3;826=;=<<27:9933958;2*\tRG:Z:2822_6\tXT:A:U\tNM:i:1\tSM:i:37\tAM:i:0\tX0:i:1\tX1:i:0\tXM:i:1\tXO:i:0\tXG:i:0\tMD:Z:53T0");

my @actual_reads = get_bam_records($tag_bam);
is_deeply [@actual_reads[0..9]], \@expected_reads, 'reads is original bam as expected';

my @expected_stripped_reads = ("IL3_2822:6:1:20:1902\t121\tfake_chr1\t36591\t37\t5M1I55M\t=\t36591\t0\tTGCCTAACATTGAAAATGTCAAGTAAAGACAGATGTGTCAACGTTACTTGACAGACCTGTA\t?>?<?>/445224;97556294/359543443595:79998649;7<=999?>?=B@==>>\tRG:Z:2822_6\tXT:A:U\tNM:i:3\tSM:i:37\tAM:i:0\tX0:i:1\tX1:i:0\tMD:Z:1T0T57",
"IL3_2822:6:1:20:1902\t181\tfake_chr1\t36591\t0\t53M1S\t=\t36591\t0\tTGGTGGGTCTAAATACACCCCCCGCCGGGTCTGAATATCAATACCATAGAAACA\t)-3;*0**%&.,.&%&,&(2;:32%*==656)2.==;<>\@9=0=@@@@@><9==\tRG:Z:2822_6\tXC:i:53",
"IL3_2822:6:1:24:477\t69\tfake_chr1\t70725\t0\t*\t=\t70725\t0\tTCAACATCTGTTTCTTCATCAAAACTGTATTGATTCCAAATGAAATTACGGATATGACCAA\t=99::==:=<:?:79;89:777;7169344434778356;/43476875351-7;@@=<:9\tRG:Z:2822_6",
"IL3_2822:6:1:24:477\t137\tfake_chr1\t70725\t37\t40M14S\t=\t70725\t0\tTTTCAGCCCGCGTTAAGAAATAGCCAGGGGAGGAATGAATCATGCGCCAGAAAA\t=>><@><;=2&246//</<8<==0&2=<<(7;8277092&2,00&*=0&,2(.&\tRG:Z:2822_6\tXC:i:40\tXT:A:U\tNM:i:2\tSM:i:37\tAM:i:0\tX0:i:1\tX1:i:0\tMD:Z:10A13A15",
"IL3_2822:6:1:21:402\t117\tfake_chr1\t104074\t0\t*\t=\t104074\t0\tTTCAAGTCCCCTTCGCACTAATCCCATTACCCTTTGTTGTATCAAGTGATTTTTCAACATC\t;:???\@81)'\$\$\$\$\$\$\$\$\$\$++'\$\$\$'\$\$\$'155338:9;89<=<<;;=>>>>??><@>>?\tRG:Z:2822_6",
"IL3_2822:6:1:21:402\t185\tfake_chr1\t104074\t37\t11S43M\t=\t104074\t0\tGCGGGCTTGATTGAGAAATTAGATGCTTGGGATGCTAAATATTCTGAAACATTA\t)%%%''+'.'0'295;3+2++<-7>::5)=66/=:35=>3;9>=:>>:;>=*/;\tRG:Z:2822_6\tXC:i:43\tXT:A:U\tNM:i:1\tSM:i:37\tAM:i:0\tX0:i:1\tX1:i:0\tMD:Z:17T25",
"IL3_2822:6:1:21:1773\t69\tfake_chr1\t170008\t0\t*\t=\t170008\t0\tTGCTGGCATGATGGCCCCCATTTCTTCTAGTTTTTATGGCAAGAAGACCCTGCTCAGATCA\t>>>\@A?@@@<<:A;:6=?6)6@>9==?'\$'<;??<<<?>:;64<<;9759<===;9>>@;;\tRG:Z:2822_6",
"IL3_2822:6:1:21:1773\t137\tfake_chr1\t170008\t37\t53M1S\t=\t170008\t0\tGACACCACCTCAGTTCCTGTACGGATGTCCACGCCATTTTCCAACATCTTCATC\t)9=1<=3:1*191&1;;/64=;=<0806==9467;1\%1,7<;1/:)0)/4=6,,\tRG:Z:2822_6\tXC:i:53\tXT:A:U\tNM:i:1\tSM:i:37\tAM:i:0\tX0:i:1\tX1:i:0\tMD:Z:47T5",
"IL3_2822:6:1:24:1620\t69\tfake_chr1\t179410\t0\t*\t=\t179410\t0\tGCGCTCCTTGGGTACCGAACACCAGACCAGATACGNGCAAAATGGCTAGACGAATGACCAA\t====:;;;98::/:998689::878686387575'\$'4459''\$21+'\$\$\$\$\$%=@>>>=<\tRG:Z:2822_6",
"IL3_2822:6:1:24:1620\t137\tfake_chr1\t179410\t37\t54M\t=\t179410\t0\tTGTGCAACTGTTAGCAGAGCATTTTGTACCAAGTGCTTCTTACTGGCTGGCAGG\t>>>==========:5====5======4:85=3;826=;=<<27:9933958;2*\tRG:Z:2822_6\tXT:A:U\tNM:i:1\tSM:i:37\tAM:i:0\tX0:i:1\tX1:i:0\tMD:Z:53T0");

my @actual_stripped_reads = get_bam_records($strip_bam);
is_deeply [@actual_stripped_reads[0..9]], \@expected_stripped_reads, 'correct tags stripped from bam';

# test as part of a pipeline
my $setup = VRPipe::PipelineSetup->create(name => 'strip_setup',
                                       datasource => VRPipe::DataSource->create(type => 'fofn', method => 'all', source => file(qw(t data improvement_datasource.fofn))->absolute),
                                       output_root => $output_dir,
                                       pipeline => $pipeline,
                                       options => { bam_tags_to_strip => $tags_to_strip });

ok handle_pipeline(), 'single-step pipeline ran ok';

finish;
