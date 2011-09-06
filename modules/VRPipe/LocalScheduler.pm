use VRPipe::Base;

class VRPipe::LocalScheduler {
    with 'MooseX::Daemonize';
    use VRPipe::Config;
    my $vrp_config = VRPipe::Config->new();
    use VRPipe::Persistent::Schema;
    
    has 'deployment' => (is => 'ro',
                         isa => 'Str',
                         default => 'testing',
                         documentation => 'testing|production (default testing): determins which database is used to track locally scheduled jobs');
    
    has 'pidbase' => (is => 'rw',
                      isa => Dir,
                      coerce => 1,
                      builder => '_default_pidbase',
                      lazy => 1);
    
    # for submit
    has 'o' => (is => 'ro',
                isa => 'Str',
                documentation => 'absolute path to stdout of scheduler (required for submit)');
    
    has 'e' => (is => 'ro',
                isa => 'Str',
                documentation => 'absolute path to stderr of scheduler (required for submit)');
    
    has 'a' => (is => 'ro',
                isa => 'Int',
                documentation => 'to specify a job array give the size of the array (default 1)',
                default => 1);
    
    
    sub BUILD {
        my $self = shift;
        my $d = $self->deployment;
        unless ($d eq 'testing' || $d eq 'production') {
            $self->throw("'$d' is not a valid deployment type; --deployment testing|production");
        }
        VRPipe::Persistent::SchemaBase->database_deployment($self->deployment);
    }
    
    method _default_pidbase {
        my $method_name = VRPipe::Persistent::SchemaBase->database_deployment.'_scheduler_output_root';
        return $vrp_config->$method_name();
    }
    
    after start {
        return unless $self->is_daemon;
        # your daemon code here ...
    }
    
    method submit (Str $cmd) {
        my $o_file = $self->o;
        my $e_file = $self->e;
        unless ($o_file && $e_file) {
            $self->throw("-o and -e are required for submit");
        }
        
        my $array_size = $self->a;
        my $lsj = VRPipe::LocalSchedulerJob->get(cmd => $cmd, array_size => $array_size);
        
        foreach my $aid (1..$array_size) {
            my $this_o = $o_file;
            $this_o =~ s/\%I/$aid/g;
            my $this_e = $e_file;
            $this_e =~ s/\%I/$aid/g;
            VRPipe::LocalSchedulerJobState->get(localschedulerjob => $lsj, aid => $aid, o_file => $this_o, e_file => $this_e);
        }
        
        my $sid = $lsj->id;
        print "Job <$sid> is submitted\n";
    }
    
    method job (Str $id) {
        my @lsjss;
        if ($id =~ /(\d+)\[(\d+)\]/) {
            
        }
        elsif ($id =~ /^\d+$/) {
            
        }
        else {
            $self->throw("bad id format '$id'");
        }
    }
}