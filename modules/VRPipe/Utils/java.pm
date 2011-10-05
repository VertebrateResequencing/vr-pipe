use VRPipe::Base;

class VRPipe::Utils::java {
    
    has 'java_exe' => (is => 'ro',
                       isa => 'Str',
                       default => 'java');
    
    method jvm_args (ClassName|Object $self: Int $memory, Str|Dir $dir?) {
        my $java_mem = int($memory * 0.9);
        my $xss = 280;
        if ($java_mem > 1000)
        {
            $xss = " -Xss${xss}m";
        }
        else
        {
            $xss = ''; # login node with small memory limit doesn't like Xss option at all
        }
        
        my $temp_dir = '';
        if ($dir)
        {
            $temp_dir = ' -Djava.io.tmpdir='.$self->tempdir(DIR => $dir);
        }
        return qq[-Xmx${java_mem}m -Xms${java_mem}m$xss$temp_dir -server -XX:+UseParallelGC -XX:ParallelGCThreads=2];
    }
}

1;
