use VRPipe::Base;

class VRPipe::Utils::picard extends VRPipe::Utils::java {
    
    has 'picard_path' => (is => 'ro',
                          isa => Dir,
                          coerce => 1,
                          default => "$ENV{PICARD}");
    
    method determine_picard_version (ClassName|Object $self:) {
        my $picard_path = $self->picard_path;
        my $version = 0;
        opendir(my $dh, $picard_path) || $self->throw("Could not open picard directory $picard_path");
        foreach (readdir $dh) 
        {
            if (/^picard-([\d\.]+)\.jar/)
            {
                $version = $1;
                last;
            }
        }
        closedir($dh);
        return $version;
    }
}

1;
