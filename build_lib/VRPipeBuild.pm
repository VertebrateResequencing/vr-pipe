package VRPipeBuild;

BEGIN { unshift(@INC, './modules') }
use base 'Module::Build';
use VRPipe::Config;

my $vrp_config = VRPipe::Config->new();
my $siteconfig_module_path = $vrp_config->config_module_path;

sub create_site_config {
    my $self = shift;
    
    my $do_config = 'y';
    if (-s $siteconfig_module_path) {
        $do_config = $self->y_n('SiteConfig module already exists; do you wish to go through setup again? [y|n]', 'n');
    }
    
    if ($do_config) {
        print "
When answering the following questions, hit return to accept the default given in [].
To say that you don't want the default or any other value, type the word 'undef' without quotes.
To specify the answer should come from an environment variable, type 'ENV{variable_name}' without the quotes, replacing 'variable_name' as desired\n\n";
        
        while (my $option = $vrp_config->next_option) {
            $option->prompt;
        }
        
        $vrp_config->write_config_module();
        warn "Wrote configuration to $siteconfig_module_path\n";
    }
}

sub ACTION_realclean {
    my $self = shift;
    unlink($siteconfig_module_path);
    $self->SUPER::ACTION_realclean( @_ );
}

1;

