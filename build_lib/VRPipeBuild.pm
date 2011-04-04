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
        $do_config = $self->y_n('SiteConfig module already exists; do you wish to go through setup again? [y/n]', 'n');
    }
    
    if ($do_config) {
        my @options = $vrp_config->get_options;
        warn "options: (@options)\n";
        
        my %config;
        
        $config{database_mysql_dbname} = $self->prompt( 'test DB name: ', $vrp_config->database_mysql_dbname);
        $config{database_mysql_username} = $self->prompt( 'test DB user: ', $vrp_config->database_mysql_username);
        
        $vrp_config->write_config_module(values => \%config);
        warn "Wrote configuration to $siteconfig_module_path\n";
    }
}

sub ACTION_realclean {
    my $self = shift;
    unlink($siteconfig_module_path);
    $self->SUPER::ACTION_realclean( @_ );
}

1;

