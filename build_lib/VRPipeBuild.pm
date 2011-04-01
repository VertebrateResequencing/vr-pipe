package VRPipeBuild;

BEGIN { unshift(@INC, './modules') }
use base 'Module::Build';
use VRPipe::Config;

my $vrp_config = VRPipe::Config->new();
my $siteconfig_module = $vrp_config->config_module;

sub create_site_config {
    my $self = shift;
    
    my $do_config = 'y';
    if (-s $siteconfig_module) {
        $do_config = $self->prompt('SiteConfig module already exists; do you wish to go through setup again? [y/n] ', 'n');
    }
    
    if ($do_config eq 'y') {
        my @options = $vrp_config->get_options;
        warn "options: (@options)\n";
        
        my %config;
        
        $config{database_mysql_dbname} = $self->prompt( 'test DB name: ', $vrp_config->database_mysql_dbname);
        $config{database_mysql_username} = $self->prompt( 'test DB user: ', $vrp_config->database_mysql_username);
        
        $vrp_config->write_config_module(values => \%config);
    }
}

sub ACTION_realclean {
    my $self = shift;
    unlink($siteconfig_module);
    $self->SUPER::ACTION_realclean( @_ );
}

1;

