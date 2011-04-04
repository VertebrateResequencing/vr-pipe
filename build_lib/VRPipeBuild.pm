package VRPipeBuild;

BEGIN { unshift(@INC, './modules') }
use base 'Module::Build';
use VRPipe::Config;
use VRPipe::Base::Configuration::Env;

my $vrp_config = VRPipe::Config->new();
my $siteconfig_module_path = $vrp_config->config_module_path;

sub create_site_config {
    my $self = shift;
    
    my $do_config = 'y';
    if (-s $siteconfig_module_path) {
        $do_config = $self->y_n('SiteConfig module already exists; do you wish to go through setup again? [y|n]', 'n');
    }
    
    if ($do_config) {
        print "\nWhen answering the following questions, hit return to accept the default given in [].
To say that you don't want the default or any other value, type the word 'undef' without quotes.
To specify the answer should come from an environment variable, type 'ENV{variable_name}' without the quotes, replacing 'variable_name' as desired\n\n";
        
        my @options = $vrp_config->get_options;
        
        foreach my $option (@options) {
            my $key = $option->{key};
            
            my $valid = '';
            if (defined $option->{valid}) {
                $valid = '['.join('|', @{$option->{valid}}).']';
            }
            
            my $default = $vrp_config->$key();
            if (ref $default) {
                my $env_value = $default->value ? "'$default'" : 'undefined';
                $default = 'ENV{'.$default->variable.'} (currently '.$env_value.')';
            }
            
            my $answer = $self->prompt($option->{question}.' '.$valid, $default);
            
            if (defined $option->{valid}) {
                my %allowed = map { $_ => 1 } @{$option->{valid}};
                while (! exists $allowed{$answer}) {
                    warn "'$answer' was not a valid answer for that question; try again:\n";
                    $answer = $self->prompt($option->{question}.' '.$valid, $default);
                }
            }
            
            if ($answer =~ /^ENV\{(.+)}/) {
                $answer = VRPipe::Base::Configuration::Env->new(variable => $1);
            }
            
            $vrp_config->$key($answer);
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

