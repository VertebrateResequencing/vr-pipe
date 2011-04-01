use VRPipe::Base;

role VRPipe::Base::Configuration::Trait::Object {
    use Data::Dumper;
    use Module::Find;
    
    has config_module => (
        is      => 'ro',
        isa     => File,
        coerce  => 1,
        lazy    => 1,
        builder => '_build_config_module'
    );
    
    has _raw_config => (
        is      => 'ro',
        isa     => 'HashRef',
        lazy    => 1,
        builder => '_build_raw_config'
    );
    
    method _build_config_module {
        return file(qw(modules VRPipe SiteConfig.pm));
    }
    
    method _build_raw_config {
        # given the path to the config module, load the corresponding module,
        # which we require to be in the VRPipe namespace
        my $config_module = $self->config_module;
        my @dirs = $config_module->dir->dir_list;
        my @parents;
        my $saw_vrpipe = 0;
        my @module_name;
        foreach my $dir (@dirs) {
            if ($dir eq 'VRPipe') {
                $saw_vrpipe = 1;
            }
            
            if ($saw_vrpipe) {
                push(@module_name, $dir);
            }
            else {
                push(@parents, $dir);
            }
        }
        $self->throw("'$config_module' did not correspond to a VRPipe module") unless $saw_vrpipe;
        
        my $module_basename = $config_module->basename;
        $module_basename =~ s/\.pm$//; 
        my $module_name = join('::', @module_name, $module_basename);
        
        try { local @INC = (dir(@parents)->stringify); eval "require $module_name"; die $@ if $@; } # require $module_name; with no eval does not work
        catch { return {} }
        return $module_name->get_config() || {};
    }
    
    method _from_config (Str $name) {
        my $hash = $self->_raw_config();
        return $hash->{$name};
    }
    
    method write_config_module (HashRef :$values) {
        my $config_module = $self->config_module();
        
        $self->throw('Cannot write a configuration module without config_module set') unless defined $config_module;
        
        my $dd = Data::Dumper->new( [ $values ], [ 'site_config' ] );
        open( my $fh, '>', $config_module ) or $self->throw("Cannot write to '$config_module': $!");
        printf $fh <<'END_HERE', $dd->Dump();
package VRPipe::SiteConfig;
use strict;
use warnings;

my $site_config;
%s;

sub get_config {
    return $site_config;
}

1;
END_HERE
        close($fh);
    }
    
    method get_options ($class:) {
        my $meta = $class->meta;
        my @names;
        foreach my $attr ($meta->get_all_attributes) {
            my $name = $attr->name;
            push(@names, $name);
        }
        return @names;
    }
}

1;