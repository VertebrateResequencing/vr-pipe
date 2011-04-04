use VRPipe::Base;

role VRPipe::Base::Configuration::Trait::Object {
    use Data::Dumper;
    
    has config_module => (
        is      => 'ro',
        isa     => 'Str',
        lazy    => 1,
        builder => '_build_config_module'
    );
    
    has config_module_path => (
        is      => 'ro',
        isa     => File,
        coerce  => 1,
        lazy    => 1,
        builder => '_build_config_path'
    );
    
    has _raw_config => (
        is      => 'ro',
        isa     => 'HashRef',
        lazy    => 1,
        builder => '_build_raw_config'
    );
    
    method _build_config_module {
        return 'VRPipe::SiteConfig';
    }
    
    method _build_config_path {
        my $module = $self->config_module;
        
        my $path = file(split('::', $module));
        $path .= '.pm';
        
        eval "require $module";
        if (exists $INC{$path}) {
            $path = $INC{$path};
        }
        else {
            # stick it in the first writable path in INC
            my $writable_inc;
            foreach my $inc (@INC) {
                if (-w $inc) { # only a permission test; doesn't really guarantee we can really write there
                    $writable_inc = $inc;
                    last;
                }
            }
            $writable_inc || $self->throw("There don't seem to be any writable paths in your INC");
            
            $path = file($writable_inc, $path);
        }
        
        return $path;
    }
    
    method _build_raw_config {
        my $config_module = $self->config_module;
        try { eval "require $config_module"; die $@ if $@; } # require $config_module; with no eval does not work
        catch { return {} }
        return $config_module->get_config() || {};
    }
    
    method _from_config (Str $name) {
        my $hash = $self->_raw_config();
        return $hash->{$name};
    }
    
    method write_config_module (HashRef :$values) {
        my $config_module = $self->config_module_path();
        
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