=head1 NAME

VRPipe::FrontEnd - shared methods of all cmd-line frontend scripts

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

The various B<vrpipe-*> scripts use this module which provides methods that
they all need to function.

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011-2012 Genome Research Limited.

This file is part of VRPipe.

VRPipe is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see L<http://www.gnu.org/licenses/>.

=cut

use VRPipe::Base;

class VRPipe::FrontEnd {
    use Getopt::Long qw(GetOptions GetOptionsFromString);
    use Perl6::Form;
    use Module::Find;
    use VRPipe::Persistent::SchemaBase;
    
    has 'description' => (is => 'rw',
                          isa => 'Str',
                          required => 1);
    
    has 'opt_spec' => (is => 'rw',,
                       isa => 'ArrayRef[ArrayRef]',
                       lazy => 1,
                       builder => '_default_opt_spec');
    
    has '_opts_hash' => (is => 'ro',
                         isa => 'HashRef',
                         default => sub { {} },
                         writer => '_set_opts',
                         traits => ['Hash'],
                         handles => { opts => 'get',
                                      option_was_set => 'defined',
                                      '_set_opt' => 'set' });
    
    has 'usage' => (is => 'ro',
                    isa => 'Str',
                    writer => '_set_usage');
    
    has 'schema' => (is => 'ro',
                     isa => 'VRPipe::Persistent::Schema',
                     lazy => 1,
                     builder => '_build_schema');
    
    has '_multiple_setups' => (is => 'rw',
                               isa => 'Bool');
    
    has 'no_user_option' => (is => 'ro',
                             isa => 'Bool',
                             default => 0);
    
    method _default_opt_spec {
        return [ [ 'deployment=s', 'Use the production or testing database', { default => 'production' } ],
                 [ 'env|e=s', 'Use options stored in an environment variable' ],
                 [ 'help|h', 'Print this usage message and exit' ] ];
    }
    
    method _build_schema {
        my $m = VRPipe::Manager->get;
        return $m->result_source->schema;
    }
    
    sub BUILD {
        my $self = shift;
        
        # we initially used Getopt::Long::Descriptive, but because it does not
        # show what kind of arg each option takes, we implement here a similar
        # interface, but do not use the actual GLD code.
        
        my $opt_spec = $self->opt_spec;
        unless (@$opt_spec <= 4 && @$opt_spec >= 3 && $opt_spec->[-1]->[0] eq 'help|h') {
            my $default = $self->_default_opt_spec;
            
            foreach my $opt_spec (@$opt_spec) {
                my ($def, $help, $extra) = @$opt_spec;
                if ($help) {
                    my ($name, $req_or_opt, $type) = split(/([=:])/, $def);
                    if ($name eq 'setup') {
                        if ($type eq 's@') {
                            unless ($self->no_user_option) {
                                $default->[4] = $default->[2];
                                $default->[2] = [ 'user|u=s', 'Only show entries for PipelineSetups created by this user; use "all" to show entries for all users', { default => getlogin || getpwuid($<) || 'vrpipe' } ];
                                $default->[3] = [ 'deactivated', 'Also show deactivated PipelineSetups', { default => 0 } ];
                            }
                            $self->_multiple_setups(1);
                        }
                    }
                }
            }
            push(@$opt_spec, [], [ 'General options:' ], @$default);
        }
        
        my @opts;
        my %defaults;
        my %convert_to_persistent;
        my $script_name = file($0)->basename;
        my $usage = '';
        my (@shorts, $has_long);
        foreach my $opt_spec (@{$self->opt_spec}) {
            my ($def, $help, $extra) = @$opt_spec;
            if ($help) {
                push(@opts, $def);
                
                my ($name, $req_or_opt, $type) = split(/([=:])/, $def);
                
                my ($long, $short) = length($name) == 1 ? (undef, $name) : (split(/\|/, $name));
                if (! $short && length($long) == 1) {
                    $short = $long;
                    undef($long);
                }
                if ($short && length($short) != 1) {
                    undef($short);
                }
                my $option = '';
                if ($short) {
                    push(@shorts, $short);
                    $option .= ' -'.$short;
                }
                if ($long) {
                    $has_long = 1;
                    $option .= ' --'.$long;
                }
                
                if ($extra) {
                    my $default = $extra->{default};
                    if (defined $default) {
                        $defaults{$long || $short} = $default;
                        $help = "[$default] ".$help;
                    }
                    
                    my $persistent_class = $extra->{persistent_object};
                    if ($persistent_class) {
                        $convert_to_persistent{$long || $short} = $persistent_class;
                    }
                    
                    next if $extra->{hidden};
                }
                
                my $value = '';
                if ($type) {
                    if ($type =~ /^s/) {
                        $value = '<str>';
                    }
                    elsif ($type =~ /^i/) {
                        $value = '<int>';
                    }
                    elsif ($type =~ /^f/) {
                        $value = '<num>';
                    }
                    $value .= $req_or_opt eq ':' ? '?' : '';
                }
                
                $usage .= form "  {[[[[[[[[[[[[[[[} {IIII} {[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[}", $option, $value, $help;
            }
            elsif ($def) {
                $usage .= form "{[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[}", $def;
            }
            else {
                $usage .= "\n";
            }
        }
        my $o = '';
        if (@shorts) {
            $o .= ' [-'.join('', sort @shorts).']';
        }
        if ($has_long) {
            $o .= ' [long options...]';
        }
        $usage = $self->description."\n$script_name$o\n".$usage;
        $self->_set_usage($usage);
        
        my %opts;
        Getopt::Long::Configure("bundling");
        $self->help unless GetOptions(\%opts, @opts);
        $self->_set_opts(\%opts);
        
        my $from_env = $self->opts('env');
        if ($from_env) {
            $self->help unless GetOptionsFromString($ENV{$from_env}, \%opts, @opts);
            $self->_set_opts(\%opts);
        }
        
        while (my ($opt, $val) = each %defaults) {
            next if $self->option_was_set($opt);
            $self->_set_opt($opt => $val);
        }
        
        my $help = $self->opts('help');
        my $deployment = $self->opts('deployment');
        unless ($deployment eq 'production' || $deployment eq 'testing') {
            $self->error("--deployment must be <production|testing>");
            $help = 1;
        }
        $self->help if $help;
        
        VRPipe::Persistent::SchemaBase->database_deployment($deployment);
        require VRPipe::Persistent::Schema;
        
        while (my ($opt, $class) = each %convert_to_persistent) {
            next unless $self->option_was_set($opt);
            my $val = $self->opts($opt);
            my @desired = ref($val) eq 'ARRAY' ? @$val : ($val);
            
            my $schema = $self->schema;
            my @found;
            foreach my $desired (@desired) {
                my $found;
                if ($desired =~ /^\d+$/) {
                    $found = $schema->resultset($class)->find({ id => $desired });
                    unless ($found) {
                        $self->die_with_error("$desired is not a valid $class id");
                    }
                }
                else {
                    $found = $schema->resultset($class)->find({ name => $desired });
                    unless ($found) {
                        $self->die_with_error("$desired is not a valid $class name");
                    }
                }
                push(@found, $found);
            }
            
            if (ref($val) eq 'ARRAY') {
                $self->_set_opt($opt => \@found);
            }
            else {
                $self->_set_opt($opt => $found[0]);
                if (@found > 1) {
                    $self->error("--$opt @desired resulted in more than one $class object; only using the first");
                }
            }
        }
    }
    
    method help {
        $self->die_with_error($self->usage);
    }
    
    method get_pipelinesetups (Bool :$inactive?) {
        my @requested_setups = $self->option_was_set('setup') ? ($self->_multiple_setups ? @{$self->opts('setup')} : ($self->opts('setup'))) : ();
        
        my @setups;
        if (@requested_setups) {
            @setups = @requested_setups;
        }
        elsif ($self->_multiple_setups) {
            my $user = $self->opts('user');
            unless (defined $inactive) {
                $inactive = $self->opts('deactivated');
            }
            my $rs = $self->schema->resultset("PipelineSetup")->search( { $user eq 'all' ? () : (user => $user), $inactive ? () : (active => 1) } );
            while (my $setup = $rs->next) {
                push(@setups, $setup);
            }
        }
        
        if ($self->_multiple_setups && ! @setups) {
            $self->die_with_error("No PipelineSetups match your settings (did you remember to specifiy --user?)");
        }
        
        return $self->_multiple_setups ? @setups : $setups[0];
    }
    
    method output (@messages) {
        chomp($messages[-1]);
        print @messages, "\n";
    }
    method error (@messages) {
        chomp($messages[-1]);
        warn @messages, "\n";
    }
    
    method die_with_error (@messages) {
        $self->error(@messages);
        die "\n";
    }
    
    method display_hash (Str $name, HashRef $hash, ArrayRef[Str] $key_order?) {
        $key_order ||= [ sort { $a cmp $b } keys %$hash ];
        $self->output("$name:");
        my ($extra_tabs) = $name =~ /^(\t+)/;
        $extra_tabs ||= '';
        foreach my $key (@$key_order) {
            next unless defined $hash->{$key};
            $self->output("$extra_tabs\t", $key, ' => ', $hash->{$key});
        }
    }
    
    method sub_modules (Str $base) {
        $base = "VRPipe::$base";
        my @modules = findsubmod $base;
        unless (@modules) {
            @modules = findsubmod "${base}s";
        }
        my @names;
        foreach my $module (@modules) {
            my ($name) = $module =~ /::([^:]+)$/;
            push(@names, $name);
        }
        return @names;
    }
    
    method make_all_objects (Str $class) {
        my @modules = $self->sub_modules($class);
        $class = "VRPipe::$class";
        foreach my $name (@modules) {
            $class->get(name => $name);
        }
    }
    
    method ask_for_object (Str :$question!, Str :$class!, Str :$column!) {
        $self->make_all_objects($class);
        my $rs = $self->schema->resultset($class);
        my @things;
        while (my $thing = $rs->next) {
            push(@things, $thing);
        }
        my %things = map { $_->$column => $_ } @things;
        my @thing_keys = sort keys %things;
        $self->output("\n");
        my %num_to_key;
        foreach my $i (0..$#thing_keys) {
            my $num = $i + 1;
            my $key = $thing_keys[$i];
            $num_to_key{$num} = $key;
            my $output = "$num. $key";
            my $obj = $things{$key};
            if ($obj->can('description')) {
                $output .= ' ['.$obj->description.']';
            }
            $self->output($output);
        }
        my $chosen_num = $self->pick_number(question => $question, max => scalar(@thing_keys));
        return $things{$num_to_key{$chosen_num}};
    }
    
    method already_exists (Str $class!, Str $key!, Str $value!) {
        my @found = $self->schema->resultset($class)->search({ $key => $value });
        if (@found) {
            return "a $class already exists with $key '$value'";
        }
        return;
    }
    
    method ask_question (Str :$question!, ArrayRef :$possibles?, Str :$allow_multiple?, Str :$default?, Bool :$required?, CodeRef :$not_allowed?, ArrayRef :$na_args?) {
        undef $possibles unless $possibles && @$possibles;
        if (defined $default && length($default) == 0) {
            undef $default;
        }
        
        print STDERR "$question";
        my %allowed;
        if ($possibles) {
            print STDERR " <", join('|', @$possibles), ">";
            %allowed = map { $_ => 1 } @$possibles;
        }
        if (defined $default) {
            print STDERR " [$default]";
        }
        print STDERR ": ";
        
        my $answer;
        do {
            $answer = <STDIN>;
            chomp($answer);
            
            if ($possibles) {
                my $valid = 1;
                my @answers = $allow_multiple ? split(/$allow_multiple/, $answer) : ($answer);
                foreach my $sub_answer (@answers) {
                    unless (exists $allowed{$sub_answer}) {
                        $valid = 0;
                    }
                }
                
                unless ($valid) {
                    if (defined $default && ! $answer) {
                        $answer = $default;
                    }
                    else {
                        undef $answer;
                        my $one_of = $allow_multiple ? 'one or more of' : 'one of';
                        print STDERR "Your answer can only consist of $one_of <", join('|', @$possibles), ">: ";
                    }
                }
            }
            elsif (! $answer && "$answer" ne "0") {
                if (defined $default) {
                    $answer = $default;
                }
                elsif ($required) {
                    undef $answer;
                    print STDERR "An answer is required. Try again: ";
                }
            }
            
            if ($answer && $not_allowed) {
                push(@$na_args, $answer);
                my $reason = &$not_allowed(@$na_args);
                if ($reason) {
                    undef $answer;
                    pop(@$na_args);
                    print STDERR "Your answer isn't allowed because $reason. Try again: ";
                }
            }
        }
        while (! defined $answer);
        
        return $answer;
    }
    
    method pick_number (Str :$question!, PositiveInt :$max!, PositiveInt :$default?) {
        return $self->ask_question(question => $question, possibles => [1..$max], required => 1, $default ? (default => $default) : ());
    }
}

1;
