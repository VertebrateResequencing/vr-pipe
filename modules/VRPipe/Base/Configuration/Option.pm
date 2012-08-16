
=head1 NAME

VRPipe::Base::Configuration::Option - configuration options

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011 Genome Research Limited.

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

class VRPipe::Base::Configuration::Option {
    use VRPipe::Base::Configuration::Env;
    
    has _attr => (
        is       => 'ro',
        does     => 'VRPipe::Base::Configuration::Trait::Attribute::ConfigKey',
        required => 1
    );
    
    has _obj => (
        is       => 'ro',
        does     => 'VRPipe::Base::Configuration::Trait::Object',
        weak_ref => 1,
        required => 1
    );
    
    method key {
        $self->_attr->name;
    }
    
    method value (StrOrEnv $set?) {
        my $method_name = $self->key;
        
        if (defined $set) {
            unless (ref $set) {
                if ($set eq 'undef') {
                    return $self->_obj->$method_name(undef);
                }
                
                my $valid = $self->valid;
                if ($valid) {
                    my %allowed = map { $_ => 1 } @{$valid};
                    if (!exists $allowed{$set}) {
                        $self->throw("'$set' was not a valid value for option '$method_name'");
                    }
                }
            }
            
            return $self->_obj->$method_name($set);
        }
        
        $self->_obj->$method_name();
    }
    
    method env (Str $var) {
        my $env = VRPipe::Base::Configuration::Env->new(variable => $var);
        $self->value($env);
        return $env;
    }
    
    method question {
        my $attr = $self->_attr;
        if ($attr->has_question) {
            return $attr->question;
        }
        else {
            my $key = $self->key;
            return "$key?";
        }
    }
    
    method question_number {
        $self->_attr->question_number;
    }
    
    method valid {
        my $attr = $self->_attr;
        return $attr->valid if $attr->has_valid;
        return;
    }
    
    method secure {
        my $attr = $self->_attr;
        return $attr->secure if $attr->has_secure;
        return;
    }
    
    # prompt-related methods stolen from Module::Build
    method _is_interactive {
        return -t STDIN && (-t STDOUT || !(-f STDOUT || -c STDOUT));
    }
    
    method _is_unattended {
        return $ENV{PERL_MM_USE_DEFAULT} || (!$self->_is_interactive && eof STDIN);
    }
    
    method _readline {
        return undef if $self->_is_unattended;
        
        # *** could use http://search.cpan.org/~phoenix/Term-ReadPassword/ReadPassword.pm
        # for secure options so that what the user types is replaced on the
        # terminal with *s. Probably not worth it though...
        
        my $answer = <STDIN>;
        chomp $answer if defined $answer;
        return $answer;
    }
    
    method _prompt (Str $valid, Str $default) {
        if ($valid) {
            $valid .= ' ';
        }
        
        local $| = 1;
        print $self->question, ' ', $valid, '[', $default, ']', ' ';
        
        if ($self->_is_unattended && !$default) {
            die <<EOF;
ERROR: This build seems to be unattended, but there is no default value
for this question. Aborting.
EOF
        }
        
        my $ans = $self->_readline();
        if (
            !defined($ans) # Ctrl-D or unattended
            or !length($ans)
          ) {              # User hit return
            $ans = $default;
        }
        return $ans;
    }
    
    method prompt {
        my $valid_ref = $self->valid;
        my $valid     = '';
        if ($valid_ref) {
            $valid = '<' . join('|', @{$valid_ref}) . '>';
        }
        
        my $default = $self->value;
        unless (defined $default) {
            $default = '';
        }
        if (ref $default) {
            my $env_value = $default->value ? "'$default'" : 'undefined';
            $default = 'ENV{' . $default->variable . '} (currently ' . $env_value . ')';
        }
        
        my $answer = $self->_prompt($valid, $default);
        
        if ($valid_ref) {
            my %allowed = map { $_ => 1 } @{$valid_ref};
            my $env_answer;
            if ($answer =~ /^ENV\{(.+)}/) {
                $env_answer = $ENV{$1};
            }
            while (!exists $allowed{$answer} && (defined $env_answer ? (!exists $allowed{$env_answer}) : 1)) {
                warn "'$answer' was not a valid answer for that question; try again:\n";
                $answer = $self->_prompt($valid, $default);
                if ($answer =~ /^ENV\{(.+)}/) {
                    $env_answer = $ENV{$1};
                }
            }
        }
        
        if ($answer =~ /^ENV\{(.+)}/) {
            my $var = $1;
            $self->env($var);
        }
        else {
            $self->value($answer);
        }
    }
}

1;
