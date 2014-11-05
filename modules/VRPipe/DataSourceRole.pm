
=head1 NAME

VRPipe::DataSourceRole - a role that must be used by all DataSources

=head1 SYNOPSIS

*** more documentation to come

=head1 DESCRIPTION

*** more documentation to come

=head1 AUTHOR

Sendu Bala <sb10@sanger.ac.uk>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011, 2013 Genome Research Limited.

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

role VRPipe::DataSourceRole {
    has 'method' => (
        is  => 'ro',
        isa => 'Str'
    );
    
    has 'source' => (
        is  => 'ro',
        isa => 'Defined'
    );
    
    has 'options' => (
        is  => 'ro',
        isa => 'HashRef'
    );
    
    has '_handle' => (
        is  => 'rw',
        isa => 'Defined',
    );
    
    has '_changed_marker' => (
        is  => 'rw',
        isa => 'Str'
    );
    
    has '_datasource_id' => (
        is  => 'ro',
        isa => Persistent
    );
    
    requires '_open_source';
    requires '_has_changed';
    requires '_update_changed_marker';
    requires 'description';
    requires 'source_description';
    requires 'method_description';
    
    method _generate_elements {
        # get a fresh handle every time
        my $handle = $self->_open_source || return;
        $self->_handle($handle);
        
        # call the datasource method
        my $method = $self->method;
        $self->can($method) || $self->throw("Invalid method '$method' for " . ref($self));
        $self->$method(%{ $self->options }, handle => $handle);
    }
    
    method get_methods {
        my $metarole     = Moose::Meta::Role->initialize(__PACKAGE__);
        my @meta_methods = $metarole->get_method_list;
        push(@meta_methods, $metarole->get_attribute_list);
        push(@meta_methods, $metarole->get_required_method_list);
        my %role_methods = map { $_ => 1 } @meta_methods;
        
        my %methods;
        my $this_classmeta = $self->meta;
        foreach my $class ($this_classmeta->linearized_isa) {
            next unless ($class =~ /VRPipe::DataSource/);
            my $classmeta = $class->meta;
            my %self_methods = map { $_->name => 1 } $classmeta->get_all_attributes; #*** really we want to skip methods on all of the roles a class uses, but I don't know how to get that list of all roles...
            $self_methods{new}     = 1;
            $self_methods{DESTROY} = 1; #*** bleugh, is there a way to detect these kinds of things automatically as well?
            map { $methods{ $_->name } = 1 } grep { $_->original_package_name eq $class && !exists $role_methods{ $_->name } && !exists $self_methods{ $_->name } && index($_->name, '_') != 0 } $classmeta->get_all_methods;
        }
        return keys %methods;
    }
    
    method method_options (Str $method) {
        $self->throw("'$method' is not a method of " . ref($self)) unless $self->can($method);
        my $method_meta = $self->meta->find_method_by_name($method);
        $self->throw("method $method in class " . ref($self) . " is not a signature method") unless $method_meta->isa('MooseX::Method::Signatures::Meta::Method');
        my $sig = $method_meta->parsed_signature;
        
        my @return;
        foreach my $kind (qw(positional named)) {
            my $check_method = 'has_' . $kind . '_params';
            next unless $sig->$check_method;
            # *** actually, we don't support positional args in method options?
            next if $kind eq 'positional';
            my $sig_method = $kind . '_params';
            foreach my $param ($sig->$sig_method) {
                my $var = $param->variable_name;
                next if $var eq '$handle';
                $var =~ s/^\$//; #*** we only handle $ sigils atm...
                my $tc_name = $param->meta_type_constraint->name;
                $tc_name =~ s/VRPipe::Base::Types:://g;
                push(@return, [$kind, $var, $param->required, $param->default_value, $tc_name]);
            }
        }
        
        return @return;
    }
    
    # a datasource method should call this when it wants to create/update
    # all the current dataelements for the source. Here is where we handle
    # withdrawing things no longer in the source as well.
    method _create_elements (ArrayRef $e_args) {
        # first make any existing elements withdrawn
        VRPipe::DataElement->search_rs({ datasource => $self->_datasource_id })->update({ withdrawn => 1 });
        
        # now create/update, setting withdrawn => 0 if not set
        VRPipe::DataElement->bulk_create_or_update(map { $_->{withdrawn} = 0 unless defined $_->{withdrawn}; $_; } @$e_args);
    }
    
    # $anti_repeat_store is just an empty hash ref that you must supply; define
    # it before a loop that calls this method
    method _start_over_elements_due_to_file_metadata_change (HashRef $result, ArrayRef $changed_details, HashRef $anti_repeat_store) {
        # 'changed' is based on file metadata changing, and the file may
        # have had its metadata applied in some other pipeline for some
        # other datasource. We'll start_from_scratch all affected
        # dataelementstates for all datasources, then update the metadata
        
        # get all the files with changed metadata
        my %affected_file_ids;
        foreach my $fm (@{ $result->{changed} }) {
            my ($vrfile) = @$fm;
            $affected_file_ids{ $vrfile->id } = 1;
        }
        my @affected_file_ids = keys %affected_file_ids;
        
        # get all the dataelements that have those files as inputs
        my %flids;
        foreach my $fid (@affected_file_ids) {
            foreach my $flid (VRPipe::FileListMember->get_column_values('filelist', { file => $fid })) {
                $flids{$flid} = 1;
            }
        }
        
        my %elements;
        foreach my $flid (keys %flids) {
            foreach my $element (VRPipe::DataElement->search({ filelist => $flid })) {
                $elements{ $element->id } = $element;
            }
        }
        
        # reset element states first
        my $class = ref($self);
        my ($type) = $class =~ /VRPipe::DataSource::(\S+)/;
        $changed_details = join(' | ', @$changed_details);
        foreach my $element (values %elements) {
            foreach my $estate ($element->element_states) {
                $estate->pipelinesetup->log_event("$type DataSource will call start_from_scratch because file metadata changed [$changed_details]", dataelement => $estate->dataelement->id);
                #*** shouldn't we work out exactly which steps used the input
                # file and along the chain?
                $estate->start_from_scratch(anti_repeat_store => $anti_repeat_store);
            }
        }
        
        #*** problems happen if we start_from_scratch some of them, but then get
        # killed before updating the metadata... also, what happens if something
        # starts running again following a start_from_scratch before we updated
        # the metadata?... *** should probably lock up the estates with timed
        # redis locks, then release them after we do the metadata update...
        
        # then change metadata in files
        foreach my $fm (@{ $result->{changed} }) {
            my ($vrfile, $new_metadata) = @$fm;
            $vrfile->add_metadata($new_metadata, replace_data => 1);
        }
    }
    
    method _vals_different ($orig, $new) {
        if (!ref($orig) && !ref($new)) {
            if ($orig ne $new) {
                return "$orig => $new";
            }
            return;
        }
        
        if (ref($orig) && !ref($new)) {
            return '[' . join(', ', @$orig) . "] => $new";
        }
        elsif (!ref($orig) && ref($new)) {
            return "$orig => " . '[' . join(', ', @$new) . ']';
        }
        
        my %orig = map { $_ => 1 } @$orig;
        my %new  = map { $_ => 1 } @$new;
        foreach my $orig (keys %orig) {
            unless (delete $new{$orig}) {
                return "$orig => undef";
            }
        }
        foreach my $new (keys %new) {
            return "undef => $new";
        }
        
        return;
    }
}

1;
