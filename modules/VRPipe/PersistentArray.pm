use VRPipe::Base;

class VRPipe::PersistentArray extends VRPipe::Persistent {
    has 'members' => (is => 'rw',
                      isa => ArrayRefOfPersistent,
                      lazy => 1,
                      builder => '_build_members',
                      predicate => '_members_populated');
    
    __PACKAGE__->make_persistent(has_many => [_members => 'VRPipe::PersistentArrayMember']);
    
    around get (ClassName|Object $self: ArrayRefOfPersistent :$members?, Persistent :$id?) {
        $self->throw("both id and members cannot be supplied to get() at the same time") if $id && $members;
        $self->throw("get() needs id or members option") unless $id || $members;
        
        if ($id) {
            # get by id
            return $self->$orig(id => $id);
        }
        else {
            # create a new row, then use the new id to create new
            # PersistentArrayMember rows for each supplied member
            my $array = $self->$orig();
            
            my $index = 0;
            foreach my $member (@{$members}) {
                VRPipe::PersistentArrayMember->get(persistentarray => $array, class => ref($member), class_id => $member->id, array_index => ++$index);
            }
            
            $array->members($members);
            
            return $array;
        }
    }
    
    method _build_members {
        # get all PersistentArrayMember rows with our id, and instantiate the
        # corresponding Persistent objects
        my @members;
        foreach my $array_member ($self->_members) {
            push(@members, $self->_array_member_to_member($array_member));
        }
        
        return \@members;
    }
    
    method _array_member_to_member (VRPipe::PersistentArrayMember $array_member) {
        my $class = $array_member->class;
        return $class->get(id => $array_member->class_id);
    }
    
    method member (PositiveInt $index) {
        if ($index == 0) {
            $self->throw("The array index supplied to member() is 1-based");
        }
        
        if ($self->_members_populated) {
            my $members = $self->members;
            return $members->[$index - 1];
        }
        else {
            foreach my $array_member ($self->_members) {
                if ($array_member->array_index == $index) {
                    return $self->_array_member_to_member($array_member);
                }
            }
        }
    }
    
    method size {
        return scalar $self->_members;
    }
}

1;