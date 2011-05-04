use VRPipe::Base;

class VRPipe::PersistentArray extends VRPipe::Persistent {
    has 'id' => (is => 'rw',
                 isa => IntSQL[16],
                 traits => ['VRPipe::Persistent::Attributes'],
                 is_auto_increment => 1,
                 is_primary_key => 1);
    
    has 'members' => (is => 'rw',
                      isa => ArrayRefOfPersistent,
                      lazy => 1,
                      builder => '_build_members');
    
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
            
            foreach my $member (@{$members}) {
                VRPipe::PersistentArrayMember->get(persistentarray => $array, class => ref($member), class_id => $member->id);
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
            my $class = $array_member->class;
            push(@members, $class->get(id => $array_member->class_id));
        }
        
        return \@members;
    }
}

1;