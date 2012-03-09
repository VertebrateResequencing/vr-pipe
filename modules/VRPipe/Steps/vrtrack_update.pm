use VRPipe::Base;

class VRPipe::Steps::vrtrack_update with VRPipe::StepRole {
    # eval these so that test suite can pass syntax check on this module when
    # VertRes is not installed
    eval "use VertRes::Utils::VRTrackFactory;";
    
    method options_definition {
        return { vrtrack_db => VRPipe::StepOption->get(description => 'the name of your VRTrack database (other connection settings are taken from the standard VRTrack environment variables)') };
    }
    method inputs_definition {
        return { };
    }
    method body_sub {
        return sub { };
    }
    method outputs_definition {
    	return { };
    }
    method post_process_sub {
        return sub { return 1; };
    }
    method description {
        return "Empty shell step for basing other vrtrack-related steps on";
    }
    method max_simultaneous {
        return 15;
    }
    
    method get_vrtrack (ClassName|Object $self: Str :$db!, Str :$mode = 'rw') {
	return VertRes::Utils::VRTrackFactory->instantiate(database => $db, mode => $mode) || $self->throw("Could not connect to the database '$db'");
    }
}

1;
