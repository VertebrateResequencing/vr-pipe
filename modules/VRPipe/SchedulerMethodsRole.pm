use VRPipe::Base;

role VRPipe::SchedulerMethodsRole {
    requires 'start_command';
    requires 'stop_command';
    requires 'submit_command';
    requires 'submit_args';
    requires 'determine_queue';
    requires 'get_1based_index';
    requires 'get_sid';
    requires 'kill_sid';
    requires 'sid_status';
}

1;