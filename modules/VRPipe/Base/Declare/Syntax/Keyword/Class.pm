use MooseX::Declare;

class VRPipe::Base::Declare::Syntax::Keyword::Class extends MooseX::Declare::Syntax::Keyword::Class {
    around import_symbols_from (Object $ctx) {
        'VRPipe::Base::UseMoose';
    }

    before add_namespace_customizations (Object $ctx, Str $package) {
        $ctx->add_preamble_code_parts(
            'use metaclass (metaclass => "Moose::Meta::Class");' #, error_class => "MooseX::Error::Exception::Class"  ## MooseX::Error::Exception::Class does not actually compile
        );
    }
    
    after add_namespace_customizations (Object $ctx, Str $package) {
        $ctx->add_preamble_code_parts(
            'use MooseX::StrictConstructor; use VRPipe::Base::Types ":all";'
        );
    }
}

1;