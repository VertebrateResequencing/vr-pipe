use MooseX::Declare;

class VRPipe::Base::Declare::Syntax::Keyword::Role extends MooseX::Declare::Syntax::Keyword::Role {
    after add_namespace_customizations (Object $ctx, Str $package) {
        $ctx->add_preamble_code_parts(
            'use MooseX::StrictConstructor; use VRPipe::Base::Types ":all";'
        );
    }
}

1;