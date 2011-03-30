use MooseX::Declare;

class VRPipe::Base extends MooseX::Declare is dirty {
    use aliased 'VRPipe::Base::Declare::Syntax::Keyword::Class', 'ClassKeyword';
    use aliased 'VRPipe::Base::Declare::Syntax::Keyword::Role',  'RoleKeyword';
    use aliased 'MooseX::Declare::Syntax::Keyword::Namespace',   'NamespaceKeyword';
    
    clean;
    
    around keywords (ClassName $self:) {
        ClassKeyword->new(identifier => 'class'),
        RoleKeyword->new(identifier => 'role'),
        NamespaceKeyword->new(identifier => 'namespace'),
    }
}

1;