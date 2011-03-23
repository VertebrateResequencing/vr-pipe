=head1 DESCRIPTION

Just testing out DBIx::Class

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

use MooseX::Declare;

class VRPipe::CD extends VRPipe::Persistent {
    use VRPipe::Base::Types qw(Varchar IntSQL);
    
    has 'cdid' => (is => 'rw',
                   isa => IntSQL[16],
                   lazy_build => 1,
                   metaclass => 'MooseX::MetaDescription::Meta::Attribute',
                   description => {
                       is_auto_increment => 1,
                       is_primary_key => 1
                   });
    
    has 'artist' => (is => 'rw',
                     lazy_build => 1,
                     isa => IntSQL[16]);
    
    has 'title' => (is => 'rw',
                    lazy_build => 1,
                    isa => Varchar[256]);
    
    has 'rank' => (is => 'rw',
                   lazy_build => 1,
                   isa => IntSQL[16]);
    
    has 'non_column_data' => (is => 'rw',
                              isa => 'Str',
                              metaclass => 'MooseX::MetaDescription::Meta::Attribute',
                              description => {
                                  is_transient => 1
                              });
    
    __PACKAGE__->make_persistent;
    __PACKAGE__->belongs_to('artist' => 'VRPipe::Artist');
    __PACKAGE__->has_many('tracks' => 'VRPipe::Track');
}

1;
