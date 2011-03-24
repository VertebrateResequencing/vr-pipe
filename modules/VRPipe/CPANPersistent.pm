=head1 DESCRIPTION

Just testing out Persistent

=head1 AUTHOR

Sendu Bala: sb10 at sanger ac uk

=cut

package VRPipe::Persistent;

### we are a subclass of an all-powerful Persistent class ###
#use Persistent::DBM;
#@ISA = qw(Persistent::DBM);
use Persistent::MySQL;
@ISA = qw(Persistent::MySQL);
use DBI;

sub initialize {   ### ALWAYS implement this method ###
    my $this = shift;
    
    ### call any ancestor initialization methods ###
    $this->SUPER::initialize(@_);
    
    ### define attributes of the object ###
    $this->add_attribute('firstname', 'ID',
                         'VarChar',  undef, 10);
    $this->add_attribute('lastname',  'ID',
                         'VarChar',  undef, 20);
    $this->add_attribute('telnum',    'Persistent',
                         'VarChar',  undef, 15);
    $this->add_attribute('bday',      'Persistent',
                         'DateTime', undef);
    $this->add_attribute('age',       'Transient',
                         'Number',   undef, 2);
    
    # create the table in MySQL if necessary
    my $mysql_check = q[SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_NAME = ?];
    my $mysql_create = q[create table test ( firstname VARCHAR(10), lastname VARCHAR(20), telnum VARCHAR(15), bday DateTime, PRIMARY KEY(firstname, lastname) )];
    my $dbh = DBI->connect($_[0], $_[1], $_[2]) or die $DBI::errstr;
    my @table_cols = @{$dbh->selectall_arrayref($mysql_check, undef, $_[3]) || []};
    if (! @table_cols) {
        warn "creating table\n";
        $dbh->do($mysql_create);
    }
    else {
        warn "table existed with ", scalar(@table_cols), " cols\n";
    }
}

sub print {  ### custom method ###
  my $this = shift;

  printf("%-10s %-10s %15s %s %2s\n",
         defined $this->firstname ? $this->firstname : 'undef',
         defined $this->lastname ? $this->lastname : 'undef',
         defined $this->telnum ? $this->telnum : 'undef',
         defined $this->bday ? $this->bday : 'undef',
         defined $this->age ? $this->age : 'undef');
}


1;
