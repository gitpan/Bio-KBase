use strict;
use warnings;

use Test::More tests => 5;
use Data::Dumper;

use Bio::KBase::CentralStore;
my $obj;
my $return;
my @id_keys;

#
#  Test 1 - Can a new object be created without parameters? 
#
$obj = Bio::KBase::CentralStore->new(); # create a new object
ok( defined $obj, "Did an object get defined" );               

#
#  Test 2 - Is the object in the right class?
#
isa_ok( $obj, 'Bio::KBase::CentralStore', "Is it in the right class" );   

#
#  Test 3 - Can the object do all of the methods
#

can_ok($obj, qw[    fids_to_functions
    genomes_to_fids
]);

#
#  Test 4 - Can a new object be created with valid parameter? 
#

my $id_server = Bio::KBase::CentralStore->new();
ok( defined $id_server, "Did an object get defined" );               
#
#  Test 5 - Is the object in the right class?
#
isa_ok( $id_server, 'Bio::KBase::CentralStore', "Is it in the right class" );   
