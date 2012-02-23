#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'Bio::KBase' ) || print "Bail out!\n";
}

diag( "Testing Bio::KBase $Bio::KBase::VERSION, Perl $], $^X" );
