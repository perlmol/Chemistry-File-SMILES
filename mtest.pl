#!/home/ivan/bin/perl

use blib;
use Chemistry::File::SMILES;

print("Please give a SMILES string\n"), exit unless $ARGV[0];

my $i=1;
my $my_parser = Chemistry::File::SMILES->new(
    add_atom => sub {shift; local $"=','; print "ATOM$i(@_)\n";$i++;},
    add_bond => sub {shift; local $"=','; print "BOND(@_)\n"}
);

print "$ARGV[0]\n";
$my_parser->parse($ARGV[0], 'mol');
