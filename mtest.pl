#!/usr/bin/perl

use blib;
use Chemistry::Smiles;

print("Please give a SMILES string\n"), exit unless $ARGV[0];

my $i=1;
my $my_parser = new Chemistry::Smiles(
    add_atom => sub {shift; local $"=','; print "ATOM$i(@_)\n";$i++;},
    add_bond => sub {shift; local $"=','; print "BOND(@_)\n"}
);

print "$ARGV[0]\n";
$my_parser->parse($ARGV[0], 'mol');
