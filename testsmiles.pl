#!/home/ivan/bin/perl

use blib;
use Chemistry::File::SMILES;

my $mol = Chemistry::Mol->parse('C1CC1OC', format => 'smiles');
print $mol->formula, "\n";
print $mol->print;
