#!/home/ivan/bin/perl

use blib;
use Chemistry::File::SMILES;

my $mol = Chemistry::Mol->parse($ARGV[0] || 'C1[13CH]C1OC', format => 'smiles');
print $mol->formula, "\n";
print $mol->print;
