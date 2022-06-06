#!/usr/bin/perl

use strict;
use Getopt::Long;
my $string;
&GetOptions("string=s" => \$string);

my $newString = $string;
$newString =~ s/\|//g;

print "$newString $string";
`mv '$string' $newString`;
