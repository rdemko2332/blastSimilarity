#!/usr/bin/perl

use strict;
use Getopt::Long;
my $string;
my $newString
&GetOptions("string=s" => \$string);

$newString = $string;
$newString =~ s/\|//g;

print "$newString $string";
`mv '$string' $newString`;
