#!/usr/bin/perl

use strict;
use Getopt::Long;

my $dbFile;
my $databaseType;
my $blastVendor;
my $blastBinDir
    
&GetOptions("dbFile=s" => \$dbFile,
	    "databaseType=s" => \$databaseType,
            "blastVendor=s" => \$blastVendor);

system("cat $dbFile > newdb.fasta");

if($blastVendor eq "ncbi"){
    if($databaseType eq "protein"){
        $databaseType = "prot";
    } elsif($databaseType eq "nucleotide"){
        $databaseType = "nucl";
    }
    $blastBinDir = "/usr/bin/ncbi-blast-2.13.0+/bin/";
    system("$blastBinDir/makeblastdb -in newdb.fasta -dbtype $databaseType");
} elsif ($blastVendor eq "wu") {
    if($databaseType eq "protein"){
        $databaseType = "-p";
    } elsif($databaseType eq "nucleotide"){
        $databaseType = "-n";
    }
    $blastBinDir = "/usr/bin/";
    system("$blastBinDir/xdformat $databaseType newdb.fasta");
} 





