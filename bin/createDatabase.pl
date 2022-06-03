#!/usr/bin/perl

use strict;
use Getopt::Long;

my $dbFile;
my $blastProgram;
my $blastVendor;
my $blastBinDir;
my $databaseType;

GetOptions("dbFile=s" => \$dbFile,
	   "blastProgram=s" => \$blastProgram	
	   );

$blastVendor ="ncbi";

system("cat $dbFile > newdb.fasta");

if($blastProgram eq "blastp" || $blastProgram eq "blastx" || $blastProgram eq "rpstblastn" || $blastProgram eq "rpsblast" || $blastProgram eq "psiblast"){
  $databaseType = "prot";
}elsif($blastProgram eq "tblastn" || $blastProgram eq "tblastx") {
  $databaseType = "nucl";  
}

$blastBinDir = "/usr/bin/ncbi-blast-2.13.0+/bin";

system("$blastBinDir/makeblastdb -in newdb.fasta -dbtype $databaseType");






