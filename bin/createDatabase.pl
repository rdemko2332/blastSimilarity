#!/usr/bin/perl

use strict;
use Getopt::Long;

my $dbFile;
my $blastProgram;
my $databaseType;

GetOptions("dbFile=s" => \$dbFile,
	   "blastProgram=s" => \$blastProgram	
	   );

my $blastVendor ="ncbi";

system("cat $dbFile > newdb.fasta");

if($blastProgram eq "blastp" || $blastProgram eq "blastx" || $blastProgram eq "rpstblastn" || $blastProgram eq "rpsblast" || $blastProgram eq "psiblast"){
  $databaseType = "prot";
}elsif($blastProgram eq "tblastn" || $blastProgram eq "tblastx" || $blastProgram eq "blastn") {
  $databaseType = "nucl";  
}

my $blastBinDir = "/usr/bin/ncbi-blast-2.13.0+/bin";

system("$blastBinDir/makeblastdb -in newdb.fasta -dbtype $databaseType");






