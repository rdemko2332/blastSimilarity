#!/usr/bin/perl

use strict;
use warnings;
use Bio::SearchIO;

# Usage information
die "Usage: $0 <BLAST-report-file> <number-of-top-hits> <output-file>\n", if (@ARGV != 3);

my ($infile,$numHits,$outfile) = @ARGV;
print "Parsing the BLAST result ...";
my $in = Bio::SearchIO->new(-format => 'blastxml', -file => $infile);
open (OUT,">$outfile") or die "Cannot open $outfile: $!";

while ( my $result = $in->next_result ) {
	# the name of the query sequence
        print OUT "query_name: " . $result->query_description() . "\n";
        # the length of the query sequence
    	print OUT "query_length: " . $result->query_length . "\n";

        # output "no hits found" if there is no hits
    	if ( $result->num_hits == 0 ) {
		print OUT "No hits found\n";
    	} else {
		my $count = 0;

                # process each hit recursively
		while (my $hit = $result->next_hit) {
			print OUT "\t" if ($count > 0);
                        # get the accession numbers of the hits
			print OUT "accession: " . $hit->description() . "\n";
                        # get the lengths of the hit sequences
                        print OUT "hit length: " . $hit->length . "\n";
                        # get the description of the hit sequences
			print OUT "description: " . $hit->description . "\n";
                        # get the E value of the hit
			print OUT "E : " . $hit->significance . "\n";
                        #get the bit score of the hit
			print OUT "bits: " . $hit->bits() . "\n";
                        #get the score of the hit
                        print OUT "score: " . $hit->raw_score() . "\n";

                        my $hspcount = 0;

                        # process the top HSP for the top hit
			while (my $hsp = $hit->next_hsp) {
                        	print OUT "\t\t\t\t\t\t\t\n", if ($hspcount > 0);
                        	# get the frame of the query sequence
				print OUT "query string: " . $hsp->query_string() . "\n";
                                print OUT "hit string: " . $hsp->hit_string() . "\n";
                                # get the start and the end of the query sequence in the alignment
				print OUT "hsp query start: " . $hsp->start('query') . "\t hsp query end: " . $hsp->end('query'). "\n";
                                # get the start and the end of the hit sequence in the alignment
				print OUT "hit start: " . $hsp->start('hit') . "\t hit end: " . $hsp->end('hit') . "\n";
                                # get the similarity value
				print OUT "Similarity Value: ";
				printf OUT "%.1f" , ($hsp->frac_conserved * 100);
				print OUT "%\n";
                                # get the identity value
				print OUT "Identity Value: ";
				printf OUT "%.1f" , ($hsp->frac_identical * 100);
		       		print OUT "%\n";
                                print OUT "Num Identities: " . $hsp->frac_identical('hit') . "\n";
                                $hspcount++;
                        }
			$count++;

                        # flow control for the number of hits needed
			last if ($count == $numHits);
		}
    	}
}
close OUT;
print " DONE!!!\n";
