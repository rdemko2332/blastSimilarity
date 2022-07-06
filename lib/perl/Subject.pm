#!/usr/bin/perl 

## perl module for dealing with blast results...
##
## currently will be meant to  only deal with blastn as is meant to do analysis
## of consistency of mRNA vs Genomic sequence of mRNA vs mRNA to determine differential
## splicing things....

## Brian Brunk 5/17/99

package Subject;

use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(new setID getID setLength getLength setQueryLength getQueryLength setDescription getDescription getLinesOfLength getUnbrokenSubstring addHSP getHSPs getPValue getTotalHSPLength  getNonoverlappingQueryMatchLength getNonoverlappingSubjectMatchLength getTotalHSPPercent getTotalIdentities getTotalMatchLength getTotalPositives  getEndsMissed  getMinQueryStart getMaxQueryEnd getMinSubjectStart getMaxSubjectEnd  getQuery5endMis getQuery3endMis getQueryEndsMissed getSubjectEndsMissed getSubject5endMis getSubject3endMis  getBestHSPLength getBestHSP getBestHSPPercent  getTotalQueryGaps getTotalSubjectGaps getLongestHsp getGoodMatch getConsistentEnds checkRepeatMaskedLength getSubjectStats getNeighborStats getSimilaritySummary  getSimilaritySpans getSummaryStats compareQueryVsSubject getHSPsFromIDArray getHSPFromID setHSPs countHSPs getSubjectDirection removeWrongOrientHSPs countGoodHSPs getGoodHSPsBySubjectStart  getGoodHSPsByQueryStart  getHSPsBySubjectStart getHSPsByQueryStart checkOrderConsistency checkOrientationConsistency checkConsistency  countContainedQueryHSPs countContainedSubjectHSPs checkQueryGaps  checkSubjectGaps getConsistentGeneModel checkGoodCGM getFirstCorrectOrientedHSP getLastCorrectOrientedHSP compareHSPConsistency getBestCGM getPercentHSPsCoveredByCGM getCGMcoverage getMaxForwardCoverage getPercentCGMCoverage  getConsistentParents toString);

use HSP;
use strict;

my $debug = 0;

##datastructure....
##{id} = identifier of subject
##{length} = length of subject
##{hsp}->[arrayIndex] ..will increment by 1 each new HSP
##     [0]=[queryStart,queryEnd]
##     [1]=[subjectStart,subjectEnd]
##     [2]=[matchLength,matchPercent,pValue,direction] 

##simple constructor...need to make more robust if inheriting for other purposes

sub new {
  my($class,$id) = @_;
  my $self = {};
  bless $self, $class;
  $self->{"id"} = $id;
  $self->{"numHSPs"} = 0;
  return $self;
}

sub setID{
  my($self,$id) = @_;
  $self->{"id"} = $id;
}

sub getID{
  my $self = shift;
  return $self->{"id"};
}

sub setLength{
  my($self,$length) = @_;
  $self->{"length"} = $length;
}

sub getLength{
  my $self = shift;
  return $self->{"length"};
}

sub setQueryLength{  ##a little silly but need to calculate the query 3'end
  my($self,$length) = @_;
  $self->{"queryLength"} = $length;
}

sub getQueryLength{
  my $self = shift;
  return $self->{"queryLength"};
}   

sub setDescription {
  my($self,$desc) = @_;
  $desc =~ s/\s+/ /g;
  $self->{'description'} = $desc;
}

sub getDescription {
  my($self,$ll,$indent) = @_;
  return $ll ? $self->getLinesOfLength($self->{'description'},$ll,$indent) : $self->{description};
}

my $spaces = "                                                  ";
sub getLinesOfLength {
  my($self,$text,$ll,$indent) = @_;
  my $len = $ll - $indent;
  my $sp = substr($spaces,0,$indent);
  my $ret;
  my $line;
  while($text){
    ($line,$text) = $self->getUnbrokenSubstring($text,$len);
    $ret .= "$sp$line\n";
  }
  chomp $ret;
  return $ret;
}

sub getUnbrokenSubstring {
  my ($self,$txt,$ll) = @_;
  return ($txt,"") if length($txt) <= $ll;
  my @t = split("",$txt);
  for($ll;$ll > 0;$ll--){
    last if $t[$ll] =~ /(\s|-)/;
  }
  $ll++ if $t[$ll] eq "-";
  return (substr($txt,0,$ll),substr($txt,$t[$ll] eq " " ? $ll+1 : $ll))
}

##each HSP gets an identifier "ID" which is the number it is in the set
##first one is 1;
##keep HSPs as a hash with the HSP id as key...
sub addHSP{
  my $self = shift;
  $self->{"numHSPs"}++;		##count the number of hsps
	my $hsp =  HSP->new($self->getLength(),$self->{"numHSPs"},@_);
  $self->{"hspHash"}->{$self->{"numHSPs"}} = $hsp;
  push(@{$self->{"hsp"}},$hsp);

}

sub getHSPs{
  my $self = shift;
  return @{$self->{"hsp"}};
}

sub getPValue{
  my $self = shift;
  return $self->{"hsp"}->[0]->getPValue();
}

## need methods that get the total matchlength, total match percent (normalized
## for length)

##don't want to sum the overlaps..get max - min - exon length
sub getTotalHSPLength{
  my $self = shift;
  if(!exists $self->{"totalHSPLength"}){
    my @hsps = $self->getHSPsBySubjectStart();
    my $first = shift @hsps;
    return 0 unless $first;
    my $start = $first->getSubjectStart(); 
    my $end = $first->getSubjectEnd(); 
    my $len = 0;
    foreach my $h (@hsps){
      next if $h->getSubjectEnd() <= $end; ##does not extend
      if($h->getSubjectStart <= $end){  ##overlaps
        $end = $h->getSubjectEnd();  #extend end ... already dealt with if new end is less
      }else{  ##there is a gap in between ..
        $len += $end - $start + 1;
        $end = $h->getSubjectEnd();
        $start = $h->getSubjectStart();
      }
    }
    $len += $end - $start + 1; # deal with the last one 
    $self->{"totalHSPLength"} = $len;
  }
  return $self->{"totalHSPLength"};
}

sub getNonoverlappingQueryMatchLength {
  my $self = shift;
  if(!exists $self->{"totalQueryMatchLength"}){
    my @hsps = $self->getHSPsByQueryStart();
    my $first = shift @hsps;
    return 0 unless $first;
    my $start = $first->getQueryStart(); 
    my $end = $first->getQueryEnd(); 
    my $len = 0;
    foreach my $h (@hsps){
      next if $h->getQueryEnd() <= $end; ##does not extend
      if($h->getQueryStart <= $end){  ##overlaps
        $end = $h->getQueryEnd();  #extend end ... already dealt with if new end is less
      }else{  ##there is a gap in between ..
        $len += $end - $start + 1;
        $end = $h->getQueryEnd();
        $start = $h->getQueryStart();
      }
    }
    $len += $end - $start + 1; # deal with the last one 
    $self->{"totalQueryMatchLength"} = $len;
  }
  return $self->{"totalQueryMatchLength"};
}

sub getNonoverlappingSubjectMatchLength {
  my $self = shift;
  return $self->getTotalHSPLength();
}

sub getTotalHSPPercent{
  my $self = shift;
  if(!exists $self->{"totalHSPPercent"}){
    my $tmpP = 0;
    my $tmpL = 0;
    foreach my $h ($self->getHSPs()){
      $tmpP += ($h->getMatchPercent() * $h->getMatchLength());
      $tmpL += $h->getMatchLength();
    }
    $self->{"totalHSPPercent"} = $tmpP / $tmpL;
  }
  return $self->{"totalHSPPercent"};
}

sub getTotalIdentities {
	my $self = shift;
	if(!exists $self->{'totalIdentities'}){
		foreach my $h ($self->getHSPs()){
			$self->{'totalIdentities'} += $h->getIdentities();
		}
	}
	return $self->{'totalIdentities'};
}

sub getTotalMatchLength {
	my $self = shift;
	if(!exists $self->{'totalMatchLength'}){
		foreach my $h ($self->getHSPs()){
			$self->{'totalMatchLength'} += $h->getMatchLength();
		}
	}
	return $self->{'totalMatchLength'};
}

sub getTotalPositives {
	my $self = shift;
	if(!exists $self->{'totalPositives'}){
		foreach my $h ($self->getHSPs()){
			$self->{'totalPositives'} += $h->getPositives();
		}
	}
	return $self->{'totalPositives'};
}


##want methods that retrieve the amount of non-matching bases at each end
##of both the subject and query...
##query can be reverse complemented so what are the semantics
##  Make it relative to the input query sequence so will be same
##  for all subjects!!
## this already is how the query is represented so need to
## reverse the subject if direction is 0

##NOTE THAT AM  NOW REPRESENTING THE ORIENTATION AND NUMBERING OF BOTH
##THE QUERY AND SUBJECT CORRECTLY SO DON'T NEED TO CORRECT NOW!!

sub getEndsMissed{
  my $self = shift;
  my($q3e,$q5e,$s3e,$s5e)  = (0,1e10,0,1e10);
  foreach my $h ($self->getHSPs()){
    $q5e = $h->getQueryStart() if $h->getQueryStart() < $q5e;
    $q3e = $h->getQueryEnd() if $h->getQueryEnd() > $q3e;
    $s5e = $h->getSubjectStart() if $h->getSubjectStart() < $s5e;
    $s3e = $h->getSubjectEnd() if $h->getSubjectEnd() > $s3e;
  }
  ##query
  $self->{"minQstart"} = $q5e; ##record these so can get to later...
  $self->{"maxQend"} = $q3e;
  $self->{"Q5endMis"} = $q5e - 1;
  $self->{"Q3endMis"} = $self->getQueryLength() - $q3e;
  ##subject
  $self->{"minSstart"} = $s5e; ##record these so can get to later...
  $self->{"maxSend"} = $s3e;
  $self->{"S5endMis"} = $s5e - 1;
  $self->{"S3endMis"} = $self->getLength() - $s3e;
}

sub getMinQueryStart{
  my $self = shift;
  if(!exists $self->{"minQstart"}){
    $self->getEndsMissed();
  }
  return $self->{"minQstart"};
}

sub getMaxQueryEnd{
  my $self = shift;
  if(!exists $self->{"maxQend"}){
    $self->getEndsMissed();
  }
  return $self->{"maxQend"};
}

sub getMinSubjectStart{
  my $self = shift;
  if(!exists $self->{"minSstart"}){
    $self->getEndsMissed();
  }
  return $self->{"minSstart"};
}

sub getMaxSubjectEnd{
  my $self = shift;
  if(!exists $self->{"maxSend"}){
    $self->getEndsMissed();
  }
  return $self->{"maxSend"};
}

sub getQuery5endMis{
  my $self = shift;
  if(!exists $self->{"Q5endMis"}){
    $self->getEndsMissed();
  }
  return $self->{"Q5endMis"}
}

sub getQuery3endMis{
  my $self = shift;
  if(!exists $self->{"Q3endMis"}){
    $self->getEndsMissed();
  }
  return $self->{"Q3endMis"}
}

sub getQueryEndsMissed{
  my($self,$slop) = @_;
  my $cnt = 0;
  if(!exists $self->{'queryEndsMissed'}){
    $cnt++ if $self->getQuery5endMis() > $slop;
    $cnt++ if $self->getQuery3endMis() > $slop;
    $self->{'queryEndsMissed'} = $cnt;
  }
  return $self->{'queryEndsMissed'};
}

sub getSubjectEndsMissed{
  my($self,$slop) = @_;
  my $cnt = 0;
  if(!exists $self->{'subjectEndsMissed'}){
    $cnt++ if $self->getSubject5endMis() > $slop;
    $cnt++ if $self->getSubject3endMis() > $slop;
    $self->{'subjectEndsMissed'} = $cnt;
  }
  return $self->{'subjectEndsMissed'};
}

sub getSubject5endMis{
  my $self = shift;
  if(!exists $self->{"S5endMis"}){
    $self->getEndsMissed();
  }
  return $self->{"S5endMis"}
}

sub getSubject3endMis{
  my $self = shift;
  if(!exists $self->{"S3endMis"}){
    $self->getEndsMissed();
  }
  return $self->{"S3endMis"}
}

sub getBestHSPLength {
  my $self = shift;
  return $self->{"hsp"}->[0]->getMatchLength();
}

sub getBestHSP {
	my $self = shift;
	return $self->{"hsp"}->[0];
}

sub getBestHSPPercent{
  my $self = shift;
  return $self->{"hsp"}->[0]->getMatchPercent();
}

sub getTotalQueryGaps{
  my $self = shift;
  if(!exists $self->{"totalQGaps"}){
    $self->{"totalQGaps"} = 0;
    my @sort = $self->getHSPsByQueryStart();
    for(my $i=0;$i<$self->countHSPs()-1;$i++){
      $self->{"totalQGaps"} += ($sort[$i+1]->getQueryStart() - $sort[$i]->getQueryEnd() < 0 ? 0 : $sort[$i+1]->getQueryStart() - $sort[$i]->getQueryEnd() - 1);
    }
  }
  return $self->{"totalQGaps"};
}

sub getTotalSubjectGaps{
  my $self = shift;
  my $tot = 0;
  if(!exists $self->{"totalSGaps"}){
    $self->{"totalSGaps"} = 0;
    my @sort = $self->getHSPsBySubjectStart();
    for(my $i=0;$i<$self->countHSPs()-1;$i++){
      $self->{"totalSGaps"} += ($sort[$i+1]->getSubjectStart() - $sort[$i]->getSubjectEnd() < 0 ? 0 : $sort[$i+1]->getSubjectStart() - $sort[$i]->getSubjectEnd() - 1);
    }
  }
  return $self->{"totalSGaps"};

}

sub getLongestHsp {
  my($self) = @_;
  if(!exists $self->{longestHsp}){
    my @sort = sort{$b->getMatchLength()  <=> $a->getMatchLength()} $self->getHSPs();
    $self->{longestHsp} = $sort[0];
  }
  return $self->{longestHsp};

}

##now want to determine if meets requirements for clustering...
##totalHSPLength >= 40, ends??, percent ID taken care of by upfront limitatins

## need to block out repeat regions....
  ##need to generalize this to do more than one HSP!!!
  ##also need to take in more than one repeat!!

sub getGoodMatch{
  my($self,$lengthCut,$percentCut,$endCut,$rep) = @_;  ##parameters...
  if(!exists $self->{"clusterMatch"}){
    $lengthCut = 40 if !defined $lengthCut;
    $percentCut = defined $percentCut ? $percentCut : 92; 
    $endCut = defined $endCut ? $endCut : 15; 
    $self->{"clusterMatch"} = 0;
    if($self->getTotalHSPPercent() >= $percentCut && $self->checkRepeatMaskedLength($lengthCut,@$rep) == 1){
      ##following  checks the end matches if desired..
      $self->{"clusterMatch"} = 1;
      if(($self->getQuery5endMis() <= $endCut && $self->getQuery3endMis() <= $endCut) || ($self->getSubject5endMis() <= $endCut && $self->getSubject3endMis() <= $endCut) || ($self->getQuery5endMis() <= $endCut && $self->getSubject3endMis() <= $endCut) || ($self->getQuery3endMis() <= $endCut && $self->getSubject5endMis() <= $endCut)){
        $self->{consistentEnds} = 1;
      }else{ 
        $self->{consistentEnds} = 0;
      }
    }else{ $self->{clusterMatch} = 0;}
  }
  return $self->{"clusterMatch"};
}

sub getConsistentEnds {
  my $self = shift;
  return $self->{consistentEnds} if exists $self->{consistentEnds};
#  print STDERR "you must call Subject->getGoodMatch() before Subject->getConsistentEnds()\n";
  return undef;
}

##only need to go until meets length then return....
##retruns 1 if OK...
##am assuming that HSPs don't overlap much if at all given blast parameters...
sub checkRepeatMaskedLength{
  my($self,$length, @rep) = @_;
  my $totLength = 0;
  foreach my $h ($self->getHSPsByQueryStart()){
    my $hLen = $h->getMatchLength(); ##start with length and subtract if repeat..
    print STDERR "HSPLength: $hLen\n" if $debug == 1;
    foreach my $r (@rep){
      if($r->getRepeatStart() <= $h->getQueryStart() && $r->getRepeatEnd() >= $h->getQueryStart()){  ##starts before and extends into or beyond the hsp
        print STDERR " HSP:\(",$h->getQueryStart(),":",$h->getQueryEnd(),"\), Repeat:\(",$r->getRepeatStart(),":",$r->getRepeatEnd(),"\)\n" if $debug == 1;
        if($h->getQueryEnd() > $r->getRepeatEnd()){ ##ends within the hsp
          $hLen = $hLen - ($r->getRepeatEnd() - $h->getQueryStart());
        }else{ ##repeat covers hsp....$hLen = 0;
          $hLen = 0;
        }
      }elsif($r->getRepeatStart() >= $h->getQueryStart() && $r->getRepeatStart() <= $h->getQueryEnd()){  ##starts within the hsp and ends within or beyond hsp
        print STDERR " HSP:\(",$h->getQueryStart(),":",$h->getQueryEnd(),"\), Repeat:\(",$r->getRepeatStart(),":",$r->getRepeatEnd(),"\)\n" if $debug == 1;
        if($r->getRepeatEnd() > $h->getQueryEnd()){  ##extends beyond hsp
          $hLen = $hLen - ($h->getQueryEnd() - $r->getRepeatStart());
        }else{ ##ends within hsp..subtract the entire repeat
          $hLen = $hLen - $r->getRepeatLength();
        }
      }
    }
    print STDERR " Result: $hLen\n" if $debug == 1;
    $totLength += $hLen;
    return 1 if $totLength >= $length; ##check after each hsp to see if meets length min
  }
  return 0;
}

sub getSubjectStats{ 
  my $self = shift;
  if(!exists $self->{"subjectStats"}){
    @{$self->{"subjectStats"}} = ($self->getID(),$self->getLength(),$self->getSubjectDirection(),$self->countHSPs(),$self->getPValue(),$self->getTotalHSPLength(),int($self->getTotalHSPPercent()),$self->getTotalQueryGaps(),$self->getQuery5endMis(),$self->getQuery3endMis(),$self->getTotalSubjectGaps(),$self->getSubject5endMis(),$self->getSubject3endMis());
  }
  return @{$self->{"subjectStats"}};
}

##summary
sub getNeighborStats{
  my $self = shift;
  ##changing this so that it returns the longest HSP rather than the sum of them all...
  return $self->getID().":".$self->getPValue().":".$self->getLongestHsp()->getMatchLength().":".int($self->getLongestHsp()->getMatchPercent()).":".$self->getConsistentEnds();
#  return $self->getID().":".$self->getPValue().":".$self->getTotalHSPLength().":".int($self->getTotalHSPPercent()).":".$self->getConsistentEnds();
}

##note that reading frame will be null for blastn type queries.
##$d is the delimiter...
sub getSimilaritySummary {
	my($self,$d) = @_;
	if(!$d){ $d = ", "; }
	my $tmp = "  Sum: ".$self->getID().$d.$self->getBestHSP()->getScore().$d.$self->getBestHSP()->getPValue().$d;
	$tmp .= $self->getMinSubjectStart().$d.$self->getMaxSubjectEnd().$d.$self->getMinQueryStart().$d.$self->getMaxQueryEnd().$d;
	$tmp .= $self->countHSPs().$d.$self->getTotalMatchLength().$d.$self->getTotalIdentities().$d.$self->getTotalPositives().$d;
	$tmp .= ($self->getBestHSP()->getDirection() ? 0 : 1).$d.$self->getBestHSP()->getFrame().$d.$self->getNonoverlappingQueryMatchLength().$d.$self->getNonoverlappingSubjectMatchLength();
        # add percentage match to shortest sequence
        $tmp .= $d.(int(($self->getLength() < $self->getQueryLength() ? $self->getNonoverlappingSubjectMatchLength() / $self->getLength() : $self->getNonoverlappingQueryMatchLength() / $self->getQueryLength()) * 10000) / 100 );
	return $tmp;
}

sub getSimilaritySpans {
	my($self,$d) = @_;
	my @tmp;
	my $count = 0;
	foreach my $h ($self->getHSPs()){
		$count++;
		push(@tmp,"   HSP$count: ".$self->getID().$d.$h->getSimilaritySpan($d));
	}
	return join("\n",@tmp);
}

sub getSummaryStats{
  my $self = shift;
  return $self->getID().": PV=".$self->getPValue().", ML=".$self->getTotalHSPLength().", \%ID=".int($self->getTotalHSPPercent()).", MaxSpan=\(query\(".$self->getMinQueryStart().",".$self->getMaxQueryEnd()."\), sub\(".$self->getMinSubjectStart().",".$self->getMaxSubjectEnd()."\)\)";
}

##compare query to subject to see if likely repeat sequencing of same clone..
sub compareQueryVsSubject{
  my $self = shift;
  my $db = 0;
  if($db == 1){
    print STDERR "compareQueryVsSubject:\n";
    print STDERR "  totalHSPPercent: ",$self->getTotalHSPPercent(),"\n";
    print STDERR "  lengthComparison: ",abs($self->getQueryLength() - $self->getLength()) / $self->getQueryLength(),"\n";
    print "  matchComparison: ",abs($self->getQueryLength() - $self->getTotalHSPLength()) / $self->getQueryLength(),"\n";
  }        
  if(!exists $self->{"QvsS"}){
    $self->{"QvsS"} = 1; ##start off with one then set to 0 if fails one test
    if($self->getTotalHSPPercent() < 95){
      $self->{"QvsS"} = 0;
    }elsif(abs($self->getQueryLength() - $self->getLength()) / $self->getQueryLength() > 0.05){
      $self->{"QvsS"} = 0;
    }elsif(abs($self->getQueryLength() - $self->getTotalHSPLength()) / $self->getQueryLength() > 0.05){
      $self->{"QvsS"} = 0;
    }
  }
  print STDERR "  RETURNS: $self->{'QvsS'}\n" if $db == 1;
  return $self->{"QvsS"};
}

sub getHSPsFromIDArray{
  my($self,@ids) = @_;
  my @ret;
  my $id;
  foreach my $id (@ids){
    push(@ret,$self->{"hspHash"}->{$id});
  }
  return @ret;
}

sub getHSPFromID{
  my($self,$id) = @_;
  return $self->{"hspHash"}->{$id};
}

sub setHSPs{
  my $self = shift;
  @{$self->{"hsp"}} = @_;
  $self->{"numHSPs"} = scalar(@{$self->{"hsp"}});
}


sub countHSPs{
  my($self) = @_;
  return $self->{"numHSPs"};
}

##gets the direction ofthe best HSP!!
sub getSubjectDirection{
  my $self = shift;
  if(!exists $self->{"subjectDirection"}){
    my @hsps = $self->getHSPs();
    $self->{"subjectDirection"} = $hsps[0]->getDirection();
  }
  return $self->{"subjectDirection"};
}

##note that want to always get the HSPs relative to the orientation of the subject
##even if getting by query locations.
##the query get reverse complemented when it is the minus (0) strand
##so need to take this into account when gettting the query locations..
sub removeWrongOrientHSPs{
  my $self = shift;
  my(@allhsps) = @_;
  my @hsps;
  my $h;
  foreach my $h (@allhsps){
    push(@hsps, $h) if $h->getDirection() == $self->getSubjectDirection();
  }
  return @hsps;
}

sub countGoodHSPs{
  my $self = shift;
  if(!exists $self->{"goodHSPsBySbjct"}){
    @{$self->{"goodHSPsBySbjct"}} = $self->removeWrongOrientHSPs($self->getHSPsBySubjectStart());
  }
  return scalar(@{$self->{"goodHSPsBySbjct"}});
}

sub getGoodHSPsBySubjectStart{
  my $self = shift;
  if(!exists $self->{"goodHSPsBySbjct"}){
    @{$self->{"goodHSPsBySbjct"}} = $self->removeWrongOrientHSPs($self->getHSPsBySubjectStart());
  }
  return @{$self->{"goodHSPsBySbjct"}};
}

sub getGoodHSPsByQueryStart{
  my $self = shift;
  if(!exists $self->{"goodHSPsByQuery"}){
    @{$self->{"goodHSPsByQuery"}} = $self->removeWrongOrientHSPs($self->getHSPsByQueryStart());
  }
  return @{$self->{"goodHSPsByQuery"}};
}

sub getHSPsBySubjectStart{
  my $self = shift;
  if(!exists $self->{"sortedHSPsBySbjct"}){
    @{$self->{"sortedHSPsBySbjct"}} = sort{$a->getSubjectStart() <=> $b->getSubjectStart()}$self->getHSPs();
    $self->{"countGoodOrientHSPs"} = scalar(@{$self->{"sortedHSPsBySbjct"}});
  }
  return @{$self->{"sortedHSPsBySbjct"}};
}


sub getHSPsByQueryStart{
  my $self = shift;
  if(!exists $self->{"sortedHSPsByQuery"}){
    @{$self->{"sortedHSPsByQuery"}} = sort{$a->getQueryStart() <=> $b->getQueryStart()}$self->getHSPs();
  }
  return @{$self->{"sortedHSPsByQuery"}};
}

sub checkOrderConsistency{
  my($self,@hsps) = @_;
  my @ss;
  my @qs;
  if(scalar(@hsps) < 1){
    @ss = $self->getGoodHSPsBySubjectStart();
    @qs = $self->getGoodHSPsByQueryStart();
  }else{
    @ss = sort{$a->getSubjectStart() <=> $b->getSubjectStart()}@hsps;
    @qs = sort{$a->getQueryStart() <=> $b->getQueryStart()}@hsps;
  }
##  print STDERR "Num HSPs: ",$self->countHSPs(),", Subjects: ",scalar(@ss),", Query: ",scalar(@qs),"\n";
  my $i;
  for(my $i = 0;$i<scalar(@qs);$i++){
    if($qs[$i]->getID() != $ss[$i]->getID()){return 0;}
  }
  return 1;
}

##check that the orientation does not change for either query or subject
sub checkOrientationConsistency{
  my $self = shift;
  my @qs = $self->getHSPs();
  my $i;
  for(my $i=0;$i<$self->countHSPs()-1;$i++){
    if($qs[$i]->getDirection() != $qs[$i+1]->getDirection()){
      return 0;			##not consistent as direction has flipped for either query or subject
    }
  }
  return 1;
}

##checks both for order and orientation
sub checkConsistency{
  my $self= shift;
  if($self->checkOrderConsistency() && $self->checkOrientationConsistency()){
    return 1;
  }
  return 0;			##one or the other of checks is 0;
}

sub countContainedQueryHSPs{
  my $self = shift;
  my $count = 0;
  my @qs = $self->getHSPsByQueryStart();
  my $i;
  for(my $i=0;$i<$self->countHSPs()-1;$i++){
    if($qs[$i+1]->getQueryEnd() <= $qs[$i]->getQueryEnd()){
      #      print "DEBUG contained: ",$qs[$i+1]->getQueryEnd()," <= ",$qs[$i]->getQueryEnd(),"\n";
      $count++;
    }
  }
  return $count;
}

sub countContainedSubjectHSPs{
  my $self = shift;
  my $count = 0;
  my @qs = $self->getHSPsBySubjectStart();
  for(my $i=0;$i<$self->countHSPs()-1;$i++){
    if($qs[$i+1]->getSubjectEnd() <= $qs[$i]->getSubjectEnd()){
      #      print "DEBUG contained: ",$qs[$i+1]->getSubjectEnd()," <= ",$qs[$i]->getSubjectEnd(),"\n";
      $count++;
    }
  }
  return $count;
}

sub checkQueryGaps{
  my $self = shift;
  
  my @qs = $self->getHSPsByQueryStart();
  my $queryGap;
  my $countOverlaps = 0;
  my $countGaps = 0;
  my $subGap;
  for(my $i=0;$i<$self->countHSPs()-1;$i++){
    $queryGap = $qs[$i+1]->getQueryStart() - $qs[$i]->getQueryEnd();
    $subGap = $qs[$i+1]->getSubjectStart() - $qs[$i]->getSubjectEnd();
    if($queryGap < -10){	##overlapps by > 10 bp so count it!!
      ##should check here to see what the overlap is of the subject....if same
      ##then is likely just a region of where two different HSPs can be started and
      ##extended and blast has taken them apart unnecessarily.
      if(abs($queryGap - $subGap) > 5){
	$countOverlaps++;
      }
    }elsif($queryGap > 10){	##gap of greater than 10 bp..
      ##check to see if similar size gap in subject...if so don't count
      if(abs($queryGap - $subGap) > 5){  
	$countGaps++;
      }
    }
  }
  return ($countOverlaps,$countGaps);
}

##in an analysis of genomic sequence vs dots, I want to check the subject gaps
##as it is the subject that is the mRNA...
##returns (overlaps,gaps) (1,1) if consistent 0 if not
sub checkSubjectGaps{
  my $self = shift;
  
  my @qs = $self->getHSPsBySubjectStart();
  my $queryGap;
  my $countOverlaps = 0;	##
  my $countGaps = 0;
  my $subGap;
  for(my $i=0;$i<$self->countHSPs()-1;$i++){
    $queryGap = $qs[$i+1]->getQueryStart() - $qs[$i]->getQueryEnd();
    ##    $queryGap = $qs[$i+1]->getQueryStart() - $qs[$i]->getQueryEnd();
    $subGap = $qs[$i+1]->getSubjectStart() - $qs[$i]->getSubjectEnd();
    if($subGap < -10){		##overlapps by > 10 bp so count it!!
      ##should check here to see what the overlap is of the subject....if same
      ##then is likely just a region of where two different HSPs can be started and
      ##extended and blast has taken them apart unnecessarily.
      if(abs(abs$queryGap - abs$subGap) > 10){
	$countOverlaps++;
      }
    }elsif($subGap > 10){	##gap of greater than 10 bp..
      ##check to see if similar size gap in subject
      
      $countGaps++;
    }
  }
  return($countOverlaps,$countGaps);
}

## checking for consistent gene model as a biologist would
## ignore temporarily those HSPs that are not consistent....could analyze later
## Check that order is same and then that end of current matched beginning of next
## if not then go onto next one and see if end of that matches beginning.. etc...
## if there is more than one possibility for next step, check to see which of them
## have order and orientatino consistency with the query sequence HSPs
## this last would happen in the case of contained (largely) HSPs...
## if there is more than one possibility for next step, check to see if each can
## be extended to the next etc...
## model should be consistent to the extents (beginning and end of first and last HSP).

##the following method is limited and has been replaced by getBestCGM
sub getConsistentGeneModel{	##check"CGM"
  my($self,$overlap) = @_;
  if($overlap eq ""){$overlap = 10};
  if(!exists $self->{"CGM"}){
    my $haveFirst = 0;
    my $subDir = $self->getSubjectDirection();
    my @sa = $self->getHSPsBySubjectStart(); ##array of subjects sorted by subject start
    for(my $i=0;$i<$self->countHSPs()-1;$i++){
      next if $sa[$i]->getDirection() != $subDir; ##must be in dominant direction
      if($haveFirst == 0){
	push(@{$self->{"CGM"}},$sa[$i]);
	$haveFirst = 1;
      }
      my $add = 1;
      while(($i + $add) < $self->countHSPs()){
	if($sa[$i+$add]->getDirection() != $subDir){ ##must be in dominant direction
	  $add++;
	  next;
	}
	if(abs($sa[$i+$add]->getSubjectStart() - $sa[$i]->getSubjectEnd()) <= $overlap){ ##$overlap bp 
	  ##this is the next one that is potentially good!!
	  ##check the query gap....
	  if(($sa[$i+$add]->getQueryStart() - $sa[$i]->getQueryEnd()) >= -$overlap){
	    ##consistent with gene model
	    ##??should check the next one to see if also consistent with this one...
	    push (@{$self->{"CGM"}},$sa[$i+$add]); ##push onto good the current one am working on...
	    last;
	  }
	}
	$add++;			##increment to next one!!
      }
      $i += ($add - 1);		##it will get incremented on next round to equal $i+$add
    }
    if($haveFirst == 0){  ##is singleton or the last HSP is the only one in subject direction...
      push(@{$self->{"CGM"}},$sa[-1]);
    }
  } 
  ##should I check that covers as much of subject as possible??
  ##end of last and beginning of first correctly oriented HSP should be max span
  ##separate method "checkGoodCGM()"
  
  return @{$self->{"CGM"}};	##returns the good gene model
}

sub checkGoodCGM{		##consistent gene model..returns 1 if good, else 0;
  my $self = shift;
  if(!exists $self->{"goodCGM"}){
    my @cgm = $self->getConsistentGeneModel();
#    print STDERR "CGM size: ",scalar(@cgm),"\n";
    ##know that the first one is correct....check to see if last HSP in cgm is the
    ##last correct HSP and must be more than one HSP in model
    my$lastC = $self->getLastCorrectOrientedHSP();
    my $lastCGM = $cgm[-1];
    if($lastC->getSubjectStart() == $lastCGM->getSubjectStart() && scalar(@cgm) > 1){
      $self->{"goodCGM"} = 1;
    }else{
      $self->{"goodCGM"} = 0;
    }
  }
  return $self->{"goodCGM"};
}

##get first hsp (on subject) in correct orientation.
sub getFirstCorrectOrientedHSP{
  my $self = shift;
  ##cache...
  if(!exists $self->{"firstHSP"}){
    my @sa = $self->getHSPsBySubjectStart(); ##array of subjects sorted by subject start
    for(my $i=0;$i<$self->{"countHSPs"};$i++){
      if($sa[$i]->getDirection() == $self->getSubjectDirection()){
	$self->{"firstHSP"} = $sa[$i];
	last;
      }
    }
  }
  return $self->{"firstHSP"};
}

sub getLastCorrectOrientedHSP{
  my $self = shift;
  print STDERR "getting last hsp\n" if $debug == 1;
  if(!exists $self->{"lastHSP"}){
    my @sa = $self->getHSPsBySubjectStart();
    for(my $i=scalar(@sa)-1;$i>=0;$i--){
      print STDERR "getLastHSP: ",$sa[$i]->toString("  ") if $debug == 1;
      if($sa[$i]->getDirection() == $self->getSubjectDirection()){
	$self->{"lastHSP"} = $sa[$i];
	last;
      }
    }
  }
  return $self->{"lastHSP"};
}

##take in two hsps and compare if consistent.
##check gaps,overlaps,orientation, ?order??
sub compareHSPConsistency{
  my($self,$h1,$h2,$overlap) = @_;
  return 0 if h1->getDirection() != h2->getDirection();
  if($overlap eq ""){$overlap = 10};
  if(abs($h2->getSubjectStart() - $h1->getSubjectEnd()) <= $overlap){ ##$overlap bp 
	  ##this is the next one that is potentially good!!
	  ##check the query gap....
    if(($self->getSubjectdirection() == 1 && ($h2->getQueryStart() - $h1->getQueryEnd()) >= -$overlap) || ($h1->getQueryEnd() - $h2->getQueryStart()) >= -$overlap){
	    ##consistent with gene model
      return 1;
    }
  }
  return 0;
}

##make a matrix of all combinations...
##how to keep trackof the data???sounds like a dynamic programming task!
sub getBestCGM{
  my($self,$overlap) = @_;
  if(! exists $self->{"bestCGM"}){
     if($overlap eq ""){$overlap = 10;}
     my @gs = $self->getGoodHSPsBySubjectStart(); ##gets only the good ones...correct orient
     my%data;
     $data{$gs[0]->getID()}->[0] = $gs[0]->getID();
     for(my $i=1;$i<scalar(@gs);$i++){
       my@par = $self->getConsistentParents($gs[$i],$overlap);
       if(scalar(@par) > 1){  ##need to get one with longest cgm
	 #      print "MULTIPLE PARENTS:\n",$self->hspsToString(@par);
	 my@st = sort{ scalar(@{$data{$b->getID()}}) <=> scalar(@{$data{$a->getID()}}) }@par;
	 #      print "  Sorted:\n",$self->hspsToString(@st);
	 my$best = $st[0];
	 if(scalar(@{$data{$st[0]->getID()}}) == scalar(@{$data{$st[1]->getID()}})){  ##have at least twothe same ...
	   for(my $a=1;$a<scalar(@st);$a++){
	  ##get the one with lowest start
	     last if scalar(@{$data{$st[$a]->getID()}}) != scalar(@{$data{$st[$a-1]->getID()}});
	     $best = $st[$a] if $st[$a]->getSubjectStart() < $best->getSubjectStart();
	   }
	 }
	 $data{$gs[$i]->getID()} = [@{$data{$best->getID()}},$gs[$i]->getID()];
       }elsif(scalar(@par) == 1){
	 $data{$gs[$i]->getID()} = [@{$data{$par[0]->getID()}},$gs[$i]->getID()];
       }else{  ##there are none so..
	 $data{$gs[$i]->getID()}->[0] = $gs[$i]->getID();
       }
     }
     ##now have the data....return the one with the longest cgm (most HSPs) or if
     ##multiple ones have the same # HSPs, then return the one that covers the
     ##most of the subject
     #  print "Keys from \%data: ",join(', ',keys(%data)),"\n";
     #  foreach my $key (keys%data){
     #    print "  $key: \(",join(', ',@{$data{$key}}),")\n";
     #  }
     #  my@good = sort{scalar(@{$data{$b}}) <=> scalar(@{$data{$a}})}keys%data;
     ##@good is the keys you idiot...not the values!!!!
     my@good = sort{scalar(@{$b}) <=> scalar(@{$a})}values%data;
     if(scalar(@good) > 1){
#       print "Best cgm: ",join(', ',@{$good[0]}),"\n";
       if(scalar(@{$good[0]}) == scalar(@{$good[1]})){
	 #      print "More than one cgm of length ",scalar(@{$good[0]}),"\n";
	 my @tmp;
	 push(@tmp,$good[0]);
	 for(my $b=1;$b<scalar(@good);$b++){
	   last if scalar(@{$good[$b]}) != scalar(@{$good[$b-1]});
	   push(@tmp,$good[$b]);
	 }
	 ##need to sort @tmp into @good by CGMcoverage
	 @good = sort{$self->getCGMcoverage($self->getHSPsFromIDArray(@{$b})) <=> $self->getCGMcoverage($self->getHSPsFromIDArray(@{$a})) }@tmp;
       }
       foreach my $gd (@good){
#	       print "\@good contains\n",$self->hspsToString($self->getHSPsFromIDArray(@{$gd}));
	 @{$self->{"bestCGM"}} = $self->getHSPsFromIDArray(@{$gd});
	 last;
       }
     }elsif(scalar(@good) == 1){
       @{$self->{"bestCGM"}} = $self->getHSPsFromIDArray(@{$good[0]});
     }else{
       @{$self->{"bestCGM"}} = undef;
     }
   }
  return @{$self->{"bestCGM"}};
}

sub getPercentHSPsCoveredByCGM{
  my $self = shift;
  return int((scalar($self->getBestCGM())/$self->countGoodHSPs())*100);
}

sub getCGMcoverage{
  my($self,@cgm) = @_;
  return int((($cgm[-1]->getSubjectEnd() - $cgm[0]->getSubjectStart())/$self->getLength())*100);
}

sub getMaxForwardCoverage{
  my($self) = @_;
  my @s = $self->getGoodHSPsBySubjectStart();
  return int((($s[-1]->getSubjectEnd() - $s[0]->getSubjectStart())/$self->getLength())*100)
}

##returns the percentof maxForward coverage that it covered by cgm
sub getPercentCGMCoverage{
  my($self,@cgm) = @_;
  return int(($self->getCGMcoverage(@cgm) / $self->getMaxForwardCoverage())*100);
}

sub getConsistentParents{
  my($self,$hsp,$overlap) = @_;
  my@par;
  foreach my $h ($self->getGoodHSPsBySubjectStart()){
    push(@par,$h) if (abs($hsp->getSubjectStart() - $h->getSubjectEnd()) <= $overlap && $self->checkOrderConsistency($hsp,$h) == 1);
  }
  return @par;
}

sub toString{
  my( $self,$long) = shift;
  my $ret = "";
  $ret = "Subject=".$self->getID().", Length=".$self->getLength().", Num HSPs=".$self->countHSPs()."\n";

  if($long){
    $ret .= ($self->checkOrderConsistency() ? " Order Consistent," : " Order InConsistent,").($self->checkOrientationConsistency() ? " Orientation Consistent," : " Orientation InConsistent,")."\n";
    $ret .= "  Containments: Query=".$self->countContainedQueryHSPs().", Subject=".$self->countContainedSubjectHSPs()."\n";
    my($subOver,$subGaps) = $self->checkSubjectGaps();
    $ret .= "  Subject Overlaps = $subOver, SubjectGaps = $subGaps\n";
    my($queOver,$queGaps) = $self->checkQueryGaps();
    $ret .= "  Query Overlaps = $queOver, Query Gaps = $queGaps\n";
  }
  $ret .= $self->hspsToString($self->getHSPsBySubjectStart());
  return $ret;
}

sub hspsToString {
  my($self,@hsps) = @_;
  my $ret = "";
  foreach my $h (@hsps){
    $ret .= $h->toString("    ");
  }
  return $ret;
}

##following is a snippet that could be used to check for Ns given thequery sequence..
##  my $querySeq = shift;  ##if pass in query sequence then check to see if blocked
##if the gap between HSPs is same for subject and query
## need to take into account orientation of query..NO.. ??

##  $querySeq =~ s/\s+//g; ##remove newlines

##if the difference is less than 5 bp check to see if blocked!!.
##what do I want to do here??  check if blocked??
##if difference is very similar then is likely just a region
##of poor quality sequence where there is not sufficient identityfor match

##	if($querySeq ne ""){  ##have a query sequence
##	  my $seq = substr($querySeq,$qs[$i]->getQueryEnd(),$queryGap);
##	  my $countN = 0;
##	  for(my $c=0;$c<length($seq);$c++){
##	    $countN++ if substr($seq,$c,1) =~ /N/i;
##	  }
##	  if($countN / $queryGap < 0.5){ ##less than half the residues are N's 
##	    $countGaps++ if $queryGap > 20;  ##count the gap..if > half are N's then don't
##	  }else{
##	    print STDERR "The gap \($queryGap bp\) is >50% Ns\n" if $debug == 1;
##	  }
##	}	

1;
