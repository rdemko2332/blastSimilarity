#!/usr/bin/perl

## perl module for dealing with blast results...
## contains the hits for one query sequence
## will parse both  --noseqs and normal blast.
## does not parse out the alignments, rather meta information for each subject and hsp

## Brian Brunk 5/17/99

package BlastAnal;

use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(new parseBlast setQueryLength getQueryLength sortSubjectsByMinQueryStart sortSubjectsByMaxQueryEnd sortSubjectsByPValue detectChimera getChimeras findChimeraBreakpoints buildChimera findSubjectsSpanningChimeraBreakpoint detectRepeats findRepeats getQuery5endMis getQuery3endMis getCoveredBySubjects analyzeSubjectsBySubject getStats getSummaryStats printStats addSubject getSubjectCount getSubjects getQueryID printThreshHoldHits printGeneModels compareGeneModelWithTotCons);

use Subject;
use HSP;
use Chimera;
use Repeat;

use strict;

my $debug = 0;

##simple constructor...need to make more robust if inheriting for other purposes
sub new {
  my $self = {};
  my $class = shift;
  my $d = shift;
  bless $self,$class;
  $debug = defined $d ? $d : 0;
  print STDERR "DEBUGGING ON\n" if $debug == 1;
  $self->{"countSubjects"} = 0;
  return $self;
}




##DataStructure returned by parseBlast...must have -noseqs on commandline
##new datastructure....
##{queryLength},{queryName}
##{subjects} = array of Subject objects...

## takes in a blast result as an array and builds the entire datastructure
sub parseBlast{
  system("touch blastAnal.log");
  open(LOG,">blastAnal.log");
  my $self = shift;
  my($minLength,$minPercent,$minPvalue,$regex,$blast,$remMaskedFromLen,$rpsblast,$minPercentLength) = @_;
  print LOG "BlastAnal->parseBlast($minLength,$minPercent,$minPvalue,$regex,$blast,$remMaskedFromLen,$rpsblast,$minPercentLength)\n";
  $self->{minLength} = $minLength;
  $self->{minPercent} = $minPercent;
  $self->{minPvalue} = $minPvalue;
  $minPercentLength =  $minPercentLength ? $minPercentLength / 100 : .0001;
  $self->{minPercentLength} = $minPercentLength;
  my %blast;
  my ($strand,$dir,$qStart,$qEnd,$sStart,$sEnd,$pValue,$matchLength,$matchPercent);
  my $sbjct;
  my $queryMatch;
  my $haveQStart = 0;
  my $haveSStart = 0;
  my $parseType;
  my $frame;
  my $score;
  my $desc;
  my $inDesc = 0;
  my $identities;
  my $positives;
  die "these are the args to parseBlast($minLength,$minPercent,$minPvalue,$regex,$rpsblast,$minPercentLength)\n" if scalar(@{$blast}) <= 10;
  foreach (@{$blast}) {
    #    print STDERR $_;
    if (/^(\S*BLAST\S+)/){   ##gets the algorighm used for the query...
      my $algo = $1;
      print LOG "algorithm = $algo\n";
      if($algo =~ /BLASTX/i){
        $parseType = 1;
      }else{
        $parseType = 0;
      }
      print LOG "ParseType='$parseType'\n";
    }
    if (/^Query=\s*(\S+)/) {
      $self->{queryName} = $1;	##query name/id
    } elsif (/^Length=(\S+)/) {
	print LOG "Letters Trigger\n";
	$self->setQueryLength($1); ##query length
	print LOG "$1\n";
    }
    if (/^\>\s*(\S+)/ || /ctxfactor/ || /^Lambda/) { ##ctxfactor catches the last one...Lambda last of rpsblast
      my $sbjctId;
      $desc = "";
      if (/^\>\s*(\S+)/){
        $sbjctId = $1 ? $1 : $2;
      }else{
        print "$self->{queryName}: Unable to match subject using regex '$regex' on:\n  $_" unless (/ctxfactor/ || /^Lambda/);
      }

      print LOG "Matching subjectId = $sbjctId\n";
      ##      print STDERR "REGEX matched $_";
      if ($haveQStart) {
        print "Have last QStart:", $sbjct->getID()," SS:$sStart,SE:$sEnd, QS:$qStart, QE:$qEnd, Length=$matchLength, Percent=$matchPercent, pValue=$pValue, Frame=$frame\n" if $debug == 1 && $sbjct;
        ##have the complete HSP...
        ##want to add only if meets minimum reqs....do this on an HSP basis...
        #	print $sbjct->getID()," $sStart,$sEnd: Length=$matchLength, Percent=$matchPercent, pValue=$pValue\n";
        if ($matchLength >= $minLength && $matchPercent >= $minPercent && $pValue <= $minPvalue) {
          print LOG "Match meets reqs\n";
          if($remMaskedFromLen){  ##want to remove Xs from the match length....
#            print STDERR "removing X from match\n";
#            print STDERR "Before: $queryMatch\n";
            $queryMatch =~ s/X//g;
#            print STDERR "After: $queryMatch\n";
            if(length($queryMatch) < 10){ ##don't remove if final length < 10
              print "RemoveMaskedResiduesFromLength: error...length = ".length($queryMatch)." following removal of X's\n";
            }else{
#              print "match before $matchLength\n";
              $matchLength = length($queryMatch);
#              print "match after $matchLength\n";
            }
          }
          $sbjct->addHSP($qStart,$qEnd,$sStart,$sEnd,$matchLength,$matchPercent, $pValue, $dir, $frame, $score,$identities,$positives) if $sbjct;
        } 
        $haveQStart = 0;				##reset for next one...
        $haveSStart = 0;
      } 
      my $shortSeq = defined $sbjct ? $sbjct->getQueryLength() < $sbjct->getLength() ? $sbjct->getQueryLength() : $sbjct->getLength() : 0;
      print LOG "QueryLength: ".$sbjct->getQueryLength().", subjectLength: ".$sbjct->getLength().", shortest sequence is $shortSeq long\n" if defined $sbjct;
      if (defined $sbjct && $sbjct->countHSPs() >= 1 && $sbjct->getTotalHSPLength() > $minPercentLength * $shortSeq ) {
        $self->addSubject($sbjct);
      }


      ##return heree if end...
      return if (/ctxfactor/ || /^Lambda/);
      
      ##      print STDERR "New Subject $sbjctId\n";
      $sbjct = Subject->new($sbjctId) if $sbjctId;
      print LOG "DING\n" if $sbjctId;
      ##lets get the description here could be on multiple lines so need to set something to take care of this..
      if(/^\>\S+\s(.*)$/){
        $desc = $1;
        $inDesc = 1;
        next;
      }
      
    }
    if (/^Length=(\S+)/) {
      if($sbjct){
        print LOG "SubLen Trigger\n";	    
        $inDesc = 0;
        $sbjct->setDescription($desc);  ##set the description
        my $sbjctLength = $1;
        $sbjctLength =~ s/,//g;
        $sbjct->setLength($sbjctLength);
        $sbjct->setQueryLength($self->{"queryLength"});	##set the query length for use in 
      }
      ##end remaining calculations
    }
    if (/^\s*Score\s*=\s*(\d+).*?=\s(\S+)/ || ($rpsblast && /^\s*Score\s*=\s*\d+\sbits\s\((\d+).*=\s(\S+)$/)) {
      my $tmpScore = $1;
      my $tmpValue = $2;
      $tmpValue =~ s/,//g;

      ##need to add here if $haveQStart != 0;
      if ($haveQStart) {
        ##have the complete HSP...
        ##want to add only if meets minimum reqs....do this on an HSP basis...
        if ($matchLength >= $minLength && $matchPercent >= $minPercent && $pValue <= $minPvalue && $sbjct) {
#        if (($minPercentLength ? ( $matchLength >= ($sbjct->getLength() < $sbjct->getQueryLength() ? $sbjct->getLength() * $minPercentLength : $sbjct->getQueryLength() * $minPercentLength))  : $matchLength >= $minLength) && $matchPercent >= $minPercent && $pValue <= $minPvalue && $sbjct) {
          print STDERR "Adding Sbjct (score) ",$sbjct->getID()," SS:$sStart,SE:$sEnd, QS:$qStart, QE:$qEnd, Length=$matchLength, Percent=$matchPercent, pValue=$pValue, Frame=$frame\n" if $debug == 1;
          print STDERR "Adding Sbjct\n" if $debug;
          if($remMaskedFromLen){  ##want to remove Xs from the match length....
#            print STDERR "removing X from match\n";
#            print STDERR "Before: $queryMatch\n";
            $queryMatch =~ s/X//g;
#            print STDERR "After: $queryMatch\n";
            if(length($queryMatch) < 10){ ##don't remove if final length < 10
              print STDERR "RemoveMaskedResiduesFromLength: error...length = ".length($queryMatch)." following removal of X's\n";
            }else{
#              print "match before $matchLength\n";
              $matchLength = length($queryMatch);
#              print "match after $matchLength\n";
            }
          }
          $sbjct->addHSP($qStart,$qEnd,$sStart,$sEnd,$matchLength,$matchPercent, $pValue, $dir, $frame, $score,$identities,$positives); 
        }
        $haveQStart = 0;				##reset for next one...
        $haveSStart = 0;
      }
      $pValue = $tmpValue;
      $score = $tmpScore;
    }
    ##get the strand in the Identities line as the following messes things up!!!
    ##get the strand...will work for all blast types...
#    if (/^\s*(Minus|Plus)\sStrand/){
#      $strand = $1;
#    }
    ##following specific to blastx
    if ($parseType && /^\s*Identities\s=\s(\d+)\/(\d+)\s\((\d+)\%\),\sPositives\s=\s(\d+).*Frame\s=\s.*(.\d)$/) {
      $identities = $1; $matchLength = $2; $matchPercent = $3; $positives = $4; $frame = $5;
      $dir = $frame =~ /^\+/ ? 1 : 0;
      ##following specific to blastn
    }elsif (/^\s*Identities\s=\s(\d+)\/(\d+)\s\((\d+)\%\),\sPositives\s=\s(\d+).*Strand\s=\s(\w+)/) {
      $identities = $1; $matchLength = $2; $matchPercent = $3; $positives = $4; $strand = $5;
      $dir = $strand eq 'Minus' ? 0 : 1;
      ##following for blastp the default if others not matched...
    }elsif (/^\s*Identities\s=\s(\d+)\/(\d+)\s\((\d+)\%\),\sPositives\s=\s(\d+)/){
      $identities = $1; $matchLength = $2; $matchPercent = $3; $positives = $4;
      $dir = 1;  ##always 1 for blastp which is what this catches...
    }elsif (/^\s*Frame\s=\s(.\d)/){ ##frame for rpsblast
      $frame = $1;
      $dir = $frame =~ /^\+/ ? 1 : 0;
    }

    if (/^Query\s+(\d+)\s(.*)\s(\d+)\s*$/) {
      print LOG "Matching query: $1 $3\n";
      if ($haveQStart == 0) {
        $queryMatch = "";
        $qStart = $1; 
        $haveQStart = 1;
      }
      $qEnd = $3;
      my $tmpQMatch = $2;
      $tmpQMatch =~ s/\s//g;  ##gets rid of any spaces
      $queryMatch .= $tmpQMatch;
    } elsif (/^Sbjct\s+(\d+)\s.*\s(\d+)\s*$/) {
      if ($haveSStart == 0) {
        $sStart = $1; 
        $haveSStart = 1;
      }
      $sEnd = $2;
    }
    $desc .= $_ if $inDesc;  ##append to description if in the description
  }
  return;
}  

sub setQueryLength {
  my($self,$l) = @_;
  $l =~ s/,//g;
  $self->{queryLength} = $l;
}

sub getQueryLength {
  my($self) = @_;
  return $self->{queryLength};
}

sub sortSubjectsByMinQueryStart{
  my $self = shift;
  if (!exists $self->{"sortByQStart"}) {
    @{$self->{"sortByQStart"}} = sort{$a->getMinQueryStart() <=> $b->getMinQueryStart()} $self->getSubjects();
  }
  return @{$self->{"sortByQStart"}};
}

sub sortSubjectsByMaxQueryEnd{
  my $self = shift;
  if (!exists $self->{"sortByQEnd"}) {
    @{$self->{"sortByQEnd"}} = sort{$a->getMaxQueryEnd() <=> $b->getMaxQueryEnd()} $self->getSubjects();
  }
  return @{$self->{"sortByQEnd"}};
}

sub sortSubjectsByPValue{
  my $self = shift;
  if (!exists $self->{"sortByPValue"}) {
    @{$self->{"sortByPValue"}} = sort{$a->getPValue() <=> $b->getPValue()} $self->getSubjects();
  }
  return @{$self->{"sortByPValue"}};
}

##detecting chimeras
############################################################
## Breakpoint = 0          * OK
## Breakpoint = 1          > [Ch, AS, OK, !]
##   Dangling Ends = none  > [Ch, AS, OK] 
##     Spanned = 1         * AS, OK
##     Spanned = 0         > [Ch, AS, OK]
##       Deep = 1          * Ch
##       Deep = 0          * Ch, AS, OK (review)
##   Dangling Ends = uni   > [Ch, As]
##     Spanned = 1         > [Ch, AS]
##       QvsSub = 1        * Ch
##       QvsSub = 0        * AS
##     Spanned = 0         > [Ch, AS]
##       Deep = 1          * Ch
##       Deep = 0          * Ch, AS (review)
##   Dangling Ends = bi    > [Ch, !]
##     Spanned = 1         > [Ch, !]
##       QvsSub = 1        * Ch
##       QvsSub = 0        * ! (review) ##not sure about this one....should it be OK
##     Spanned = 0         * Ch
############################################################    

##return values: 0=OK, 1=chimera, 2=review no dangling, 3=review 1 dangling, 4=review 2 dangling ends.
##add to ichimera obj...
sub detectChimera{
  my($self,$slop) = @_;					##$slop is the amount of slop in the breakpoint de to blast
  ##loop through each possible chimera breakpoint and test if chimera
  my $qvssub = 1;
  my @span;
  foreach my $c ($self->findChimeraBreakpoints($slop)) {
    
    @span = $self->findSubjectsSpanningChimeraBreakpoint($c->getBreakpoint(),$slop);
    $c->{"numberSpanned"} = scalar(@span);

    ##first dangling == 0
    if ($c->getDangling($slop) == 0) { ##dangle by at least $slop
      if (scalar(@span) > 0) {
        $c->{"isChimera"} = 0;	##spans and no dangling ends therefore OK
      } else {
        if ($c->getDepth() == 1) {
          $c->{"isChimera"} = 1;
        } else {
          $c->{"isChimera"} = 2; ##will requiere review
        }
      }
      
      ##dangling on only one end...
    } elsif ($c->getDangling($slop) == 1) {
      if (scalar(@span) > 0) {
        $qvssub = 1;						##are very similar....looking for one that isn't...
        foreach my $sp (@span) {
          if ($sp->compareQueryVsSubject() == 0) {
            ##good subject spanning....looks lke should be alternative splicing
            $c->{"isChimera"} = 0;
            $qvssub = 0;
          }
        }
        $c->{"isChimera"} = 1 if $qvssub == 1; 
      } else {
        if ($c->getDepth() == 1) {
          $c->{"isChimera"} = 1;
        } else {								#will need to review....
          $c->{"isChimera"} = 3;
        }
      }
      
      ##dangling on two ends...
    } elsif ($c->getDangling($slop) == 2) {
      
      if (scalar(@span) > 0) {
        $qvssub = 1;						##are very similar....looking for one that isn't...
        foreach my $sp (@span) {
          if ($sp->compareQueryVsSubject() == 0) {
            ##unclear here....need to review:
            $c->{"isChimera"} = 4;
            $qvssub = 0;
            ##            last;
          }
        }
        $c->{"isChimera"} = 1 if $qvssub == 1; 
      } else {
        $c->{"isChimera"} = 1;  ##is a chimera
      }
      
    } else {  
      print "ERROR: getDangling returned invalid value\n";
    }
  }
  ##now need to go through each of chimeras and output results...
  my $maxChim = 0;
  my $isChimera = 0;
  my $errString = "";
  if (exists $self->{"chimeras"}) {
    ##    $errString = "Chimeras for $self->{'queryName'}\n";
    my $cnt = 1;
    foreach my $chim (@{$self->{"chimeras"}}) {
      last if $chim eq "";			##this is just a holder that indicates that the findChimera
      ##method has been run and there are not chimeras...
      $maxChim = $chim->{'isChimera'} if $chim->{'isChimera'} > $maxChim;
      $isChimera = 1 if $chim->{'isChimera'} == 1;
      $errString .= "  Chimera $cnt: $chim->{'isChimera'}, danglingEnds=".$chim->getDangling($slop).", numberSpanned=$chim->{'numberSpanned'}, depth=".$chim->getDepth().", numberLeft=".$chim->countLeft().", numRight=".$chim->countRight().", breakpoint=".$chim->getBreakpoint()."\n";
      #      $errString .= $chim->toString();
      $cnt++;
    }
  }
  return ($isChimera,$maxChim,$errString);
}

sub getChimeras {
  my($self,$slop) = @_;
  if(!exists $self->{'chimeras'}){
    $self->detectChimera($slop);
  }
  return @{$self->{'chimeras'}};
}

## returns array of Chimera objs
sub findChimeraBreakpoints{
  my($self,$slop) = @_;
  if ($self->{'findChimera'} != 1) {
    my @l = $self->sortSubjectsByMaxQueryEnd();
    my @r = reverse($self->sortSubjectsByMinQueryStart()); ##want right one reversed
    my $lindex = 0;
    my $rindex = 0;
    for (my $li = $lindex;$li < scalar(@l);$li++) {
      last if $l[$li]->getMaxQueryEnd() - $r[$rindex]->getMinQueryStart() > $slop;
      for (my $ri = $rindex;$ri < scalar(@r);$ri++) {
        last if $l[$li]->getMaxQueryEnd() - $r[$ri]->getMinQueryStart() > $slop;
        if (abs($l[$li]->getMaxQueryEnd() - $r[$ri]->getMinQueryStart()) <= $slop) {
          ##create new chimera object and pass arrays into extendBrkpt()
          ($li,$ri) = $self->buildChimera($li,$ri,$slop); ##pass back the index of last l and r
          last;
        }
      }
    }
    $self->{'findChimera'} = 1;	##indicates that have done the analysis
  }
  if (!exists $self->{'chimeras'}) {
    return;
  }
  return @{$self->{"chimeras"}};
}

sub buildChimera{
  my($self,$lindex,$rindex,$slop) = @_;
  my @l = $self->sortSubjectsByMaxQueryEnd();
  my @r = reverse($self->sortSubjectsByMinQueryStart()); ##want right one reversed
  my $chim = Chimera->new();
  ##now add the current and any others to the chimera
  ## the first time in loop will add the current one to $chim
  ##first do left
  my $ledge = $l[$lindex]->getMaxQueryEnd();
  for (my $li = $lindex;$li < scalar(@l);$li++) {
    if ($l[$li]->getMaxQueryEnd() - $ledge <= $slop) {
      $chim->addLeftSubject($l[$li]);
    } else {
      $lindex = $li;
      last;
    }
    $lindex = $li;
  }
  ##then do right
  my $redge = $r[$rindex]->getMinQueryStart();
  for (my $ri = $rindex;$ri < scalar(@r);$ri++) {
    if ($redge - $r[$ri]->getMinQueryStart() <= $slop) {
      $chim->addRightSubject($r[$ri]); 
    } else {
      $rindex = $ri;
      last;
    }
    $rindex = $ri;
  }
  push(@{$self->{"chimeras"}},$chim);
  return($lindex,$rindex);
}

##returns array of subjects
##$slop is amount by which it needs to extend beyond the brkpt.
sub findSubjectsSpanningChimeraBreakpoint{
  my($self,$loc,$slop) = @_;
  my @span;
  foreach my $s ($self->sortSubjectsByMinQueryStart()) {
    last if $s->getMinQueryStart() > $loc - $slop;
    if ($s->getMaxQueryEnd() > $loc + $slop) {
      ##      print "Spans $loc: ",$s->getMinQueryStart(), "-",$s->getMaxQueryEnd(),"\n";
      push(@span,$s);
    }
  }
  return @span;
}

########################################
## Detecting unblocked repeats...
## If No stack - OK
## 
## If Stack
## 
##    If No overhanging ends in query then OK (query is subsumed by subjects)
##       -will block the query sequences and throw out sequences that do not have
##             >50 bp of un-blocked sequence
## 
##    If Overhanging ends in query (may be only one if repeat is at end)
## 
##       If 1 HSP then likely REPEAT => ignore all matches that don't extend
##          outside this region.
##          
##          Overhanging ends in subjects = 1 => definitely REPEAT
## 
##          Overhanging ends in subjects = 0 => ??very unlikely scenario..not clear
## 
##       If >1 HSP (could be multiple copies of same repeat or perhaps alt. splicing
##         OR could it be a valid match that also has a repeat in it...no should in
##         this case match entirely and have just 1 HSP)
## 
##          If subject locations overlap/superimpose then REPEAT 
##             =ignore all matches that do not extend outside these locations
##                     between them would be "outside" the locations
##          
##          If subject locations distinct then likely alternative splicing
##             OR alternatively the query and subject could contain two distinct
##               different repeats (I will be blocking for repeats so the likelihood
##               that there are two different repeats that are both novel is small).
## 
##             ?do I want to do further analysis to detemine consistency with
##                 alternative splicing model before ruling OK...?
## ########################################

##this method is the call to find repeats...should return an array of subjects to
##ignore...

## find 5' edge (minQueryStart) and gather up then check to see if there is a 
## correlation with 3' (maxQueryEnd)
## edge...if so is potential repeat and need to check the above rules..

#what to return here....region to ignore when determining good matches...

sub detectRepeats {
  my($self,$slop) = @_;

  if ($self->getSubjectCount() < 5 ) { ##need at least 5 to detecxt repeat...
    return;
  }

  my @rep = $self->findRepeats($slop); ##array of repeats....

  if ($debug == 1) {
    print STDERR "detectRepeats returns ",scalar(@rep), " repeats:\n";
    foreach my $r (@rep) {
      print STDERR $r->toString(1),"\n";
    }
  }

  ##want to merge repeats if they are overlapping...
  my @sort = sort{$a->getRepeatStart() <=> $b->getRepeatStart()} @rep;
  undef @rep;										##use for merged ones...
  my $merge = 0;
  my $nr;
  for (my $i = 0;$i<scalar(@sort) - 1;$i++) {
    if ($sort[$i+1]->getRepeatStart() < $sort[$i]->getRepeatEnd()) {
      ##are overlapping...merge
      if ($merge == 0) {
        $nr = Repeat->new();
        $merge = 1;
      }
      foreach my $s ($sort[$i]->getSubjects()) {
        $nr->addSubject($s);
      }
    } elsif ($merge == 1) {
      $merge = 0;
      foreach my $s ($sort[$i]->getSubjects()) {
        $nr->addSubject($s);
      }
      push(@rep,$nr);
    } else {
      ##does not overlap so just add to rep
      push(@rep,$sort[$i]);
    }
  }
  ##need to do last if $merge == 1;
  if ($merge == 1) {						##need to add the i+1 (last) case
    foreach my $s ($sort[-1]->getSubjects()) {
      $nr->addSubject($s);
    }
    push(@rep,$nr);
  }
  ##repeats merged and should be non-overlapping

  if ($debug == 1) {
    print STDERR "detectRepeats returns ",scalar(@rep), " repeats:\n";
    foreach my $r (@rep) {
      print STDERR $r->toString(1),"\n";
    }
  }
  return @rep;
}

##finds all repeat stacks..
# populates array of subjects that are repeats ... @{$self->{'repeats'}}

sub findRepeats{
  my($self,$slop) = @_;
  my $rep;
  if (!exists $self->{'repeats'}) {
    my $firstElement = 1;
    my $inRepeat = 0;
    my @s = $self->sortSubjectsByMinQueryStart();
    ###NOTE:  want to only take repeat from a single point rather than each minQuery start
    my $lEdge = $s[0]->getMinQueryStart(); ##this is the first edge
    for (my  $i = 0;$i < scalar(@s) - 1;$i++) { 
      if ($s[$i+1]->getMinQueryStart() - $lEdge <= $slop) {
        if ($firstElement == 1) {
          $rep = Repeat->new(); ##create new Repeat object..
          $firstElement = 0;
          $inRepeat = 1;
        }
        print STDERR "REPEAT: adding subject ",$s[$i]->getID(),"\n" if $debug == 1;
        $rep->addSubject($s[$i]);
        if ($i == scalar(@s) - 2) {	##need to add next one as am in last loop
          $rep->addSubject($s[$i+1]);
          my $verified = $rep->verifyRepeat($slop); ##does all the verification that this is a repeat
          push(@{$self->{'repeats'}},$rep) if $verified == 1;
        }
        
      } elsif ($inRepeat == 1) { ##need this to get the $i+1 as am working one ahead
        $lEdge = $s[$i]->getMinQueryStart(); ##new left edge
        $firstElement = 1;			##now looking for the first element
        $inRepeat = 0;
        ##add the last one nd then verify repeats....Repeat method..
        ##if really a repeat then add to @$self->{'repeats'}
        $rep->addSubject($s[$i]);
        print STDERR "Repeat contains ",$rep->countSubjects(), " subjects\n" if $debug == 1;
        $rep->verifyRepeat($slop); ##does all the verification that this is a repeat
        push(@{$self->{'repeats'}},$rep) if($rep->countSubjects() > 1 && $rep->getVerified() == 1);
      } else {
        $lEdge = $s[$i]->getMinQueryStart(); ##new left edge
        $firstElement = 1;			##now looking for the first element
      }
    }
    
    if (!exists $self->{'repeats'}) {	##there aren't any....make 0
      $self->{'repeats'}->[0] = "";
    }
  }
  if ($self->{'repeats'}->[0] eq "") {
    return;
  }
  return @{$self->{'repeats'}};
}

##NOTE:  the repeat object is initially populated with things that have a common
## left edge.  It is the job of the Repeat to then verify that this constitutes
## a repeat and remove any subjects that don't meet ehe rules for a repeat stack above
#3will need a 'verified' bit to know if the repeat bj is verified...else verify
## before doing anyhing else.


sub getQuery5endMis{
  my $self = shift;
  if (!exists $self->{"query5pMis"}) {
    $self->{"query5pMis"} = 10000; ##arbritrarily large as want smallest
    foreach my $s ($self->getSubjects()) {
      if ($self->{"query5pMis"} > $s->getQuery5endMis()) {
        $self->{"query5pMis"} = $s->getQuery5endMis();
      } 
    }
  }
  return $self->{"query5pMis"};
}

sub getQuery3endMis{
  my $self = shift;
  if (!exists $self->{"query3pMis"}) {
    $self->{"query3pMis"} = 10000; ##arbritrarily large as want smallest
    foreach my $s ($self->getSubjects()) {
      if ($self->{"query3pMis"} > $s->getQuery3endMis()) {
        $self->{"query3pMis"} = $s->getQuery3endMis();
      } 
    }
  }
  return $self->{"query3pMis"};
}

##returns 0 if leftover at ends and 1 if covered to ends
sub getCoveredBySubjects{
  my $self = shift;
  my $left = shift;							##amount to be leftover....
  if (! $left) { $left = 20; }
  if ($self->getQuery5endMis() <= $left && $self->getQuery3endMis() <= $left) {
    return 1;
  }
  return 0;
}

sub analyzeSubjectsBySubject{
  my $self = shift;
  return if $self->getSubjectCount() == 0;
  my @singletons;
  my %multiples;
  my $haveImmuno = 0;						##keeps track if is immunoglobulin gene...
  my $order = 0;								##order
  my $orientation = 0;					##orientation
  my $subOver = 0;							##number of subject overlaps
  my $subGaps = 0;							##number of subject Gaps
  my @multiples;
  foreach my $s ($self->getSubjects()) {
    my $id = $s->getID();
    #    if($id <= 100009){
    #      $haveImmuno = 1;
    #      next;  ##don't want to include these in the analysis
    #    }
    if ($s->countHSPs() == 1) {	##don't print out singletons HSPS
      push(@singletons,$id);
      next;
    }
    push(@multiples,$id);
    print "\n\>DT.",$id,", Length=",$s->getLength(),", Num HSPs=",$s->countHSPs(),($s->checkConsistency() ? " Order and Orient Consistent" : " Order or Orient InConsistent"),"\n" if $debug == 1;
    ($subOver,$subGaps) = $s->checkSubjectGaps();
    print "\tSubject Overlaps = $subOver, SubjectGaps = $subGaps\n" if $debug == 1;
    $order = $s->checkOrderConsistency();
    $orientation = $s->checkOrientationConsistency();
    print STDERR $s->toString() if $order + $orientation != 2;
    $multiples{$id} = [$s->countHSPs(),$s->countContainedQueryHSPs(),$order,$orientation,$subOver == 0 ? 1 : 0,$subGaps == 0 ? 1 : 0]; ##this is the data for multiples!
  }
  print "\nSingletons: ",scalar(@singletons),"\n" if $debug == 1;

  #  if($haveImmuno == 1){print "MATCHES immunoglobulin genes\n";}
  return(\@singletons,\%multiples,$haveImmuno);
}

sub getStats{
  my $self = shift;
  my $minNum = shift;						##minimum number of HSPs
  my $ret;
  $minNum = $minNum ? $minNum : 1;
  $ret = "QueryName = ".$self->{queryName}.", QueryLength = ".$self->{queryLength}."\n";
  if ($self->{"countSubjects"} == 0) {
    $ret = "\tNo Hits\n";
    return $ret;
  }
  foreach my $s ($self->getSubjects()) {
    next if $s->countHSPs() < $minNum; ##don't print out singletons HSPS
    $ret .= $s->toString();
  }
  return $ret;
}

sub getSummaryStats {
  my($self,$doDesc) = @_;
  my $ret = "QueryName = ".$self->{queryName}.", QueryLength = ".$self->{queryLength}."\n";
  if ($self->{"countSubjects"} == 0) {
    $ret .= "\tNo Hits\n";
    return $ret;
  }
  foreach my $s ($self->getSubjects()) {
    $ret .= "  ".$s->getSummaryStats()."\n";
    $ret .= $s->getDescription(80,4)."\n" if $doDesc;
  }
  return $ret;
}

sub printStats {
  my($self,$minNum) = @_;
  print $self->getStats($minNum);
}

sub addSubject{
  my($self,$sbjct) = @_;
  push(@{$self->{subjects}},$sbjct);
  $self->{"countSubjects"}++;
}

sub getSubjectCount{
  my $self = shift;
  return $self->{"countSubjects"};
}

sub getSubjects{
  my $self = shift;
  return if $self->{"countSubjects"} == 0;
  return @{$self->{subjects}};
}

sub getQueryID{
  my $self = shift;
  return $self->{"queryName"};
}

sub printThreshHoldHits {
  my $self = shift;
  return if $self->{"countSubjects"} == 0;
  foreach my $s ($self->getSubjects()) {
    print $self->getQueryID(), "\t",$s->getID(),"\n";
  }
}

sub printGeneModels{
  my $self = shift;
  my @cgm;
  print "QueryName = $self->{queryName}, QueryLength = $self->{queryLength}\n";
  if ($self->{"countSubjects"} == 0) {
    print "\tNo Hits\n";
    return;
  }
  foreach my $s ($self->getSubjects()) {
    print "\n",$s->toString();
    @cgm = $s->getConsistentGeneModel();
    print "  CONSISTENT GENE MODEL: ",($s->checkGoodCGM() ? "Yes" : "No"), ", ",scalar(@cgm), " HSPs\n";
    foreach my $h (@cgm) {
      print "    ",$h->toString();
    }
  }
}

sub compareGeneModelWithTotCons{
  my $self = shift;
  #  print "\nQueryName = $self->{queryName}, QueryLength = $self->{queryLength}\n";
  my $cgmgt1 = 0;
  my $corCGM = 0;
  my $totCor = 0;
  my $totSubs = 0;
  my @vals;
  foreach my $s ($self->getSubjects()) {
    next if $s->countHSPs() < 2; ##don't do singletons
    $totSubs++;
    #    print "\n",$s->toString();
    #    print $s->getID(),": ",$s->countHSPs()," HSPs\n";
    #    print "  CONSISTENT GENE MODEL: ",($s->checkGoodCGM() ? "Yes" : "No"), ", ",scalar($s->getConsistentGeneModel()), " HSPs\n";
    my($subOver,$subGaps) = $s->checkSubjectGaps();
    @vals = ($s->checkOrderConsistency(),$s->checkOrientationConsistency(),$subOver == 0 ? 1 : 0,$subGaps == 0 ? 1 : 0);
    #    print "  Old Analysis: Order=$vals[0], Orient=$vals[1], Over=$vals[2], Gaps=$vals[3]\n";
    $totCor++ if $vals[0]+$vals[1]+$vals[2]+$vals[3] == 4;
    $corCGM += $s->checkGoodCGM();
    $cgmgt1++ if scalar($s->getConsistentGeneModel()) > 1;
  }
  return($totSubs,$corCGM,$totCor,$cgmgt1);
}


1;
