#!/usr/bin/perl

## module for storing subject data and analyzing potential repeats

## a Repeat is an array of subjects that stack with ends that are distinct


## Brian Brunk 11/6/99

package Repeat;

use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(addSubject getSubjects resetSubjects getSubjectsByQueryStart getSubjectsByQueryEnd countSubjects verifyRepeat  doSubjectsSuperimpose countAveHSPs checkQueryEnds checkSubjectEnds pruneToTallestStack getVerified  getRepeatBorders getRepeatStart getRepeatEnd getRepeatLength toString);

use strict;

my $debug = 0;  ##flip debugging on and off...

sub new {
  my($class) = @_;  
  my $self = {};
  bless $self, $class;
  $self->{'countSubjects'} = 0;
  return $self;
}

##store subjects as hash so can easily delete...
sub addSubject{
  my($self,$sbj) = @_;
  push(@{$self->{"subjects"}}, $sbj);
  $self->{'countSubjects'}++;
}

sub getSubjects{
  my($self) = @_;
  return @{$self->{"subjects"}};
}

sub resetSubjects{
  my $self = shift;
  undef $self->{'subjects'};
  $self->{'countSubjects'} = 0;
}

sub getSubjectsByQueryStart {
  my$self = shift;
  if(!exists $self->{'subByStart'}){
    @{$self->{'subByStart'}} = sort{$a->getMinQueryStart() <=> $b->getMinQueryStart()}$self->getSubjects();
  }
  return @{$self->{'subByStart'}};
}

sub getSubjectsByQueryEnd {
  my$self = shift;
  if(!exists $self->{'subByEnd'}){
    @{$self->{'subByEnd'}} = sort{$a->getMaxQueryEnd() <=> $b->getMaxQueryEnd()}$self->getSubjects();
  }
  return @{$self->{'subByEnd'}};
}

sub countSubjects{
  my $self = shift;
  return $self->{'countSubjects'};
}

##verifies that this is a repeat using above rules...need to remove subjects
##that do not meet rules so that repeat just has one stack with well defined
##edges.  Should be quite tall and unlikely to generate more than one stack

## steps
## If Stack
## 1.  find tallest stack and remove those that do not stack
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

##set verified to 1 if is repeat else set to 0...that way can keep stack for
##further analysis if necessary

sub verifyRepeat{
  my($self,$slop) = @_;
  if(!exists $self->{'verified'}){  ##if verified set to 1 else set to 0
    $self->{'verified'} = 0;
    ##need to now build right end of stack..get the tallest one..
    $self->pruneToTallestStack($slop);
    if($self->countSubjects() == 0){
      $self->{'verified'} = 0;
      return 0;
    }
    ##check for overhanging ends in query...only need to check one!
    if($self->checkQueryEnds($slop) == 0){  ##has 0 overhanging ends...subsumed
      $self->{'verified'} = 0;
      return 0;
    }
    ##if makes it to here, has at least one overhanging query end...don't
    ## care whether one or two.
    ##check to see how many hsps...
    ##?should all be checked??
    if($self->countAveHSPs() == 1){ ##check to make certain that subject ends over
      if($self->checkSubjectEnds($slop) >= 1){  ##almost certainly a repeat stack
        $self->{'verified'} = 1;
        return 1;
      }else{  ##happens if subjects subsumed.....
        $self->{'verified'} = 0;
##        print STDERR "ERROR: one HSP but overhanging subject ends < 1\n";
      }
    }else{  ##there are more than one hsps...do subject hsps superimpose??
      if($self->doSubjectsSuperimpose($slop) == 1){
        $self->{'verified'} = 1;
        return 1;
      }else{
        ##could do further analyses here but is likely not a repeat....
        $self->{'verified'} = 0;
        return 0;
      }
    }
  }
  return $self->{'verified'};
}

##check each subject that has more than one HSP and look for the locations
##of each hsp to superimpose with the ends within $slop
##if one of the subjects has locations that superimpose then the repeat does!!
sub doSubjectsSuperimpose{
  my($self,$slop) = @_;
  my $sup = 0;
  my $haveOne = 0;
  if(!exists $self->{'subSuperImpose'}){
    $self->{'subSuperImpose'} = 0;
    foreach my $s ($self->getSubjects()){
      next if $s->countHSPs() == 1;
      my @h = $s->getHSPsBySubjectStart();
      for(my $i=0;$i<$s->countHSPs() - 1;$i++){
        if($h[$i+1]->getSubjectStart() - $h[$i]->getSubjectStart() <= $slop && abs($h[$i+1]->getSubjectEnd() - $h[$i]->getSubjectEnd()) <= $slop){
          $self->{'subSuperImpose'} = 1;
          push(@{$self->{'superImposeExemplar'}}, $h[$i]);  ##this is the sequence that is the examplar that superimposes...use to get the locations for determining what to block out
          $haveOne = 1;
        }elsif($haveOne == 1){
          $haveOne = 0;
          push(@{$self->{'superImposeExemplar'}}, $h[$i]);  
        }
          
      }
      return 1 if (exists $self->{'superImposeExemplar'});
    }
  }
  return $self->{'subSuperImpose'};
}

sub countAveHSPs{
  my $self = shift;
  my $tot = 0;
  
  foreach my $s ($self->getSubjects()){
    $tot += $s->countHSPs();
  }
  my $ave = $tot/$self->countSubjects();
  if($ave % 1 < 0.5){
    return int($ave);
  }else{
    return int($ave + 1);
  }
}

sub checkQueryEnds{
  my($self,$slop) = shift;
  if($self->countSubjects() > 0){
    print STDERR "checkQueryEnds: ",$self->{'subjects'}->[0]->getQueryEndsMissed($slop),"\n" if $debug == 1;
    return $self->{'subjects'}->[0]->getQueryEndsMissed($slop);
  }
}

##need to check each of the subjects for this one...and summarize
sub checkSubjectEnds{
  my($self,$slop) = shift;
  my $tot = 0;
  foreach my $s ($self->getSubjects()){
    $tot += $s->getSubjectEndsMissed($slop);
  }
  print STDERR "checkSubjectEnds: ",$tot / $self->countSubjects(),"\n" if $debug == 1;
  return $tot / $self->countSubjects();
}

##find the tallest stack (using right edge) and delete members not in stack
sub pruneToTallestStack{
  my($self,$slop) = @_;
  if($slop eq ""){die "USAGE: pruneToTallestStack\(\$slop\)\n";}
  print STDERR "pruning: \$slop = $slop\n" if $debug == 1;
  my %tmp;  ##for storing the arrays of potential subjects.
  my @s = $self->getSubjectsByQueryEnd();
  print STDERR "pruning: ",$self->countSubjects()," and ",scalar(@s)," in array\n" if $debug == 1;
  my $edge = $s[0]->getMaxQueryEnd();
  my $firstElement = 1;
  my $inRepeat = 0;
  my $cnt = 0;
  for(my $i = 0;$i < scalar(@s) - 1;$i++){
    if($s[$i+1]->getMaxQueryEnd() - $edge <= $slop){
      if($firstElement == 1){
        print STDERR "Have first element...\n" if $debug == 1;
        $cnt++;
        $firstElement = 0;
        $inRepeat = 1;
      }
      print STDERR "Adding: ",$s[$i]->getSummaryStats(),"\n" if $debug == 1;
      push(@{$tmp{$cnt}},$s[$i]);
      if($i == scalar(@s) - 2){  ##in last round...add last element
        print STDERR "Adding last one: ",$s[$i]->getSummaryStats(),"\n" if $debug == 1;
        push(@{$tmp{$cnt}},$s[$i+1]);
      } 
    }elsif($inRepeat == 1){
      $inRepeat = 0;
      $firstElement = 1; ##ready for first element
      $edge = $s[$i]->getMaxQueryEnd();  ##new edge
      ##add the last one nd then verify repeats....Repeat method..
      ##if really a repeat then add to @$self->{'repeats'}
      print STDERR "Adding last: ",$s[$i]->getSummaryStats(),"\n" if $debug == 1;
      push(@{$tmp{$cnt}},$s[$i]);
    }else{
      $firstElement = 1;
      $edge = $s[$i]->getMaxQueryEnd();  ##new left edge
    }
  } 
  ##now get tallest stack and prune to this stack.....
  my($tall,$ntall,@rest) = sort{scalar(@{$tmp{$b}}) <=> scalar(@{$tmp{$a}})}keys %tmp;
  if($debug == 1){
    print STDERR "TALLEST STACK:\n";
    foreach my $s (@{$tmp{$tall}}){
      print STDERR " ",$s->getID(),": \(",$s->getMinQueryStart(),",",$s->getMaxQueryEnd(),"\), ",$s->countHSPs()," HSPs\n";
    }

    if(exists $tmp{$ntall} && exists $tmp{$tall}){
      print STDERR "Pruned to tallest stack: ",$self->countSubjects(), " total, ",scalar(@{$tmp{$tall}}), " tallest, ",scalar(@{$tmp{$ntall}})," next\n" if $debug == 1;
    }
  }

  ##??do I need a minimum size for the stack here???
  ##also some difference between tallest and next??...if same then what??
  ##I think that a minimum size should be > 5 or even 10 as by definition a repeat
  ##will have many menbers....
  $self->resetSubjects(); ##removes all subjects....rebuild if tall stack 
  if(exists $tmp{$tall} && scalar(@{$tmp{$tall}}) >= 5){  ##ignore difference between tallest and next for now
    foreach my $s (@{$tmp{$tall}}){
      $self->addSubject($s);
    }
  }
  
}

## returns 1 if verified, 0 if not
sub getVerified {
  my $self = shift;
  if(!exists $self->{'verified'}){
    $self->verifyRepeat(15);  ##use slop of 15 for default if call getVerified
                              ##without first calling verifyRepeat
  }
  return $self->{'verified'};
}

##alternatively, could return the middle and then use $slop when ruling out or in
#3sequences based upon locations of HSPs
sub getRepeatBorders{
  my $self= shift;
  my $max3 = 0;
  my $min5 = 1e10;
##need to deal with > 1 hsp case..
  
  ##for one HSP...perhaps this is all I need for now...will miss those things 
  ##just in the middle but...
  foreach my $s ($self->getSubjects()){
    $min5 = $s->getMinQueryStart() if $s->getMinQueryStart < $min5;
    $max3 = $s->getMaxQueryEnd() if $s->getMaxQueryEnd() > $max3;
  }
  $self->{'minRepeatBorder'} = $min5;
  $self->{'maxRepeatBorder'} = $max3;
}

sub getRepeatStart{
  my $self = shift;
  if(!exists $self->{'minRepeatBorder'}){
    $self->getRepeatBorders();
  }
  return $self->{'minRepeatBorder'};
}

sub getRepeatEnd{
  my $self = shift;
  if(!exists $self->{'maxRepeatBorder'}){
    $self->getRepeatBorders();
  }
  return $self->{'maxRepeatBorder'};

}

sub getRepeatLength{
  my $self = shift;
  return $self->getRepeatEnd() - $self->getRepeatStart();
}

##if $type == 1 trhen just summary...
sub toString{
  my($self,$type) = @_;
  if($type eq ""){
    $type = 0;
  }
  my $string = "REPEAT: ".$self->countSubjects()." subjects, \(".$self->getRepeatStart().",".$self->getRepeatEnd()."\)\n";
  if($type != 1){
    foreach my $s ($self->getSubjects()){
      $string .= "  ".$s->getID(). ": \(".$s->getMinQueryStart().",".$s->getMaxQueryEnd()."\), ".$s->countHSPs()." HSPs\n";
    }
  }
  return $string;
}
1;
