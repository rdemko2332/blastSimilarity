#!/usr/bin/perl

## perl module for dealing with blast results...
## note that the parseBlast method requires the blast be run with -noseqs option
##
## currently will be meant to  only deal with blastn as is meant to do analysis
## of consistency of mRNA vs Genomic sequence of mRNA vs mRNA to determine differential
## splicing things....

## NOTE:  Will always represtent the locations with the start < end...the direction indicates  
## whether one or the other (query/subject) should then be reversed..

## Brian Brunk 5/17/99

package HSP;

use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(new _inititialize getScore getFrame getIdentities getPositives getID getSubjectStart getSubjectEnd getCorrectedQueryStart getCorrectedQueryEnd getQueryStart getQueryEnd getMatchLength getMatchPercent getPValue getDirection getSimilaritySpan toString);

use strict;

sub new {
  my $class = shift;  ##removes the class name from constructor...so will not be passed to initialize
  my $self = {};
  bless $self, $class;
  $self->_inititialize(@_);
  return $self;
}

sub _inititialize{
  my($self,$subjectLength,$ID,$qStart,$qEnd,$sStart,$sEnd,$matchLength,$matchPercent,$pValue,$direction,$frame,$score,$identities,$positives) = @_;
##  $self->{'subjectLength'} = $sbjctLength; ##shouldn't need to store as am revcomping int he init method.

  $self->{"ID"} = $ID;
  $self->{"matchLength"} = $matchLength;
  $self->{"matchPercent"} = $matchPercent;
  $self->{"pValue"} = $pValue;
  $self->{"direction"} = $direction;
  $self->{'frame'} = $frame;
  $self->{'score'} = $score;
  $self->{'identities'} = $identities;
  $self->{'positives'} = $positives;

  ##note that can infer direction always by the query start and end:
  ##if direction is minus then reverse the numbers!!!! 
  $self->{"direction"} = 1;
  if($qEnd < $qStart){
    ##for query just reverse the start asnd end as is already numbered correctly
    $self->{"qStart"} = $qEnd;
    $self->{"qEnd"} = $qStart;
    $self->{"direction"} = 0;

  }else{
    $self->{"qStart"} = $qStart;
    $self->{"qEnd"} = $qEnd;
  }
  if($sEnd < $sStart){
    ##for subject just reverse the start asnd end as is already numbered correctly (for tblastn)
    $self->{"sStart"} = $sEnd;
    $self->{"sEnd"} = $sStart;
    $self->{"direction"} = 0;

  }else{
    $self->{"sStart"} = $sStart;
    $self->{"sEnd"} = $sEnd;
  }
}

sub getScore { my $self = shift; return $self->{'score'}; }
sub getFrame { my $self = shift; return $self->{'frame'}; }
sub getIdentities { my $self = shift; return $self->{'identities'}; }
sub getPositives { my $self = shift; return $self->{'positives'}; }

sub getID{
  my $self = shift;
  return $self->{"ID"};
}

sub getSubjectStart{
  my $self = shift;
  return $self->{"sStart"};
}

sub getSubjectEnd{
  my $self = shift;
  return $self->{"sEnd"};
}

##return orientation corrected start
##have already fixed this in init....change so if miss method call will do correct thing
sub getCorrectedQueryStart{
  my $self = shift;
  return $self->{"qStart"};
}

##return orientation corrected end 
sub getCorrectedQueryEnd{
  my $self = shift;
  return $self->{"qEnd"};
}

sub getQueryStart{
  my $self = shift;
  return $self->{"qStart"};
}

sub getQueryEnd{
  my $self = shift;
  return $self->{"qEnd"};
}

sub getMatchLength{
  my $self = shift;
  return $self->{"matchLength"};
}

sub getMatchPercent{
  my $self = shift;
  return $self->{"matchPercent"};
}

sub getPValue{
  my $self = shift;
  return $self->{"pValue"};
}

sub getDirection{
  my $self = shift;
  return $self->{"direction"};
}

sub getSimilaritySpan {
  my($self,$d) = @_;
  if(!$d) { $d = ", "; }
  my $tmp = $self->getIdentities().$d.$self->getPositives().$d.$self->getMatchLength().$d.$self->getScore().$d.$self->getPValue().$d;
  $tmp .= $self->getSubjectStart().$d.$self->getSubjectEnd().$d.$self->getQueryStart().$d.$self->getQueryEnd().$d;
  $tmp .= ($self->getDirection() ? 0 : 1).$d.$self->getFrame();
  return $tmp;
}

sub toString{
  my($self,$space) = @_;
  if(!defined $space){$space = "";}
  return "$space".$self->getID()." MatLen=".$self->getMatchLength()." MatPer=".$self->getMatchPercent().", Query=\(".$self->getQueryStart().",".$self->getQueryEnd()."\), Sbjct=\(".$self->getSubjectStart().",".$self->getSubjectEnd()."\), ".$self->getDirection()."\n";
}


1;
