#!/usr/bin/perl

## module for storing subject data and analyzing potential chimeras

## will contain initially subjects that cluster around a potential
## breakpoint...analyses will be run on these to determine if chimera 
## and then the information retruned if is chimera...

## each Chimera object represents 1 breakpoint.  The determination that this
## represents a breakpoint is the responsibility of the BlastAnal object

## Brian Brunk 11/6/99

package Chimera;

use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(new addLeftSubject addRightSubject getLeftSubjects getRightSubjects sortSubjects getBreakpoint getDepth countLeft countRight getDangling getLeftDangling getRightDangling toString);

use strict;

sub new {
  my($class) = @_;  
  my $self = {};
  bless $self, $class;
  return $self;
}

##convention will be that the "left" side will abut to breakpoint with 3' end
##and "right" side with 5' end thus 5'left3'|5'right3'

sub addLeftSubject{
  my($self,$sbj) = @_;
  push(@{$self->{"left"}},$sbj);
}

sub addRightSubject{
  my($self,$sbj) = @_;
  push(@{$self->{"right"}},$sbj);
}

sub getLeftSubjects{
  my($self,$sbj) = @_;
  return @{$self->{"left"}};
}

sub getRightSubjects{
  my($self,$sbj) = @_;
  return @{$self->{"right"}};
}

##are already sorted
sub sortSubjects{
  my $self = shift;
  @{$self->{"sLeft"}} = sort{$a->getMaxQueryEnd() <=> $b->getMaxQueryEnd()} $self->getLeftSubjects();
  @{$self->{"sRight"}} = sort{$b->getMinQueryStart() <=> $a->getMinQueryStart()} $self->getRightSubjects();
}

##breakpoint is halfway between min leftqueryend and max rightquerystart
sub getBreakpoint{
  my $self = shift;
  if(!exists $self->{"breakpoint"}){
##    if(!exists $self->{"sLeft"}){
##      $self->sortSubjects();  ##   are already sorted.. 
##}
    $self->{"breakpoint"} = $self->{"left"}->[0]->getMaxQueryEnd() + int(($self->{"right"}->[0]->getMinQueryStart() - $self->{"left"}->[0]->getMaxQueryEnd()) / 2);
  }
##  print STDERR "getBreakpoint: $self->{'breakpoint'}\n";
  return $self->{"breakpoint"};
}

##rules here??
##returns 
##  0 - <4 total
##  1 if >=2 on each side
##  1 if 1 on one side and >= 4 on other
sub getDepth{
  my $self = shift;
  if($self->countLeft() >= 2 && $self->countRight() >= 2){
    return 1;
  }elsif(($self->countLeft() + $self->countRight()) > 4){ ##requires at least 1 + 4..
    return 1;
  }
  return 0;
}

sub countLeft{
  my $self = shift;
  return scalar(@{$self->{"left"}});
}

sub countRight{
  my $self = shift;
  return scalar(@{$self->{"right"}});
}

#return 0,1,2 for number of sides dangling...
sub getDangling{
  my($self,$maxEnd) = @_;
  if($maxEnd eq ""){return "ERROR usage: getDangling(int amountDangling)\n";}
  my $ret = 0;
  $ret++ if $self->getLeftDangling($maxEnd) > 0;
  $ret++ if $self->getRightDangling($maxEnd) > 0;
  return $ret;
}

sub getLeftDangling{
  my($self,$maxEnd) = @_;
  if(!exists $self->{"lDangling"}){
    $self->{"lDangling"} = 0;
    foreach my $s ($self->getLeftSubjects()){
      $self->{"lDangling"}++ if $s->getSubject3endMis() > $maxEnd;
    }
  }
  return $self->{"lDangling"};
}

sub getRightDangling{
  my($self,$maxEnd) = @_;
  if(!exists $self->{"rDangling"}){
    $self->{"rDangling"} = 0;
    foreach my $s ($self->getRightSubjects()){
      $self->{"rDangling"}++ if $s->getSubject5endMis() > $maxEnd;
    }
  }
  return $self->{"rDangling"};
}

sub toString{
  my $self = shift;
  my $s = " Left:\n";
  foreach my $l ($self->getLeftSubjects()){
    $s .= "  ".$l->getID().", ".$l->getMinQueryStart()."-".$l->getMaxQueryEnd()." \(".$l->getSubject5endMis().",".$l->getSubject3endMis()."\)\n";
  }
  $s .= " Right:\n";
  foreach my $l ($self->getRightSubjects()){
    $s .= "  ".$l->getID().", ".$l->getMinQueryStart()."-".$l->getMaxQueryEnd()." \(".$l->getSubject5endMis().",".$l->getSubject3endMis()."\)\n";
##    $s .= "  \(" . join(', ', $l->getSubjectStats()) . "\)\n";
  }
  return $s;
}

1;
