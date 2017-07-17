#!/usr/bin/perl -w

# Contact: lh3
# Version: 0.1.0

use strict;
use warnings;
use Getopt::Std;

my $IN;
my $OUT;
open (IN, "< col1");
open (OUT, "> column.in");
my @dat;
while (<IN>) {
  chomp;
  push(@dat,$_);
}
my @l0=split(/:/,$dat[0]);
my @l1=split(/:/,$dat[1]);
my $len=@l0;
for (my $i=0;$i<$len;$i++) {
  if ($l0[$i] ne $l1[$i]) {
    print OUT "$i\n";
    last;
  }
}
