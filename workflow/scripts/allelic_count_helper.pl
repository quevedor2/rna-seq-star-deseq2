#!/usr/bin/perl
# perl allelic_count_helper.pl <analysis> <filename> <wild> <alt_col>
#
# perl allelic_count_helper.pl setlines <filename> <min_n>
#    Will report the lines in an input CSV file that have a minimum
#    number of samples that have coverage.
#    Pseudocode, Print Line 1: "0,0,1,1,0" if min_n == 2
#
# perl allelic_count_helper.pl getlines <filename> <filename2>
#    Will take the output from the `setlines` analysis as <filename>,
#    and will extract those lines from the target file as <filename2>
#
# perl allelic_count_helper.pl categorize <filename> <ref_col> <alt_col> <min_cov>
#    Will go through the GATK CollectAllelicCount outpt per sample
#    and classify the Ref/Alt counts as either REFHOM (1), HET (2),
#    or ALTHOM(3).
use List::Util qw/sum/;

#### subroutines ####
sub GetCovLines {
  my($dat, $line, $min_n) = @_;
  chomp $dat;
  my @spl=split(/,/, $dat);

  # Count # of sample with coverage (>0)
  my @nocov;
  foreach(@spl){
    push(@nocov, $_ != 0);
  }
  print($line, "\n") if(sum(@nocov) >= $min_n);
}

sub CategorizeAD {
  my($dat, $ref_col, $alt_col, $min_cov) = @_;

  my @spl=split(/\s/, $dat);
  my $ref = $spl[$ref_col];
  my $alt = $spl[$alt_col];

  ## Categorize the depth per SNP
  my $sums = $ref + $alt;
  if($sums <= $min_cov){
    print "0\n"; ## NO COVERAGE
  } else {
    my $frac = $ref / $sums;
    if ($frac == 1) {
      print "1\n"; ## REF.HOMOZYGOUS
    } elsif ($frac == 0) {
      print "3\n"; ## ALT.HOMOZYGOUS
    } else {
      print "2\n"; ## HETEROZYGOUS
    }
  }
}


#### main ####
my ($analysis, $filename, $wild, $alt_col, $min_cov) = @ARGV;
open(FH, "<", $filename) or die $!;

if($analysis eq 'getlines'){
  my %linenumbers = ();
  while (<FH>) {
    chomp;
    $linenumbers{$_} = 1;
  }
  open(TH, "<", $wild) or die $!;
  $. = 0;
  while (<TH>) {
    print if defined $linenumbers{$.};
  }
} else {
  my $line=1;
  while(<FH>){
    if($analysis eq 'setlines'){
      my $x = GetCovLines($_, $line, $wild);
    } elsif($analysis eq 'categorize'){
      my $x = CategorizeAD($_, $wild, $alt_col, $min_cov);
    }
    $line++;
  }
}
exit;
