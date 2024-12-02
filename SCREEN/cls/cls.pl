use strict;
use File::Spec::Functions;
use File::Copy;

print "Initial script for CLS-EXX testing\n";
print "Run this script from [RUNDIR]/CLS and then run the executable corex.x\n";


## Write cls.inp control file based off of SCREEN/hfinlist ##
my @edgeList;
my %ZNL;
my %Z;
my $file = catfile( updir(), "SCREEN", "hfinlist" );
open IN, "<", $file or die "Failed to open $file\n$!";

while( my $line = <IN> ) {
  $line =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/ or die "Failed to parse $file: $line";
  my $str = sprintf "%s %i %i %i %i", $1, $6, $2, $3, $4;
  push @edgeList, $str;
  $str = sprintf "z%3.3in%2.2il%2.2i", $2, $3, $4;
  $ZNL{ $str } = 1;
  $Z{ $2 } = 1;
}
close IN;

open OUT, ">", "cls.inp" or die;
print OUT scalar @edgeList;
print OUT "\n";
foreach (@edgeList) {
  print OUT $_ . "\n";
}
close OUT;
########

## Copy prjfile and gk files from OPF ##
my $zdiag = catfile( updir(), "OPF", "zpawinfo", "" );
#print $zdiag . "\n";
foreach( keys %ZNL ) {
  my $ZNL = $_;
  my @zfiles = glob "$zdiag" . "gk*" . $ZNL; 
  foreach my $file (@zfiles) {
    copy "$file", ".";
  }
}

foreach( keys %Z ) {
  my $file = catfile( updir(), "OPF", "zpawinfo", "prjfilez" ) . sprintf "%3.3i", $_;
  #print $file . "\n";
  copy $file, ".";
}
  
my @files = ( "Pquadrature", "sphpts", "kmesh.ipt", "avecsinbohr.ipt" );
foreach (@files ) {
  copy catfile( updir(), "CNBSE", $_ ), ".";
}

foreach (@edgeList) {
  $_ =~ m/^(\S+)\s+(\d+)/ or die;
  my $file = catfile( updir(), "PREP", "BSE", sprintf( "parcksv.%2s%4.4i", $1, $2) );
  copy $file, ".";
}


