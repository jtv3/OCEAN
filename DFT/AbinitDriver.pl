#!/usr/bin/perl
# Copyright (C) 2015 - 2020 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;
use POSIX qw(ceil);
use Cwd 'abs_path';
use File::Compare;
use File::Spec::Functions;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/AbinitDriver\.pl/;
#  $ENV{"OCEAN_BIN"} = $1;
  $ENV{"OCEAN_BIN"} = abs_path( $1 );
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
if (!$ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }
if (! $ENV{"OCEAN_ABINIT"} ) {$ENV{"OCEAN_ABINIT"} = $ENV{"OCEAN_BIN"} . "/abinit"; }
if (! $ENV{"OCEAN_CUT3D"} ) {$ENV{"OCEAN_CUT3D"} = $ENV{"OCEAN_BIN"} . "/cut3d"; }

####################################
# Executables to be run
# kgen.x
# PP
# ABINIT
# avecs?

my $RunKGen = 0;
my $RunPP = 0;
my $RunABINIT = 0;
my $screenRUN = 0;
my $bseRUN = 0;

my @GeneralFiles = ("para_prefix", "ser_prefix", "calc" );

my @KgenFiles = ("nkpt", "k0.ipt", "qinunitsofbvectors.ipt", "screen.nkpt", "screen.k0", "dft.split");
my @BandFiles = ("nbands", "screen.nbands");
my @AbinitFiles = ( "rscale", "rprim", "ntype", "natoms", "typat",
    "verbatim", "coord", "taulist", "ecut", "etol", "nrun", "wftol", 
    "fband", "occopt", "ngkpt", "abpad", "nspin", "smag", "metal", "degauss", 
    "dft.calc_stress", "dft.calc_force", "tot_charge", "dft");
my @PPFiles = ("pplist", "znucl", "ppdir");
my @OtherFiles = ("epsilon", "screen.mode");

my $AbiVersion = 0;
my $AbiMinorV = -1;
my $AbiSubV = -1;

foreach (@PPFiles) {
  if ( -e $_ ) {
    if ( `diff -q $_ ../Common/$_` ) {
      $RunPP = 1;
      print "$_ differs\n";
      last;
    }
  }
  else {
    $RunPP = 1;
    print "$_ not found\n";
    last;
  }
}
unless( $RunPP )
{
  if( -e '../Common/psp8.pplist' )
  {
    if ( `diff -q psp8.pplist ../Common/psp8.pplist` )
    {
      $RunPP = 1;
    }
  }
}
unless ($RunPP) {
  $RunPP = 1;
  if (open STATUS, "pp.stat" ) {
    if (<STATUS> == 1) { $RunPP = 0; }
  }
  close STATUS;
}
if ($RunPP != 0) {
  `rm -f pp.stat`;
} 

if ( $RunPP ) {
  $RunABINIT = 1;
}
else {
  foreach (@AbinitFiles) {
    if ( -e $_ ) {
      if ( `diff -q $_ ../Common/$_` ) {
        $RunABINIT = 1;
        last;
      }
    }
    else {
      $RunABINIT = 1;
      last;
    }
  }
}
unless ($RunABINIT) {
 $RunABINIT = 1;
  if (open STATUS, "abinit.stat" ) {
    if (<STATUS> == 1) { $RunABINIT = 0; }
  }
  close STATUS;
}
if ($RunABINIT) {
  print "Differences found for density run. Clearing all old data\n";
#  die;
  my @dirlisting = <*>;
  foreach my $file (@dirlisting) {
    chomp($file);
    `rm -r $file`;
  }
  $RunPP = 1;
}
else {
  `touch old`;
  open IN, "density.log" or die "Failed to open density.log\n$!";
  while( my $line = <IN> )
  {
    if( $line =~ m/\.Version\s(\d+)\.(\d+)\.(\d)/ )
    {
      $AbiVersion = $1;
      $AbiMinorV = $2;
      $AbiSubV = $3;
      last;
    }
  }
  close IN;

}

open OUT, ">core" or die;
print OUT "1\n";
close OUT;


#unless ($RunABINIT || $RunPP || $RunKGen ) {
#  print "Nothing needed in ABINIT Stage\n";
#  open GOUT, ">abinitstage.stat" or die;
#  print GOUT "1";
#  close GOUT;
#  exit 0;
#}

open GOUT, ">abinitstage.stat" or die;
print GOUT "0";
close GOUT;

foreach (@GeneralFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}
foreach (@KgenFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}
foreach (@BandFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}

open PARA, "para_prefix" or die "$!";
my $para_prefix = <PARA>;
chomp($para_prefix);
close PARA;

open SER, "ser_prefix" or die "$!";
my $ser_prefix = <SER>;
chomp($ser_prefix);
close SER;

#if  ($RunKGen) {
#  foreach (@KgenFiles) {
#    system("cp ../Common/$_ .") == 0 or die;
#  } 
#`echo "1" > kgen.stat`;
#}

#if ($RunPP) {
  foreach (@PPFiles) {
    system("cp ../Common/$_ .") == 0 or die;
  } 
#`echo "1" > pp.stat`;
#}

open IN, "calc" or die "Failed to open calc\n";
<IN> =~m/(\w+)/ or die "Failed to parse calc\n";
my $calc = $1;
close IN;

my $old_screen_mode;
if( -e "screen.mode" )
{
  open IN, "screen.mode" or die "Failed to open screen.mode\n$!";
  <IN> =~m/(\w+)/ or die "Failed to parse screen.mode\n";
  $old_screen_mode = $1;
  close IN;
}
else
{
  $old_screen_mode = '';
}

#if ($RunABINIT) {
  foreach (@AbinitFiles, @OtherFiles) {
    system("cp ../Common/$_ .") == 0 or die;
  } 
#}


#############################################
my $ecut = `cat ecut`;
chomp($ecut);
open OUT, ">ecutRy" or die;
print OUT "$ecut Ry\n";
close OUT;

my $degauss = `cat degauss`;
chomp( $degauss );
$degauss = $degauss / 2 ;
my $tsmear = "tsmear $degauss\n";


# test screen.nkpt, screen.nbands
open NKPT, "screen.nkpt" or die "Failed to open screen.nkpt\n";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse. screen.nkpt\n";
my @screennkpt = ($1, $2, $3);
close NKPT;
open NKPT, "nkpt" or die "Failed to open nkpt\n";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse. nkpt\n";
my @nkpt = ($1, $2, $3);
close NKPT;
my $screennbands;
my $nbands;
open NBANDS, "screen.nbands" or die "Failed to open screen.nbands\n";
<NBANDS> =~ m/(\d+)/ or die "Failed to parse screen.nbands\n";
$screennbands = $1;
close NBANDS;
open NBANDS, "nbands" or die "Failed to open nbands\n";
<NBANDS> =~ m/(\d+)/ or die "Failed to parse nbands\n";
$nbands = $1;
close NBANDS;
open NSPN, "nspin" or die "Failed to open nspin\n";
<NSPN> =~ m/(\d)/ or die "Failed to parse nspin\n";
my $nspn = $1;
close NSPN;
open IN, "metal" or die "Failed to open metal\n";
my $metal = 1;
if( <IN> =~ m/false/i )
{
  $metal = 0;
}
close IN;

open IN, "dft.calc_stress" or die "Failed to open dft.calc_stress\n$!";
my $calc_stress;
if ( <IN> =~ m/true/i ) 
{
  $calc_stress = 1;
}
else
{
  $calc_stress = 0;
}
close IN;

open IN, "dft.calc_force" or die "Failed to open dft.calc_force\n$!";
my $calc_force;
if ( <IN> =~ m/true/i )
{
  $calc_force = 2;
}
else
{
  $calc_force = 0;
}
close IN;


if ( $nkpt[0] + $nkpt[1] + $nkpt[2] == 0 ) {
  `cp nkpt screen.nkpt`;
  @screennkpt = @nkpt;
  if ( $screennbands == 0 ) {
    `cp nbands screen.nbands`;
    $screennbands = $nbands;
  }
  elsif ( $nbands > $screennbands ) {
    die "screen.nbands must be larger than nbands\b";
  }
}


# test the directory for the SCREENING run first
#my $screenDIR = sprintf("%03u%03u%03u", $screennkpt[0], $screennkpt[1], $screennkpt[2] );
my $screenDIR = "SCREEN";
if ( -d $screenDIR ) {
  chdir $screenDIR;
  if (-e "abinit.stat") {
    if ( `diff -q nbands ../screen.nbands`) {
      open NBANDS, "nbands" or die "Failed to open `pwd`/nbands\n";
      <NBANDS> =~ m/(\d+)/ or die "Failed to parse nbands\n";
      my $tmpnbands = $1;
      close NBANDS;
      $screenRUN = 1 if ( $tmpnbands < $screennbands);
    }
    $screenRUN = 1 if ( `diff -q k0.ipt ../screen.k0` );
    $screenRUN = 1 if ( `diff -q nkpt ../screen.nkpt` );
  }
  else {
    $screenRUN = 1;
  }
  chdir "../"
}
else {
  $screenRUN = 1;
}
open IN, "screen.mode" or die "Failed to open screen.mode";
<IN> =~m/(\w+)/ or die "Failed to parse screen.mode\n";
my $screen_mode = $1;
close IN;

open IN, "screen.mode" or die "Failed to open screen.mode\n";
<IN> =~m/(\w+)/ or die "Failed to parse screen.mode\n";
my $screen_mode = $1;
close IN;
if( $calc =~ m/val/i )
{
  $screenRUN = 0 unless( $screen_mode =~ m/grid/i );
}
if( $screenRUN == 0 && $screen_mode =~ m/grid/i )
{
  unless( $old_screen_mode =~ m/grid/i )
  {
    print "Need screening for valence: $old_screen_mode\n";
    $screenRUN = 1;
  }
}

if ($screenRUN == 1) {
  print "Need to run for SCREENING\n";
#  die;
  `rm -rf $screenDIR`;
  mkdir $screenDIR;
}
else {
  `mkdir -p  $screenDIR`;
  `touch $screenDIR/old`;
}

# test the directory for the NBSE run
my $bseDIR = sprintf("%03u%03u%03u", $nkpt[0], $nkpt[1], $nkpt[2] );
if ( -d $bseDIR) {
  chdir $bseDIR;
  if (-e "abinit.stat") 
  {
    foreach ( "nkpt", "k0.ipt", "dft.split" )
    {
      if( compare( "$_", "../$_") != 0 )
      {
        $bseRUN = 1;
        print "$_ differs\n";
        last;
      }
    }
    if ( `diff -q nbands ../nbands`) 
    {
      open NBANDS, "nbands" or die "Failed to open `pwd`/nbands\n";
      <NBANDS> =~ m/(\d+)/ or die "Failed to parse nbands\n";
      my $tmpnbands = $1;
      close NBANDS;
      $bseRUN = 1 if ( $tmpnbands < $nbands);
    }
  }
  else {
    $bseRUN = 1;
  }
  if( $bseRUN == 0 )
  {
    $bseRUN = 2 if ( `diff -q qinunitsofbvectors.ipt ../qinunitsofbvectors.ipt` );
  }
  chdir "../"
}
else {
  $bseRUN = 1;
}

# To better support RIXS/valence screen and bse wfns calcs are always separate
#if ( $screenRUN == 1 && $nkpt[0] == $screennkpt[0] && $nkpt[1] == $screennkpt[1] && $nkpt[2] == $screennkpt[2] ) {
#  $bseRUN = 0;
#}

if ($bseRUN == 1 ) {
  print "Need run for the BSE\n";
  #die;
  `rm -rf $bseDIR`;
  mkdir $bseDIR;
}
elsif( $bseRUN == 2 ) 
{
  print "Need run occupied states for the BSE\n"; 
}
else
{
  `touch $bseDIR/old`;
}

#if ($RunKGen) {
#  print "Running kgen.x\n";
#  `cp nkpt kmesh.ipt`;
#  system("$ENV{'OCEAN_BIN'}/kgen.x") == 0 or die "KGEN.X Failed\n";
#  `echo "1" > kgen.stat`;
#}

if ($RunPP) {
  if( -e '../Common/psp8.pplist' )
  {
    `cp ../Common/psp8.pplist .`;
    open IN, "psp8.pplist" or die;
    open OUT, ">", "finalpplist";
    while( my $line = <IN> )
    {
      chomp $line;
      $line = catdir( updir, "Common", "psp", $line );
      $line = abs_path( $line );
      print OUT $line . "\n";
    }
    close OUT;
    close IN;
  }
  else
  {
    system("$ENV{'OCEAN_BIN'}/pp.pl znucl pplist finalpplist") == 0
      or die "Failed to run pp.pl\n";
  }
  `echo "1" > pp.stat`;
}

#############################################
# Run abinit
# Determin what type of run is being done, par or seq

my $AbinitType = "seq";

if ($RunABINIT) {
  `echo symmorphi 0 > abfile`;
  `echo autoparal 1 >> abfile`;
  `echo chksymbreak 0 >> abfile`;
  `echo 'acell ' >> abfile`;
  `cat rscale >> abfile`;
  `echo rprim >> abfile`;
  `cat rprim >> abfile`;
  `echo 'ntypat ' >> abfile`;
  `cat ntype >> abfile`;
  `echo 'znucl ' >> abfile`;
  `cat znucl >> abfile`;
  `echo 'natom ' >> abfile`;
  `cat natoms >> abfile`;
  `echo 'typat ' >> abfile`;
  `cat typat >> abfile`;
  `cat coord >> abfile`;
  `cat taulist >> abfile`;
  `echo 'ecut ' >> abfile`;
  `cat ecutRy >> abfile`;
  `echo 'nstep ' >> abfile`;
  `cat nrun >> abfile`;
  `echo 'diemac ' >> abfile`;
  `cat epsilon >> abfile`;
  `cat verbatim >> abfile`;
  `echo 'occopt ' >> abfile`;
  `cat occopt >> abfile`;
  `echo "$tsmear" >> abfile`;
  `echo 'npfft 1' >> abfile`;
  `echo  'charge ' >> abfile`;
  `cat tot_charge >> abfile`;
  `echo  'nsppol ' >> abfile`;
  `cat nspin >> abfile`;
  if( $nspn == 2 )
  {
    `echo 'spinat' >> abfile`;
    `cat smag >> abfile`;
  }

#if ($AbinitType eq "par" ) {
#  die "not an option\n";
#}
#########################
# seq Abinit run
#else {

### Clean ####
  `rm -f density.out`;
  `rm -f SCx_DEN`;
  `rm -f SCx_EIG`;
  `rm -f SCx_WFK`;
#  `rm -f Run.????.out`;
### Done cleaning ###


  open FILES, ">denout.files";
  print FILES "inai.denout\n"
            . "density.out\n"
            . "SC\n"
            . "SCx\n"
            . "Scxx\n";
  close FILES;
  `cat finalpplist >> denout.files`;
  
  `cat abfile > inai.denout`;
  `echo 'fband ' >> inai.denout`;
  `cat fband >> inai.denout`;
  `echo prtden 1 >> inai.denout`;
  `echo prtpot 1 >> inai.denout`;
  `echo kptopt 1 >> inai.denout`;
  `echo 'ngkpt ' >> inai.denout`;
  `cat ngkpt >> inai.denout`;
  `echo 'toldfe ' >> inai.denout`;
  `cat etol >> inai.denout`;
  `echo "optstress $calc_stress" >> inai.denout`; 
  `echo "optforces $calc_force" >> inai.denout`;
#  `echo prtdos 3 >> inai.denout`;
#  `echo prtdosm 1 >> inai.denout`;

  $para_prefix =~ m/(\d+)/;
  my $max_cpus = $1;
  $max_cpus = 1 if( $max_cpus < 1 );

  my $new_ncpu = 1;
  my $test_prefix;

  if( $max_cpus > 1 ) 
  {
    print "Testing SCF Run parallelism\n";
    if( $ser_prefix =~ m/!/ )
    {
      $test_prefix = $para_prefix;
      $test_prefix =~ s/\d+/1/;
      print "    Serial prefix: $test_prefix\n";
    }
    else
    {
      $test_prefix = $ser_prefix;
    }
    `cp inai.denout par_test.in`;
    `echo paral_kgb 1 >> par_test.in`;

    `echo max_ncpus $max_cpus >> par_test.in`;
    open FILES, ">par_test.files";
    print FILES "par_test.in\n"
              . "par_test.out\n"
              . "SC\n"
              . "SCx\n"
              . "Scxx\n";
    close FILES;
    `cat finalpplist >> par_test.files`;
    system("$test_prefix $ENV{'OCEAN_ABINIT'} < par_test.files > pt.log 2> pt.err");

    open IN, "pt.log" or die "Failed to open pt.log\n$!";
    while( my $line = <IN> )
    {
      if( $line =~ m/\.Version\s(\d+)\.(\d+)\.(\d)/ )
      {
        $AbiVersion = $1;
        $AbiMinorV = $2;
        $AbiSubV = $3;
        last;
      }
    }
    close IN;


    if( ( $AbiVersion < 7 ) || ( $AbiVersion == 7 && $AbiMinorV < 6 ) )
    {
      print "Abinit is pre-7.6: $AbiVersion.$AbiMinorV\n";
      `rm -rf par_test.out`;
      `sed '/autoparal/d' inai.denout > tmp`;
      `mv tmp inai.denout`;
      `echo istwfk *1 >> inai.denout`;
      `cp inai.denout par_test.in`;
      `echo paral_kgb -$max_cpus >> par_test.in`;
      `echo paral_kgb 1 >> inai.denout`;

      system("$test_prefix $ENV{'OCEAN_ABINIT'} < par_test.files > pt.log 2> pt.err");

      open IN, "pt.log" or die "Failed to open pt.log\n$!";
      while( my $line = <IN> )
      {
        last if $line =~ m/nproc     npkpt/;
      }
      my $npfft = 0;
      while( $npfft != 1 )
      {
        <IN> =~ m/^\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse pt.log\n";
        $new_ncpu = $1;
        $npfft = $5;
      }

    }
    else
    {
      open IN, "pt.log" or die "Failed to open pt.log\n$!";
      while( my $line = <IN> )
      {
        last if $line =~ m/npimage\|/;
      }
    
      <IN>;
      <IN> =~ m/^\s+(\d+)\|\s+(\d+)\|\s+(\d+)\|\s+(\d+)\|\s+(\d+)/ or die "Failed to parse pt.log\n";
      $new_ncpu = $1*$2*$3*$4*$5;
      close IN;
    }
  }  # max_cpus > 1

  $test_prefix = $para_prefix;
  $test_prefix =~ s/\d+/$new_ncpu/;
  print "Self-Consistent Density Run\n";
  print "$test_prefix $ENV{'OCEAN_ABINIT'} < denout.files > density.log 2> density.err\n";
  system("$test_prefix $ENV{'OCEAN_ABINIT'} < denout.files > density.log 2> density.err") == 0
    or die "Failed to run initial density stage\n$test_prefix $ENV{'OCEAN_ABINIT'}\n";
  `echo 1 > den.stat`;

  `ln -s SCx_DEN SCx_DS0_DEN`;

  # If we didn't do parallelism tests, then we haven't checked the version number
  if( $max_cpus == 1 )
  {
    open IN, "density.log" or die "Failed to open density.log\n";
    while( my $line = <IN> )
    {
      if( $line =~ m/\.Version\s(\d+)\.(\d+)\.(\d)/ )
      {
        $AbiVersion = $1;
        $AbiMinorV = $2;
        $AbiSubV = $3;
        last;
      }
    }
    close IN;
  }


  print "Abi Version $AbiVersion.$AbiMinorV.$AbiSubV\n";

  open CUTIN, ">cut3d.in" or die "Failed to open cut3d.in for writing.\n$!\n";
  if( $AbiVersion <= 7 )
  {
    if( $nspn == 2 ) 
    {
      print CUTIN "SCx_DEN\n1\n0\n6\nrhoofr\n0\n";
    }
    else
    {
      print CUTIN "SCx_DEN\n1\n6\nrhoofr\n0\n";
    }
  }
  else
  {
    if( $nspn == 2 )
    {
      print CUTIN "SCx_DEN\n0\n6\nrhoofr\n0\n";
    }
    else
    {
      print CUTIN "SCx_DEN\n6\nrhoofr\n0\n";
    }
  }
  close CUTIN;

  open CUTIN, ">cut3d2.in" or die "Failed to open cut3d2.in for writing.\n$!\n";
  if( $AbiVersion <= 7 )
  {
    if( $nspn == 2 )
    {
      print CUTIN "SCx_POT\n1\n0\n6\npotofr\n0\n";
    }
    else
    {
      print CUTIN "SCx_POT\n1\n6\npotofr\n0\n";
    }
  }
  else
  {
    if( $nspn == 2 )
    {
      print CUTIN "SCx_POT\n0\n6\npotofr\n0\n";
    }
    else
    {
      print CUTIN "SCx_POT\n6\npotofr\n0\n";
    }
  }
  close CUTIN;


  my $useThis_prefix = $ser_prefix;
  if( $ser_prefix =~ m/!/ ) 
  {
    $useThis_prefix = $para_prefix;
    $useThis_prefix =~ s/\d+/1/;
    print "    Serial prefix: $useThis_prefix\n";
  }


  print "$useThis_prefix $ENV{'OCEAN_CUT3D'} < cut3d.in > cut3d.log 2> cut3d.err\n";
  system("$useThis_prefix $ENV{'OCEAN_CUT3D'} < cut3d.in > cut3d.log 2> cut3d.err") == 0
      or die "Failed to run cut3d\n";

  print "$useThis_prefix $ENV{'OCEAN_CUT3D'} < cut3d2.in > cut3d2.log 2> cut3d2.err\n";
  system("$useThis_prefix $ENV{'OCEAN_CUT3D'} < cut3d2.in > cut3d2.log 2> cut3d2.err") == 0
      or die "Failed to run cut3d\n";
  
  `echo "1" > abinit.stat`;

}

open LOG, "density.log";
my $vb;
while (<LOG>) {
  if ($_ =~ m/nband\s+(\d+)/) {
    $vb = $1;
    last;
  } 
}   
close LOG;


if ( $screenRUN ) {
  print "SCREENING run\n";
  chdir $screenDIR;   
  `cp ../abfile .`;
 # copy all files over
  foreach ( @GeneralFiles, @AbinitFiles, @PPFiles, @OtherFiles) {
    system("cp ../$_ .") == 0 or die "Failed to copy $_\n";
  }
  foreach ( "screen.nkpt", "screen.nbands", "screen.k0", "qinunitsofbvectors.ipt", "finalpplist", "core" ) {
    system("cp ../$_ .") == 0 or die "Failed to copy $_\n";
  }
  `cp screen.k0 k0.ipt`;
  `cp screen.nkpt nkpt`;
  `cp screen.nbands nbands`;
 # run KGEN
  print "Running kgen2.x\n";
  `cp screen.nkpt kmesh.ipt`;
  `echo 0.0 0.0 0.0 > qinunitsofbvectors.ipt`;
  system("$ENV{'OCEAN_BIN'}/kgen2.x") == 0 or die "KGEN.X Failed\n";
  `echo "1" > kgen.stat`;

 
  open ABPAD, "abpad" or die;
  my $abpad = <ABPAD>;
  close ABPAD;
  $screennbands += $abpad; 
  
  `echo "nband $screennbands" >> abfile`;
  `echo "nbdbuf $abpad" >> abfile`;
  `echo 'iscf -2' >> abfile`;
  `echo 'tolwfr ' >> abfile`;
  `cat wftol >> abfile`;
  `echo getden -1 >> abfile`;
  `echo kptopt 0 >> abfile`;
  `echo "istwfk *1" >> abfile`;
  
  my $nfiles = 1;
  for (my $runcount = 1; $runcount <= $nfiles; $runcount++ ) 
  {
    my $abfilename = sprintf("inabinit.%04i", $runcount );
    `cp abfile "$abfilename"`;
    my $kptfile =  sprintf("kpts.%04i", $runcount);
    `cat $kptfile >> $abfilename`;

    my $deninFiles = sprintf("denin.files.%04i", $runcount );
    open FILES, ">$deninFiles";
    my $Runout = sprintf('Run.%04i.out', $runcount );
    my $RUN = sprintf('RUN%04i', $runcount);
    print FILES "$abfilename\n"
              . "$Runout\n"
              . "../SCx\n"
              . "$RUN\n"
              . "SCxx" . $runcount . "\n";
    close FILES;
    `cat finalpplist >> $deninFiles`;

    $para_prefix =~ m/(\d+)/;
    my $max_cpus = $1;
    $max_cpus = 1 if( $max_cpus < 1 );
    open FILES, ">par_test.files";
    print FILES "par_test.in\n"
              . "par_test.out\n"
              . "../SC\n"
              . "SCx\n"
              . "Scxx\n";
    close FILES;
    `cat finalpplist >> par_test.files`;

    my $test_prefix;
    if( $ser_prefix =~ m/!/ )
    {
      $test_prefix = $para_prefix;
      $test_prefix =~ s/\d+/1/;
#      print "    Serial prefix: $test_prefix\n";
    }
    else
    {
      $test_prefix = $ser_prefix;
    }

    my $new_ncpu;
    if( $max_cpus > 1 ) 
    {
      if( ( $AbiVersion < 7 ) || ( $AbiVersion == 7 && $AbiMinorV < 6 ) )
      {
        print "Abinit is pre-7.6: $AbiVersion.$AbiMinorV\n";
        `sed '/autoparal/d' $abfilename > tmp`;
        `mv tmp $abfilename`;
        `sed /iscf/d  $abfilename > par_test.in`;

        `echo paral_kgb 1 >> $abfilename`;
        `echo paral_kgb -$max_cpus >> par_test.in`;

        system("$test_prefix $ENV{'OCEAN_ABINIT'} < par_test.files > pt.log 2> pt.err");

        open IN, "pt.log" or die "Failed to open pt.log\n$!";
        while( my $line = <IN> )
        {
          last if $line =~ m/nproc     npkpt/;
        }
        my $npkpt = 0;
        while( <IN> =~ m/^\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ )
        {
          $npkpt = $2 if( $2 > $npkpt );
        }
        $new_ncpu = $npkpt;

        close IN;
      }
      else
      {
        `cp $abfilename par_test.in`;
        `echo paral_kgb 1 >> par_test.in`;
        `echo max_ncpus $max_cpus >> par_test.in`;

        system("$test_prefix $ENV{'OCEAN_ABINIT'} < par_test.files > pt.log 2> pt.err");

        open IN, "pt.log" or die "Failed to open pt.log\n$!";
        while( my $line = <IN> )
        {
          last if $line =~ m/npimage\|/;
        }
        <IN>;
        <IN> =~ m/^\s+(\d+)\|\s+(\d+)\|\s+(\d+)\|\s+(\d+)\|\s+(\d+)/ or die "Failed to parse pt.log\n";
        $new_ncpu = $1*$2*$3*$4*$5;
        close IN;
      }
    }
    else
    {
      $new_ncpu = 1;
    }
    my $test_prefix = $para_prefix;
    $test_prefix =~ s/\d+/$new_ncpu/;

    my $denin = sprintf("denin.files.%04i", $runcount);
    print "$test_prefix $ENV{'OCEAN_ABINIT'} < $denin > ABINIT.$runcount.log";
    system("$test_prefix $ENV{'OCEAN_ABINIT'} < $denin > ABINIT.$runcount.log 2> ABINIT.$runcount.err") == 0 or
      die "$!\n";
    print "\n";
  }


  my $natoms = `cat natoms`;
  my $fband = `cat fband`;
  $screennbands = `cat screen.nbands`;
  my $true_vb = $vb - ceil( $natoms*$fband );
  print "$vb\t$true_vb\n";
  my $cb = sprintf("%.0f", $vb - 2*$natoms*$fband);
  $cb = 1 if ($cb < 1);
  open BRANGE, ">brange.ipt" or die;
  if( $metal == 1 )
  {
    print BRANGE "1  $vb\n"
               . "$cb $screennbands\n";
  }
  else
  {
    print "1 $true_vb\n";
    print BRANGE "1 $true_vb\n";
    $true_vb++;
    print BRANGE "$true_vb  $screennbands\n";
    print "$true_vb  $screennbands\n";
  }
  close BRANGE;

  open STATUS, ">abinit.stat" or die;
  print STATUS "1";
  close STATUS;

####
  chdir "../";
}

if ( $bseRUN ) {
  chdir $bseDIR;
  print "BSE run\n";
  `cp ../abfile .`;
 # copy all files over
  foreach ( @GeneralFiles, @AbinitFiles, @PPFiles, @OtherFiles) {
    system("cp ../$_ .") == 0 or die "Failed to copy $_\n";
  }
  foreach ( "nkpt", "nbands", "k0.ipt", "qinunitsofbvectors.ipt", "finalpplist", "core", "dft.split" ) {
    system("cp ../$_ .") == 0 or die "Failed to copy $_\n";
  }
 # run KGEN
  print "Running kgen2.x\n";
  `cp nkpt kmesh.ipt`;
  system("$ENV{'OCEAN_BIN'}/kgen2.x") == 0 or die "KGEN.X Failed\n";
  `echo "1" > kgen.stat`;

 
  open ABPAD, "abpad" or die;
  my $abpad = <ABPAD>;
  close ABPAD;
  my $temp_band = $nbands + $abpad; 
  my $valBandString = "fband " .  `cat fband`;
  
#  `echo "nband $temp_band" >> abfile`;
#  `echo "nbdbuf $abpad" >> abfile`;
  `echo 'iscf -2' >> abfile`;
  `echo 'tolwfr ' >> abfile`;
  `cat wftol >> abfile`;
  `echo getden -1 >> abfile`;
  `echo kptopt 0 >> abfile`;
  `echo "istwfk *1" >> abfile`;

#  open NRUNS, "Nfiles" or die;
#  my $nfiles = <NRUNS>;
#  close NRUNS;
  my $nfiles = 1;
  if( open IN, "dft.split" )
  {
    if( <IN> =~ m/t/i )
    {
      open IN, "qinunitsofbvectors.ipt" or die "Failed to open qinunitofbvectors\n$!";
#      <IN> =~ m/([+-]?\d+\.?\d*([eE][+-]?\d+)?)\s+([+-]?\d+\.?\d*([eE][+-]?\d+)?)\s+([+-]?\d+\.?\d*([eE][+-]?\d+)?)/
      <IN> =~ m/([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+/
                  or die "Failed to parse qinunitsofbvectors.ipt\n";
      print "$1\t$3\t$5\n";
      my $fake_qmag = abs($1) + abs($3) + abs($5);
      close IN;
      if( $fake_qmag > 0.000000001 )
      {
        $nfiles = 2;
        print "DFT will be split up\n";
        if( $bseRUN == 2 )
        {
          $nfiles = 1;
          print "  conduction bands re-used\n";
        }
      }
      else
      {
        print "DFT split requested, but q=0! (This is not a problem)\n";
      }
    }
    close IN;
  }



  for (my $runcount = 1; $runcount <= $nfiles; $runcount++ )
  {
    my $abfilename = sprintf("inabinit.%04i", $runcount );
    `cp abfile "$abfilename"`;

    if( $runcount < $nfiles ) 
    {
      chomp $valBandString;
      `echo $valBandString >> $abfilename`;
    }
    else
    {
      `echo "nband $temp_band" >> $abfilename`;
      `echo "nbdbuf $abpad" >> $abfilename`;
    }

    my $kptfile =  sprintf("kpts.%04i", $runcount);
    `cat $kptfile >> $abfilename`;
  
    my $deninFiles = sprintf("denin.files.%04i", $runcount );
    open FILES, ">$deninFiles";
    my $Runout = sprintf('Run.%04i.out', $runcount );
    my $RUN = sprintf('RUN%04i', $runcount);
    print FILES "$abfilename\n"
              . "$Runout\n"
              . "../SCx\n"
              . "$RUN\n"
              . "SCxx" . $runcount . "\n";
    close FILES;
    `cat finalpplist >> $deninFiles`;

    $para_prefix =~ m/(\d+)/;
    my $max_cpus = $1;
    $max_cpus = 1 if( $max_cpus < 1 );
    open FILES, ">par_test.files";
    print FILES "par_test.in\n"
              . "par_test.out\n"
              . "../SC\n"
              . "SCx\n"
              . "Scxx\n";
    close FILES;
    `cat finalpplist >> par_test.files`;

    my $test_prefix;
    if( $ser_prefix =~ m/!/ )
    {
      $test_prefix = $para_prefix;
      $test_prefix =~ s/\d+/1/;
#      print "    Serial prefix: $test_prefix\n";
    }
    else
    {
      $test_prefix = $ser_prefix;
    }

    my $new_ncpu;
    if( $max_cpus > 1 ) 
    {
      if( ( $AbiVersion < 7 ) || ( $AbiVersion == 7 && $AbiMinorV < 6 ) )
      {
        print "Abinit is pre-7.6: $AbiVersion.$AbiMinorV \n";
        `sed '/autoparal/d' $abfilename > tmp`;
        `mv tmp $abfilename`;
    #    `cp $abfilename par_test.in`;
        `sed /iscf/d  $abfilename > par_test.in`;

        `echo paral_kgb 1 >> $abfilename`;
        `echo paral_kgb -$max_cpus >> par_test.in`;

        system("$test_prefix $ENV{'OCEAN_ABINIT'} < par_test.files > pt.log 2> pt.err");

        open IN, "pt.log" or die "Failed to open pt.log\n$!";
        while( my $line = <IN> )
        {
          last if $line =~ m/nproc     npkpt/;
        }
        my $npkpt = 0;
        while( <IN> =~ m/^\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ )
        {
          $npkpt = $2 if( $2 > $npkpt );
        }
        $new_ncpu = $npkpt;

        close IN;
      }
      else
      {
        `cp $abfilename par_test.in`;
        `echo paral_kgb 1 >> par_test.in`;
        `echo max_ncpus $max_cpus >> par_test.in`;

        system("$test_prefix $ENV{'OCEAN_ABINIT'} < par_test.files > pt.log 2> pt.err");

        open IN, "pt.log" or die "Failed to open pt.log\n$!";
        while( my $line = <IN> )
        {
          last if $line =~ m/npimage\|/;
        }
        <IN>;
        <IN> =~ m/^\s+(\d+)\|\s+(\d+)\|\s+(\d+)\|\s+(\d+)\|\s+(\d+)/ or die "Failed to parse pt.log\n";
        $new_ncpu = $1*$2*$3*$4*$5;
        close IN;
      }
    }
    else
    {
      $new_ncpu = 1;
    }
    my $test_prefix = $para_prefix;
    $test_prefix =~ s/\d+/$new_ncpu/;



    my $denin = sprintf("denin.files.%04i", $runcount);
    print "$para_prefix $ENV{'OCEAN_ABINIT'} < $denin > ABINIT.$runcount.log";
    system("$para_prefix $ENV{'OCEAN_ABINIT'} < $denin > ABINIT.$runcount.log 2> ABINIT.$runcount.err") == 0 or
      die "$!\n";
    print "\n";
  }



  my $natoms = `cat natoms`;
  my $fband = `cat fband`;
  my $true_vb = $vb - ceil( $natoms*$fband );
  print "$vb\t$true_vb\n";
  my $cb = sprintf("%.0f", $vb - 2*$natoms*$fband);
  $cb = 1 if ($cb < 1);
  open BRANGE, ">brange.ipt" or die;
  if( $metal == 1 )
  {
    print BRANGE "1  $vb\n"
               . "$cb $nbands\n";
  }
  else
  {
    print "1 $true_vb\n";
    print BRANGE "1 $true_vb\n";
    $true_vb++;
    print BRANGE "$true_vb  $nbands\n";
    print "$true_vb  $nbands\n";
  }
  close BRANGE;

  
  open STATUS, ">abinit.stat" or die;
  print STATUS "1";
  close STATUS;
####
  chdir "../";
}


my $fermi = 'no';

open IN, "density.out" or die "Failed to open density.out\n$!";
while( my $line = <IN> )
{
  if( $line  =~  m/Fermi \(or HOMO\) energy \(hartree\) =\s+([+-]?\d+\.?\d+)/ )
    {
      $fermi = $1 * 2;
      my $eVfermi = $fermi * 13.60569253;
      print "Fermi level found at $eVfermi eV\n";
    }
#    if( $line =~ m/number of electrons\s+=\s+(\d+)/ )
#    {
#      $nelectron = $1;
#    }
  }
close IN;
die "Fermi level not found in scf.out\n" if( $fermi eq 'no' ) ;
#  die "Number of electrons not found in scf.out\n" if( $nelectron eq 'no' );

open FERMI, ">efermiinrydberg.ipt" or die "Failed to open efermiinrydberg\n$!";
print FERMI "$fermi\n";
close FERMI;



print "Abinit stage complete\n";


exit 0;

