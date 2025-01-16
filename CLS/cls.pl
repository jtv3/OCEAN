#!/usr/bin/perl
# Copyright (C) 2025 OCEAN collaboration
#   
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#   
    
use strict;
use File::Copy;
use File::Spec::Functions;
use File::Compare;
use Cwd 'abs_path';
    
use POSIX;
  
require JSON::PP;
use JSON::PP;
use Storable qw(dclone);
use Scalar::Util qw( looks_like_number );

use Data::Dumper;
use Time::HiRes qw( gettimeofday tv_interval );

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/cls\.pl/;
  $ENV{"OCEAN_BIN"} = abs_path( $1 );
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################

my @spdf = ( 's', 'p', 'd', 'f' );

my $json = JSON::PP->new;
$json->canonical([1]);
$json->pretty([1]);

###
my $dataFile = catfile( updir(), "Common", "postDefaultsOceanDatafile" );
my $commonOceanData;
if( open( my $in, "<", $dataFile ))
{
  local $/ = undef;
  $commonOceanData = $json->decode(<$in>);
  close($in);
}
else
{
  die "Failed to open config file $dataFile\n$!";
}

# Early exit if not core-level calculation
my $earlyExit = 0;
$earlyExit = 1 unless( $commonOceanData->{'cls'}->{'enable'} );
$earlyExit = 1 if( $commonOceanData->{'calc'}->{'mode'} eq 'val' );
if( $earlyExit != 0 )
{
  print "EXIT EARLY\n";
#  exit 0;
}


###

$dataFile = catfile( updir(), "DFT", "dft.json" );
my $dftData;
if( open( my $in, "<", $dataFile ))
{
  local $/ = undef;
  $dftData = $json->decode(<$in>);
  close($in);
}
else
{
  die "Failed to open config file $dataFile\n$!";
}

$dataFile = catfile( updir(), "OPF", "opf.json" );
my $opfData;
if( open( my $in, "<", $dataFile ))
{
  local $/ = undef;
  $opfData = $json->decode(<$in>);
  close($in);
} else {
  die "Failed to open config file $dataFile\n$!";
}

$dataFile = catfile( updir(), "SCREEN", "screen.json" );
my $screenData;
if( open( my $in, "<", $dataFile ))
{ 
  local $/ = undef;
  $screenData = $json->decode(<$in>);
  close($in);
} else { 
  die "Failed to open config file $dataFile\n$!";
} 

# Done loading JSON of previous sections

$dataFile = "cls.json";
my $clsData = {};
if( open( my $in, "<", $dataFile ))
{ 
  local $/ = undef;
  $clsData = $json->decode(<$in>);
  close($in);
} 

#### 
# Two components of CLS, W and V_X

if( $commonOceanData->{'cls'}->{'W'} ) {
  runVOffset(  $commonOceanData, $dftData, $clsData);
  runWOffset(  $commonOceanData, $screenData, $clsData);
}
else  {
  runVOffset(  $commonOceanData, $dftData, $clsData);
  runWOffset(  $commonOceanData, $screenData, $clsData);
}


exit 0;
################
################

# 
sub runVOffset
{
  my ( $cod, $dft, $cls ) = @_;

  my $t0 = [gettimeofday];
  print "V offset\n";

  $cls->{'V'} = {} unless exists $cls->{'V'};
  $cls->{'V'}->{'site'} = {} unless exists $cls->{'V'}->{'site'};
  $cls->{'V'}->{'edge'} = {} unless exists $cls->{'V'}->{'edge'};


  my @ibe = indxByElement( $cod );
#  foreach ( indxByElement( $cod ) ) {
#    print "$_ \n";
#  }

  # Do we need to calc total potential?
  # 1. Do the previous calculations exist for these atomic sites?
  # 2. Did the SCF change?

  # List of sites (ignores edge)
  my %uniqueAtomicSite;
  foreach my $edge (@{$cod->{'calc'}->{'edges'}}) {
    print "$edge \n";
    my $s = (split ' ', $edge )[0];
    $uniqueAtomicSite{ $s } = 1;  # First element of the string
#    $uniqueAtomicSite{ @{ split ' ', $edge }[0] } = 1;  # First element of the string
  }

  my @runVtot;
  foreach my $site (keys %uniqueAtomicSite) { 
    my ( $el, $i ) = split ' ', $ibe[$site-1];
    $cls->{'V'}->{'site'}->{$el} = {} unless exists( $cls->{'V'}->{'site'}->{$el} );
    if( exists $cls->{'V'}->{'site'}->{$el}->{$i} ) {
      unless( $dft->{'scf'}->{'hash'} eq $cls->{'V'}->{'site'}->{$el}->{$i}->{'SCF hash'} ) {
        push @runVtot, $site;
        print "$site 0 \n";
      }
    } else {
      push @runVtot, $site;
        print "$site 1 \n";
    }
  }

  if( scalar @runVtot > 0 ) {
    unless( -d "pot" ) {
      mkdir "pot" or die "$!";
    }

    copy( catfile( updir(), "DFT", "potofr" ), catfile( "pot", "rhoofr" ) ) or die $!;
    copy( catfile( updir(), "DFT", "nfft.pot" ), catfile( "pot", "nfft" ) ) or die $!;

    chdir "pot" or die "$!";

    writeAvecs( $cod->{'structure'} );
    writeBvecs( $cod->{'structure'} );
    writeSitelistNew(  $cod, \@runVtot );
    
#    writeSitelist( $screen->{'general'} );
#    writeXYZ( $screen->{'general'} );

    system("$ENV{'OCEAN_BIN'}/rhoofg.x") == 0  or die "Failed to run rhoofg.x\n";
    system("wc -l rhoG2 > rhoofg") == 0 or die "$!\n";
    system("sort -n -k 6 rhoG2 >> rhoofg") == 0 or die "$!\n";

    open OUT, ">", "avg.ipt" or die $!;
    print OUT "500 0.01\n";
    close OUT;

    print "$cod->{'computer'}->{'para_prefix'} $ENV{'OCEAN_BIN'}/mpi_avg.x > mpi_avg.log 2>&1\n";
    system("$cod->{'computer'}->{'para_prefix'} $ENV{'OCEAN_BIN'}/mpi_avg.x > mpi_avg.log 2>&1" );
    if ($? == -1) {
        print "failed to execute: $!\n";
        die;
    }
    elsif ($? & 127) {
        printf "mpi_avg died with signal %d, %s coredump\n",
        ($? & 127),  ($? & 128) ? 'with' : 'without';
        die;
    }
    else {
      my $errorCode = $? >> 8;
      if( $errorCode != 0 ) {
        die "CALCULATION FAILED\n  mpi_avg exited with value $errorCode\n";
      }
      else {
        printf "mpi_avg exited successfully with value %d\n", $errorCode;
      }
    }

    foreach my $site (@runVtot) {
      my ( $el, $i ) = split ' ', $ibe[$site-1];
      $cls->{'V'}->{'site'}->{$el}->{$i}->{'SCF hash'} = $dft->{'scf'}->{'hash'};
    }

    my $jsonFile = catfile( updir(), "cls.json" );
    open OUT, ">", $jsonFile or die "$!";
    print OUT $json->encode( $cls );
    close OUT;

    chdir updir();
  }


  # Now that the real-space projection is done, make a complete list
  #  (In the future when multiple edges (e.g. K + L23 ) are supported)
  @runVtot = ();
  foreach my $site (@{$cod->{'calc'}->{'edges'}}) {
    my ($i, $n, $l) = split ' ', $site;
    my ($el, $j ) = split ' ', $ibe[$i-1];
    my $z = $cod->{'structure'}->{'znucl'}[$cod->{'structure'}->{'typat'}[$i-1]-1];
    my $nl = sprintf "%1i%1s", $n, $spdf[$l];
    $cls->{'V'}->{'edge'}->{$el} = {} unless( exists( $cls->{'V'}->{'edge'}->{$el} ) );
    $cls->{'V'}->{'edge'}->{$el}->{$j} = {} unless( exists( $cls->{'V'}->{'edge'}->{$el}->{$j} ) );
    $cls->{'V'}->{'edge'}->{$el}->{$j}->{$nl} = {} unless( exists( $cls->{'V'}->{'edge'}->{$el}->{$j}->{$nl} ) );
    my $s = sprintf "%2s %3i %1i %1i %04i", $el, $z, $n, $l, $j;
#    if( exists $cls->{'V'}->{'edge'}->{$el}->{$j}->{$nl}  ) {
      unless( $dft->{'scf'}->{'hash'} eq $cls->{'V'}->{'edge'}->{$el}->{$j}->{$nl}->{'SCF hash'} ) {
        push @runVtot, $s;
      }
#    } else {
#      push @runVtot, $s;
#    }
  }

  foreach (@runVtot) {
    print $_ . "\n";
  }


  if( scalar @runVtot > 0 ) {

    my $pot_factor = 1;
    if( $cod->{'general'}->{'program'} eq 'abi' ) {
      $pot_factor = 2;
    }

    my @Vtot = projectVtot( \@runVtot );

    for( my $i = 0; $i < scalar @Vtot; $i ++ ) {
      my @hfin = split ' ', $runVtot[$i];
      my $el = $hfin[0];
      my $j = $hfin[4];
      my $n = $hfin[2];
      my $l = $hfin[3];
      my $nl = sprintf "%1i%1s", $n, $spdf[$l];
      $cls->{'V'}->{'edge'}->{$el}->{$j}->{$nl}->{'SCF hash'} = $dft->{'scf'}->{'hash'};
      $cls->{'V'}->{'edge'}->{$el}->{$j}->{$nl}->{'pot'} = $Vtot[$i] * $pot_factor;

    }

    my $jsonFile = "cls.json" ;
    open OUT, ">", $jsonFile or die "$!";
    print OUT $json->encode( $cls );
    close OUT; 

  }
}

sub runWOffset
{
  my ( $cod, $screen, $cls ) = @_;

  my $t0 = [gettimeofday];
  print "W offset\n";

  $cls->{'W'} = {} unless exists $cls->{'W'};

  my @ibe = indxByElement( $cod );

  my @runWtot;
  foreach my $site (@{$cod->{'calc'}->{'edges'}}) {
    my ($i, $n, $l) = split ' ', $site;
    my ($el, $j ) = split ' ', $ibe[$i-1];
    my $z = $cod->{'structure'}->{'znucl'}[$cod->{'structure'}->{'typat'}[$i-1]-1];
    my $nl = sprintf "%1i%1s", $n, $spdf[$l];

    $cls->{'W'}->{$el} = {} unless( exists $cls->{'W'}->{$el});
    $cls->{'W'}->{$el}->{$j} = {} unless( exists $cls->{'W'}->{$el}->{$j});
    $cls->{'W'}->{$el}->{$j}->{$nl} = {} unless( exists $cls->{'W'}->{$el}->{$j}->{$nl});
    my $s = sprintf "%2s %3i %1i %1i %04i", $el, $z, $n, $l, $j;
#    my $s = sprintf "%2s%04i", $el, $j;
    my $clean = 0;

    if( $screen->{'screen'}->{'hash'} eq $cls->{'W'}->{$el}->{$j}->{$nl}->{'hash'} ) {
      if( exists $cls->{'W'}->{$el}->{$j}->{$nl}->{'pot'} ) {
        foreach my $r (@{$screen->{'screen'}->{'shells'}}) {
          my $rad = sprintf "03.2f", $r;
          unless( exists $cls->{'W'}->{$el}->{$j}->{$nl}->{'pot'}->{$rad} ) {
            $clean = 1;
            last;
          }
        }
      } else {
        $clean = 1;
      }
    } else {
      $clean = 1;
    }
    if( $clean ) {
      $cls->{'W'}->{$el}->{$j}->{$nl}->{'pot'} = {};
      push @runWtot, $s;
    }
  }

  return if( scalar @runWtot == 0 );
  my @rad;
  foreach my $r (@{$screen->{'screen'}->{'shells'}}) {
    push @rad,sprintf( "%03.2f", $r);
  }

  foreach my $s (@runWtot) {
    my @edge = split ' ', $s;
    my $el = $edge[0];
    my $j = $edge[4];
    my $nl = sprintf "%1i%1s", $edge[2], $spdf[$edge[3]];
    my @W = projectW( $s, \@rad );
    for( my $i = 0; $i < scalar @rad; $i++) {
      $cls->{'W'}->{$el}->{$j}->{$nl}->{'pot'}->{$rad[$i]} = $W[$i]
    }
    $cls->{'W'}->{$el}->{$j}->{$nl}->{'hash'} = $screen->{'screen'}->{'hash'};
  }
  
  my $jsonFile = "cls.json" ;
  open OUT, ">", $jsonFile or die "$!";
  print OUT $json->encode( $cls );
  close OUT; 

}


sub writeAvecs
{
  my $structureRef = $_[0];

  open OUT, ">", "avecsinbohr.ipt" or die "Failed to open avecsinbohr.ipt\n$!";
  for( my $i = 0; $i < 3; $i++ ) {
    printf  OUT "%.16g  %.16g  %.16g\n", $structureRef->{'avecs'}[$i][0],
                                $structureRef->{'avecs'}[$i][1],
                                $structureRef->{'avecs'}[$i][2];

  }
  close OUT;
}

sub writeBvecs
{
  my $structureRef = $_[0];
  open OUT, ">", "bvecs" or die "Failed to open bvecs\n$!";
  for( my $i = 0; $i < 3; $i++ ) {
    printf  OUT "%.16g  %.16g  %.16g\n", $structureRef->{'bvecs'}[$i][0],
                                $structureRef->{'bvecs'}[$i][1],
                                $structureRef->{'bvecs'}[$i][2];

  }
  close OUT;
}

sub writeSitelistNew #(  $cod->{'structure'}, \@runVsite );
{
  my ( $genRef, @sitelist ) = ( $_[0], @{$_[1]} );

  my @siteCounter;
  open OUT, ">", "sitelist.new" or die "$!";
  print OUT scalar @sitelist . "\n";
  my @output;
  my $natom = scalar @{$genRef->{'structure'}->{'typat'}};
  for ( my $i = 0; $i < $natom; $i ++ ) {
    my $t = $genRef->{'structure'}->{'typat'}[$i];
    $siteCounter[$t] ++ ;
    push @output, sprintf "%2s %8i %.16g %.16g %.16g\n", $genRef->{'structure'}->{'elname'}[$t-1], 
                $siteCounter[$t], $genRef->{'structure'}->{'xred'}[$i][0], 
                $genRef->{'structure'}->{'xred'}[$i][1], $genRef->{'structure'}->{'xred'}[$i][2];
  }

  foreach my $s (@sitelist ) {
    print OUT $output[$s-1];
  }
  close OUT;
}

sub indxByElement
{
  my ( $genRef ) = @_;

  my @indxByEl;
  my @siteCounter;
  my $natom = scalar @{$genRef->{'structure'}->{'typat'}};
  for ( my $i = 0; $i < $natom; $i ++ ) {
    my $t = $genRef->{'structure'}->{'typat'}[$i];
    $siteCounter[$t] ++ ;
    push @indxByEl, sprintf "%2s %04i", $genRef->{'structure'}->{'elname'}[$t-1],
                $siteCounter[$t];
  }

  return @indxByEl;
}

sub writeXYZ 
{
  my $genRef = $_;
  open OUT, ">", "xyz.wyck" or die $!;
  print OUT scalar @{$genRef->{'wyck'}};
  print OUT "\n";
  foreach (@{$genRef->{'wyck'}}) {
    print OUT $_ . "\n";
  }
  close OUT;
}

sub writeSitelist
{
  my $genRef = $_;
  open OUT, ">", "sitelist" or die $!;
  print OUT scalar @{$genRef->{'sitelist'}};
  print OUT "\n";
  foreach (@{$genRef->{'sitelist'}}) {
    print OUT $_ . "\n";
  }
  close OUT;
}

sub projectVtot
{
  my @edges = @{$_[0]};

  my @V;
  foreach my $edge ( @edges ) {
    print "$edge \n";

    my @hfin = split ' ', $edge;
    my $cf = sprintf "coreorbz%03in%02il%02i", $hfin[1], $hfin[2], $hfin[3];
    my $coreFile = catfile( updir(), "OPF", "zpawinfo", $cf );
    open IN, "<", $coreFile or die "Failed to open $coreFile\n$!";
    <IN>;
    my @rad;
    my @wvfn;
    while ( my $line = <IN> )
    {
      $line =~ m/(\S+)\s+(\S+)/;
      push @rad, $1;
      push @wvfn, $2;
    }
    close IN;

    my $rat = $rad[-1]/$rad[0];
    my $dl = log( $rat )/ ($#rad);
    my $xrat = exp( $dl );

    my $xr1 = sqrt( $xrat ) - sqrt( 1.0/$xrat);
    my $rmin = $rad[0] / $xrat;


    my $sum = 0;
    my $sum2 = 0;
    for( my $i = 0; $i <= $#rad; $i++ )
    {
      my $temp_rad = $rmin * $xrat**($i+1);
      $sum += $rad[$i] * $rad[$i] * $rad[$i] * $xr1 * $wvfn[$i]**2;
      $sum2 += $rad[$i] * $xr1 * $wvfn[$i]**2;
    }
    print "$sum2\n";

    my $filename = catfile( "pot", sprintf("avg%2s%04i", $hfin[0], $hfin[4]) );
    open IN, $filename or die "Failed to open $filename\n$!";
    my @prad;
    my @pot;
    while( my $line = <IN> )
    {
      $line =~ m/\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+/;
      push @prad, $1;
      push @pot, $2;
    }
    close IN;

    $sum = 0;
    my $j = 0;
    for( my $i = 0; $i <= $#rad; $i++ )
    {
      while( $rad[$i] > $prad[$j] )
      {
        die if( $j > $#prad );
        $j++;
      }
      my $interp = @pot[$j];
      if( $j > 0 && $j < $#prad )
      {
        my $run = $prad[$j+1]-$prad[$j];
        my $slope = ($pot[$j+1]-$pot[$j]) / $run;
        $interp += $slope * ( $rad[$i] - $prad[$j] );
      }
      $sum += $rad[$i] * $xr1 * $wvfn[$i]**2 * $interp;
    }

    push @V, $sum;
  }

  return @V;
}

sub projectW
{
  my $edge = $_[0];
  my @shells = @{$_[1]};

  my @Wpot;
  
  my @hfin = split ' ', $edge;
  my $el = $hfin[0];
  my $z = $hfin[1];
  my $n = $hfin[2];
  my $l = $hfin[3];
  my $indx = $hfin[4];
  
  my $cf = sprintf "coreorbz%03in%02il%02i", $z, $n, $l;
  my $coreFile = catfile( updir(), "OPF", "zpawinfo", $cf );
  open IN, "<", $coreFile or die "Failed to open $coreFile\n$!";
  <IN>;
  my @rad;
  my @wvfn;
  while ( my $line = <IN> )
  {
    $line =~ m/(\S+)\s+(\S+)/;
    push @rad, $1;
    push @wvfn, $2;
  }
  close IN;

  my $rat = $rad[-1]/$rad[0];
  my $dl = log( $rat )/ ($#rad);
  my $xrat = exp( $dl );
  my $xr1 = sqrt( $xrat ) - sqrt( 1.0/$xrat);
  my $rmin = $rad[0] / $xrat;

  my $zel = sprintf "z%s%04d", $el, $indx;
  my $znl = sprintf "n%02il%02i", $n, $l;


  foreach my $radius (@shells) {
    my $rad_dir = sprintf("zR%03.2f", $radius );
    my $filename = catfile( updir(), 'SCREEN', $zel, $znl, $rad_dir, 'ropt' );
    open IN, $filename or die "Failed to open $filename\n$!";


    my @prad;
    my @pot;
    while( my $line = <IN> )
    {
      $line =~ m/(\S+)\s+\S+\s+\S+\s+(\S+)/;
      push @prad, $1;
      push @pot, $2;
    }
    close IN;

    my $sum = 0;
    my $j = 0;
    for( my $i = 0; $i <= $#rad; $i++ )
    {
      while( $rad[$i] > $prad[$j] )
      {
        die if( $j > $#prad );
        $j++;
      }
      my $interp = @pot[$j];
      if( $j > 1 && $j < $#prad - 1 )
      {
        $interp = $pot[$j-2] * ($rad[$i]-$prad[$j-1])*($rad[$i]-$prad[$j])*($rad[$i]-$prad[$j+1])
                            / (($prad[$j-2]-$prad[$j-1])*($prad[$j-2]-$prad[$j])*($prad[$j-2]-$prad[$j+1]) )
                + $pot[$j-1] * ($rad[$i]-$prad[$j])*($rad[$i]-$prad[$j+1])*($rad[$i]-$prad[$j-2])
                            / (($prad[$j-1]-$prad[$j])*($prad[$j-1]-$prad[$j+1])*($prad[$j-1]-$prad[$j-2]) )
                + $pot[$j]   * ($rad[$i]-$prad[$j-1])*($rad[$i]-$prad[$j+1])*($rad[$i]-$prad[$j-2])
                            / (($prad[$j]-$prad[$j-1])*($prad[$j]-$prad[$j+1])*($prad[$j]-$prad[$j-2]) )
                + $pot[$j+1] * ($rad[$i]-$prad[$j-1])*($rad[$i]-$prad[$j])*($rad[$i]-$prad[$j-2])
                            / (($prad[$j+1]-$prad[$j-1])*($prad[$j+1]-$prad[$j])*($prad[$j+1]-$prad[$j-2]) );
      }
      elsif( $j == $#prad - 1 || $j == 1 || $j == 0 )
      {
        my $run = $prad[$j+1]-$prad[$j];
        my $slope = ($pot[$j+1]-$pot[$j]) / $run;
        $interp += $slope * ( $rad[$i] - $prad[$j] );
      }
      $sum += $rad[$i] * $xr1 * $wvfn[$i]**2 * $interp;
    }
    push @Wpot, $sum;
  }

  return @Wpot;
}
  

