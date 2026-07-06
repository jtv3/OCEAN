#!/usr/bin/perl
# Copyright (C) 2021 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
use strict;
require JSON::PP;
JSON::PP->import;
use File::Copy;


my $input_filename = $ARGV[0];
my $config_filename = $ARGV[1];
my $type_filename = $ARGV[2];

unless(  -e $input_filename )
{
  die "Could not find file $input_filename!\n";
}

unless( -e $config_filename )
{
  die "Could not find config file $config_filename!\n";
}

my $json = JSON::PP->new;
my $config;
if( open( my $in, "<", $config_filename ))
{
  local $/ = undef;
  $config = $json->decode(<$in>);
  close($in);
}
else
{
  die "Failed to open config file $config_filename\n$!";
}

my $typeDef;
if( open( my $in, "<", $type_filename ) )
{ 
  local $/ = undef;
  $typeDef = $json->decode(<$in>);
  close($in);
}
else
{ 
  die "Failed to open $type_filename\n$!";
}

validateTypeSpecs( $typeDef, '' );




my %decoder = (
  'calc' => 'calc.mode',
  'dft' => 'dft.program',
  'nkpt' => 'bse.kmesh',
  'ngkpt' => 'dft.den.kmesh',
  'ngkpt.auto' => 'dft.den.auto',
  'photon_q' => 'calc.photon_q',
  'dft.split' => 'dft.bse.split',
  'dft.qe_redirect' => 'dft.redirect',
  'nbands' => 'bse.nbands',
  'dft_energy_range' => 'bse.dft_energy_range',
  'obf_energy_range' => 'nope.obf_energy_range',
  'obkpt' => 'nope.obkpt',
  'obf.nbands' => 'nope.obf_nbands',
  'trace_tol' => 'nope.trace_tol',
  'acc_level' => 'nope.acc_level',
  'k0' => 'bse.kshift',
  'fband' => 'dft.fband',
  'occopt' => 'dft.occopt',
  'mixing' => 'dft.mixing',
  'acell' => 'structure.rscale',
  'rprim' => 'structure.rprim',
  'ntypat' => 'nope.ntype',
  'typat' => 'structure.typat',
  'znucl' => 'structure.znucl',
  'zsymb' => 'structure.zsymb',
  'pp_list' => 'psp.pp_list',
  'pp_database' => 'psp.pp_database',
  'ecut.quality' => 'psp.ecut_quality',
  'natom' => 'nope.natoms',
  'coord' => 'nope.structure.coord',
  'xred' => 'structure.xred',
  'ecut' => 'dft.ecut',
  'diemac' => 'structure.epsilon',
  'toldfe' => 'dft.toldfe',
  'tolwfr' => 'dft.tolwfr',
  'nstep' => 'dft.nstep',
  'dft.startingwfc' => 'dft.startingwfc',
  'dft.diagonalization' => 'dft.diagonalization',
  'dft.ndiag' => 'dft.ndiag',
  'dft.functional' => 'dft.functional',
  'dft.exx.qmesh' => 'dft.exx.qmesh',
  'dft.nscf.poolsize' => 'dft.bse.poolsize',
  'verbatim' => 'dft.verbatim',
  'para_prefix' => 'computer.para_prefix',
  'ser_prefix' => 'computer.ser_prefix',
  'abpad' => 'dft.abpad',
  'scfac' => 'bse.core.scfac',
  'screen.shells' => 'screen.shells',
  'opf.hfkgrid' => 'opf.shirley.hfkgrid',
  'opf.fill' => 'opf.shirley.fill',
  'opf.opts' => 'opf.shirley.opts',
  'opf.program' => 'opf.program',
  'screen.nkpt' => 'screen.kmesh',
  'screen.k0' => 'screen.kshift',
  'screen.nbands' => 'screen.nbands',
  'caution' => 'opf.shirley.caution',
  'nedges' => 'nope.nedges',
  'edges' => 'calc.edges',
  'cnbse.nbuse' => 'nope.nbuse',
  'cnbse.xmesh' => 'bse.xmesh',
  'cnbse.rad' => 'bse.core.screen_radius',
  'metal' => 'structure.metal',
  'cksshift' => 'nope.cksshift',
  'cksstretch' => 'nope.cksstretch',
  'cnbse.niter' => 'bse.core.haydock.niter',
  'haydock_convergence' => 'bse.core.haydock.converge.thresh',
  'cnbse.spect_range' => 'bse.core.plot.range',
  'cnbse.broaden' => 'bse.core.broaden',
  'cnbse.strength' => 'bse.core.strength',
  'cnbse.solver' => 'bse.core.solver',
  'cnbse.gmres.elist' => 'bse.core.gmres.elist',
  'cnbse.gmres.erange' => 'bse.core.gmres.erange',
  'cnbse.gmres.nloop' => 'bse.core.gmres.nloop',
  'cnbse.gmres.gprc' => 'bse.core.gmres.gprc',
  'cnbse.gmres.ffff' => 'bse.core.gmres.ffff',
  'cnbse.write_rhs' => 'bse.core.write_rhs',
  'cnbse.gw.control' => 'bse.core.gw.control',
  'bse.gw.cstr' => 'nope.bse_gw_cstr',
  'bse.gw.vstr' => 'nope.bse_gw_vstr',
  'bse.gw.gap' => 'nope.gwgap',
  'degauss' => 'dft.degauss',
  'ibrav' => 'nope.ibrav',
  'isolated' => 'nope.isolated',
  'noncolin' => 'dft.noncolin',
  'prefix' => 'nope.prefix',
  'ppdir' => 'psp.ppdir',
  'dft.calc_stress' => 'dft.calc_stress',
  'dft.calc_force' => 'dft.calc_force',
  'spinorb' => 'dft.spinorb',
  'work_dir' => 'nope.wordir',
  'tmp_dir' => 'dft.tmp_dir',
  'den.kshift' => 'dft.den.kshift',
  'core_offset' => 'screen.core_offset.enable',
  'ham_kpoints' => 'nope.ham_kpoints',
  'nbse.niter' => 'bse.val.haydock.niter',
  'nbse.backf' => 'bse.val.backf',
  'nbse.aldaf' => 'bse.val.aldaf',
  'nbse.qpflg' => 'bse.val.qpflg',
  'nbse.bwflg' => 'bse.val.bwflg',
  'nbse.bande' => 'bse.val.bande',
  'nbse.bflag' => 'bse.val.bflag',
  'nbse.lflag' => 'bse.val.lflag',
  'nbse.convergence' => 'nope.convergence',
  'nbse.decut' => 'bse.val.decut',
  'nbse.se_rs' => 'nope.se.rs',
  'nbse.se_metal' => 'nope.se.metal',
  'nbse.se_niter' => 'nope.se.niter',
  'nbse.spect_range' => 'bse.val.plot.range',
  'tot_charge' => 'dft.tot_charge',
  'nspin' => 'dft.nspin',
  'smag' => 'dft.smag',
  'ldau' => 'dft.ldau.Hubbard_U',
  'qe_scissor' => 'nope.qe_scissor',
  'nphoton' => 'nope.nphoton',
  'ser_bse' => 'nope.ser_bse',
  'spin_orbit' => 'bse.core.spin_orbit',
  'screen_energy_range' => 'screen.dft_energy_range',
  'screen.grid.scheme' => 'screen.grid.scheme',
  'screen.grid.rmode' => 'screen.grid.rmode',
  'screen.grid.ninter' => 'nope.screen_grid_ninter',
  'screen.grid.shells' => 'screen.grid.shells',
  'screen.grid.xyz' => 'nope.screen_grid_xyz',
  'screen.grid.rmax' => 'screen.grid.rmax',
  'screen.grid.nr' => 'nope.screen_grid_nr',
  'screen.grid.ang' => 'screen.grid.ang',
  'screen.grid.deltar' => 'screen.grid.deltar',
  'screen.lmax' => 'screen.grid.lmax',
  'screen.grid.nb' => 'nope.screen_grid_nb',
  'screen.final.rmax' => 'screen.final.rmax',
  'screen.final.dr' => 'screen.final.dr',
  'screen.model.dq' => 'screen.model.SLL.dq',
  'screen.model.qmax' => 'screen.model.SLL.qmax',
  'screen.legacy' => 'nope.screen_legacy',
  'screen.augment' => 'screen.augment',
  'screen.wvfn' => 'nope.screen_wvfn',
  'screen.convertstyle' => 'screen.convertstyle',
  'screen.inversionstyle' => 'screen.inversionstyle',
  'screen.mode' => 'screen.mode',
  'bse.wvfn' => 'nope.bse_wvfn',
  'hamnum' => 'nope.hamnum',
  'echamp' => 'bse.core.gmres.echamp',
  'bshift' => 'nope.bshift' );


open IN, "<", $input_filename or die "Failed to open input file $input_filename\n$!\n";

my $rawInputFile = '';
my $inputString = '';
while( my $line = <IN> )
{
  $rawInputFile .= $line;
  # if there are comment characters -- #, *, or ! --
  #   remove them and everything following
  $line =~ s/[#\*\\!].*/ /;
  # just pad with spaces, not that inefficient
  $line =~ s/\{/ \{ /;
  $line =~ s/\}/ \} /;
  $inputString .= $line;
}
close IN;


my @inputFile = split ' ', $inputString;

if( 0 ){
foreach my $i (@inputFile)
{
  print "$i\n";
}
}

my %inputHash;
my @inputOrder;
my %seenRawKey;
my $i = 0;
my $tag = 1;
my $curly = 0;
my $key = '';
my $val = '';
my $errorBuffer = '';
while( $i < scalar @inputFile )
{
  if( $tag == 1 )
  {
    $errorBuffer .= $inputFile[$i] . "\n";
    die "Misplaced braces when expecting a tag\n>>>>\n$errorBuffer<<<<<\n" if( $inputFile[$i] =~ m/\{|\}/ );
#    print "$inputFile[$i] >>>> ";
    $tag = 0;
    $key = $inputFile[$i];
    if( exists $seenRawKey{ $key } )
    {
      die "Duplicate input flag: $key\n";
    }
    $seenRawKey{ $key } = 1;
    push @inputOrder, $key;
    $val = '';
    $errorBuffer = '';
    $errorBuffer .= $inputFile[$i-1] . "\n" if( $i > 0 );
    $errorBuffer .= $inputFile[$i] . "\n";
  }
  elsif ( $curly == 0 ) 
  {
    $errorBuffer .= $inputFile[$i] . "\n";
    if( $inputFile[$i] =~ m/\{/ )
    {
      $curly = 1;
    }
    elsif( $inputFile[$i] =~ m/\}/ )
    {
      die "Close brace when not expected\n>>>>\n$errorBuffer<<<<<\n";
    }
    else
    {
#      print "$inputFile[$i]\n";
      $tag = 1;
      $val .= $inputFile[$i] . " ";
    }
  } else
  {
    $errorBuffer .= $inputFile[$i] . "\n";
    die "Second open {\n>>>>\n$errorBuffer<<<<<\n" if( $inputFile[$i] =~ m/\{/ );
    if( $inputFile[$i] =~ m/\}/ )
    {
      $curly = 0;
      $tag = 1;
#      print "\n";
    }
    else
    {
      $val .= $inputFile[$i] . " ";
#      print "$inputFile[$i] ";
    }
  }
  $i++;
  # there will always be a trailing space
#  chop( $val );
  $inputHash{ $key } = $val;
}

my $haveLegacy = 0;
my $haveLegacyCorePlot = 0;
my $haveLegacyValPlot = 0;

INPUT: foreach my $key ( @inputOrder ) 
{
  unless( findInputKey( $config, $key ) )
  { 
    $haveLegacy = 1;
    print "Unrecognized input flag: $key\n  Attempting legacy conversion\n";
    last INPUT;
  }
}

my %seenInputKey;

if( $haveLegacy == 1 )
{
  if( exists $inputHash{ 'ppdir' } ) {
    if( $inputHash{ 'ppdir' } =~ m/^\s*'(.+)'\s*$/ ) {
      print $inputHash{ 'ppdir' };
      $inputHash{ 'ppdir' } = $1;
      print "  " .$inputHash{ 'ppdir' } . "\n";
    }
    if( $inputHash{ 'ppdir' } =~ m/^\.\.\/$/ ) {
      print $inputHash{ 'ppdir' };
      $inputHash{ 'ppdir' } = './';
      print "  " .$inputHash{ 'ppdir' } . "\n";
    }
  }
  foreach my $key ( @inputOrder )
  {
    my $lckey = lc($key);
#    $key = lc($key) unless( exists $decoder{$key} );
#    die "Unrecognized input flag: $key\n No recovery possible!" unless( exists $decoder{$lckey} );
    if( exists $decoder{$lckey} ) {
      my $newKey = $decoder{ $lckey };
      if( exists $seenInputKey{ $newKey } )
      {
        die "Duplicate input after legacy conversion: " 
          . $seenInputKey{ $newKey } . " and $key both set $newKey\n";
      }
      $seenInputKey{ $newKey } = $key;
      print "$key : $newKey  $inputHash{ $key }\n";
      $inputHash{ $newKey } = $inputHash{ $key };
      delete( $inputHash{ $key } ) unless( $key eq $newKey );
      $rawInputFile =~ s/$key/$newKey/;
      my @newKey = split /\./, $newKey;
      my $ref = $config;
      next if( $newKey[0] eq 'nope' );
      for( my $i = 0; $i < scalar @newKey; $i++ )
      {
        if( exists $ref->{$newKey[$i]} )
        {
          $ref = $ref->{$newKey[$i]};
        }
        else
        {
          die "Unrecognized input flag: $key\n  Legacy conversion failed\n";
        }
      }
      if( $lckey eq 'cnbse.spect_range' ) {
        $haveLegacyCorePlot = 1;
      } elsif ( $lckey eq 'nbse.spect_range' ) {
        $haveLegacyValPlot = 1;
      }
    } else {
      unless( findInputKey( $config, $key ) )
      {
        die "Unsupported input flag. Neither new nor legacy:  $key\n";
      }
      if( exists $seenInputKey{ $key } )
      {
        die "Duplicate input after legacy conversion: " 
          . $seenInputKey{ $key } . " and $key both set $key\n";
      }
      $seenInputKey{ $key } = $key;
      print "Comment: Mixed new and legacy input:  $key\n";
    }
  }

  my $newInputFile = $input_filename . ".mod3";
  open OUT, ">", $newInputFile or die "Failed to open $newInputFile \n$!";
  print OUT $rawInputFile;
  close OUT;
}

print "Storing parsed data\n\n";
# If we made it here all the keys are valid
my %suppliedInputKey;
foreach my $key ( keys %inputHash )
{
  my $value = $inputHash{ $key };
  print "$key $value\n";
  my @newKey = split /\./, $key;
  next if( $newKey[0] eq 'nope' );
  $suppliedInputKey{ $key } = 1;


  my $type = $typeDef;
  for( my $i = 0; $i < scalar @newKey; $i++ )
  {
    $type = $type->{$newKey[$i]} ;
  } 

  my $hashref = $config;
  for( my $i = 0; $i < scalar @newKey - 1; $i++ )
  {
    $hashref = $hashref->{$newKey[$i]};
  }

  my $regex;
  my ($baseType, $constraints) = parseTypeSpec( $type );
  $regex = '^\s*(-?\d+)\s*$' if( $baseType =~ m/i/ );
  # Full-token floating point match:
  #   ^\s* and \s*$ allow only optional leading/trailing whitespace.
  #   -? allows an optional minus sign.
  #   (?:\d+(?:\.\d*)?|\.\d+) accepts either digits with an optional decimal
  #     point and optional following digits, or a leading decimal point followed
  #     by digits. This allows 1, 1., 1.0, and .1, but rejects bare ".".
  #   (?:[eEdD][+-]?\d+)? accepts an optional Fortran/C exponent with e, E, d,
  #     or D, an optional sign, and at least one exponent digit.
  $regex = '^\s*(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?)\s*$' if( $baseType =~ m/f/ );
  $regex = '^([\w\S\s]+)$' if ( $baseType =~ m/s|S/ );


  # if array
  if( $baseType =~ m/a/ )
  {
    my @rawArray = split ' ', $value;
    foreach my $i (@rawArray)
    {
      die "Failed to match: $i of type $type\n" unless( $i =~ m/$regex/ );
    }
    if( $baseType =~ m/[if]/ )
    {
      for( my $i = 0; $i < scalar @rawArray; $i++ )
      { 
        $rawArray[$i] =~ s/[dD]/e/; # won't matter for int
        $rawArray[$i] *= 1 }
    }
    elsif( $baseType =~ m/s/ )
    {
      for( my $i = 0; $i < scalar @rawArray; $i++ )
      { $rawArray[$i] = lc $rawArray[$i] }
    }

    # If we did legacy translation, patch up the incompatibilities
    if( $key =~ m/bse.core.plot.range/ )
    {
      if( $haveLegacyCorePlot ) 
      {
        my $points = shift @rawArray;
        $config->{'bse'}->{'core'}->{'plot'}->{'points'} = $points if( $key =~ m/core/ );
      }
    } elsif( $key =~ m/bse.val.plot.range/ ) 
    {
      if( $haveLegacyValPlot ) {
        my $points = shift @rawArray;
        $config->{'bse'}->{'val'}->{'plot'}->{'points'} = $points if( $key =~ m/val/ );
      }
    } elsif( $haveLegacy == 1 && $key =~m/pp_list/ )
    {
      $config->{'psp'}->{'source'} = 'manual' if( scalar @rawArray > 0 );
    }
    # End legacy fix
    validateArrayLength( $key, $type, $constraints, \@rawArray );
    $hashref->{$newKey[-1]} = [@rawArray];
  }
  else
  {
#    if( $key =~ m/tol/ ) {
#      print "STOP" . $value . "\n" . $regex . "\n";
#      $value =~ s/d/e/;
#      $value =~ m/$regex/;
#      print $1 . "\n";
#      print $1*1 . "\n";
#      exit 1;
#    }
    unless( $value eq ' ' )
    {
      $value =~ s/^\s+//;
      $value =~ s/\s+$//;
    }
    # If we did legacy translation, patch up the incompatibilities
    if( $key =~ m/structure.epsilon/ )
    {
      if( $value =~ m/dfpt/ )
      {
         $value = 0;
         $config->{'dft'}->{'epsilon'}->{'method'} = 'dfpt';
      }
      else
      {
         $config->{'dft'}->{'epsilon'}->{'method'} = 'input';
      }
    }
    elsif( $key =~ m/core_offset/ )
    {
      if( $value =~ m/\d/ )
      {
        $config->{'screen'}->{'core_offset'}->{'energy'} = $value;
        $value = 'true';
      }
    }
    # end fix  
    if( $baseType =~ m/b/ )
    {
      my $boolValue = lc $value;
      if( $boolValue eq 't' || $boolValue eq 'true' || $boolValue eq '.t.'
          || $boolValue eq '.true.' || $boolValue eq '1' )
      {
        $value = $JSON::PP::true
      }
      elsif( $boolValue eq 'f' || $boolValue eq 'false' || $boolValue eq '.f.'
             || $boolValue eq '.false.' || $boolValue eq '0' )
      {
        $value = $JSON::PP::false
      }
      else
      {
        die "Failed to match: $value of type $type\n";
      }
    }
    else
    {
      if( $value =~ m/$regex/ )
      {
        $value = $1;
        $value =~ s/[dD]/e/ if( $baseType =~ m/f/ );
        $value *= 1 if( $baseType =~ m/[if]/ );
        $value = lc $value if( $baseType =~ m/s/ );
      }
      else
      {
        die "Failed to match: $value of type $type\n";
      }
    }
    $hashref->{$newKey[-1]} = $value;
  }
}

validateRelationalArrayLengths( $config, $typeDef, \%suppliedInputKey );

my $enable = 1;
$json->canonical([$enable]);
$json->pretty([$enable]);
open OUT, ">", "parsedInputFile" or die "Failed to open parsedInputFile\n$!";
print OUT $json->encode($config);
close OUT;


copy( "parsedInputFile", "oceanDatafile") ;


##### REMOVE IN FUTURE
open OUT, ">", "dft" or die;
print OUT $config->{'dft'}->{'program'} . "\n";
close OUT;

open OUT, ">", "calc" or die;
print OUT $config->{'calc'}->{'mode'} . "\n";
close OUT;



sub findInputKey
{
  my ($config, $key) = @_;
  my @newKey = split /\./, $key;
  my $ref = $config;
  for( my $i = 0; $i < scalar @newKey; $i++ )
  {
    return 0 unless( ref( $ref ) eq 'HASH' && exists $ref->{$newKey[$i]} );
    $ref = $ref->{$newKey[$i]};
  }
  return 1;
}


# Split an oparse.type.json leaf into the original compact base type and any
# comma-separated constraint suffixes. For example:
#   af:len=3
# becomes base type "af" and constraints { len => 3 }.
sub parseTypeSpec
{
  my ($type) = @_;
  my ($baseType, $constraintString) = split /:/, $type, 2;
  my %constraints;

  if( defined $constraintString && length $constraintString )
  {
    foreach my $constraint ( split /,/, $constraintString )
    {
      $constraint =~ m/^(\w+)=(.+)$/
        or die "Malformed type constraint '$constraint' in $type\n";
      my $name = $1;
      my $value = $2;

      if( $name eq 'len' || $name eq 'signlen' )
      {
        die "$name constraint '$constraint' in $type requires an array type\n"
          unless( $baseType =~ m/a/ );
        $value =~ m/^\d+$/
          or die "Invalid numeric type constraint '$constraint' in $type\n";
        die "Invalid signlen constraint '$constraint' in $type\n"
          if( $name eq 'signlen' && $value == 0 );
        $constraints{$name} = $value * 1;
      }
      elsif( $name eq 'oneof' || $name eq 'lenmatch' )
      {
        die "$name constraint '$constraint' in $type requires an array type\n"
          unless( $baseType =~ m/a/ );
        $value =~ m/^[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*$/
          or die "Invalid path type constraint '$constraint' in $type\n";
        $constraints{$name} = $value;
      }
      elsif( $name eq 'lenmul' )
      {
        die "$name constraint '$constraint' in $type requires an array type\n"
          unless( $baseType =~ m/a/ );
        $value =~ m/^[1-9]\d*\*[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*$/
          or die "Invalid lenmul type constraint '$constraint' in $type\n";
        $constraints{$name} = $value;
      }
      else
      {
        die "Unknown type constraint '$name' in $type\n";
      }
    }
  }

  return ($baseType, \%constraints);
}


# Walk the full type tree once at startup so malformed constraints in
# oparse.type.json fail even when the associated input is not supplied.
sub validateTypeSpecs
{
  my ($typeRef, $prefix) = @_;
  if( ref( $typeRef ) eq 'HASH' )
  {
    foreach my $key ( keys %$typeRef )
    {
      my $newPrefix = length $prefix ? "$prefix.$key" : $key;
      validateTypeSpecs( $typeRef->{$key}, $newPrefix );
    }
  }
  else
  {
    parseTypeSpec( $typeRef );
  }
}


# Enforce constraints that can be checked for one parsed array at a time.
# This runs after legacy compatibility fixups, so legacy plot ranges have
# already had the old leading "points" value removed before len=2 is checked.
sub validateArrayLength
{
  my ($key, $type, $constraints, $arrayRef) = @_;
  my $length = scalar @$arrayRef;

  if( exists $constraints->{'len'} )
  {
    my $expected = $constraints->{'len'};
    die "Invalid array length for $key: got $length values, expected len=$expected ($type)\n"
      unless( $length == $expected );
  }

  if( exists $constraints->{'signlen'} )
  {
    my $positiveLength = $constraints->{'signlen'};
    die "Invalid array length for $key: got 0 values, expected signlen=$positiveLength ($type)\n"
      if( $length == 0 );

    if( $arrayRef->[0] < 0 )
    {
      die "Invalid array length for $key: got $length values, expected 1 because first value is negative ($type)\n"
        unless( $length == 1 );
    }
    elsif( $arrayRef->[0] > 0 )
    {
      die "Invalid array length for $key: got $length values, expected signlen=$positiveLength because first value is positive ($type)\n"
        unless( $length == $positiveLength );
    }
    else
    {
      die "Invalid array length for $key: first value is zero, expected positive or negative first value ($type)\n";
    }
  }
}


# Fetch a dot-separated path from the final config tree.
sub getConfigValue
{
  my ($config, $key) = @_;
  my @path = split /\./, $key;
  my $ref = $config;
  foreach my $part ( @path )
  {
    return undef unless( ref( $ref ) eq 'HASH' && exists $ref->{$part} );
    $ref = $ref->{$part};
  }
  return $ref;
}


# Gather constraints that need final config values instead of one raw input array.
sub collectRelationalArrayConstraints
{
  my ($typeRef, $prefix, $oneOfRef, $checksRef) = @_;

  if( ref( $typeRef ) eq 'HASH' )
  {
    foreach my $key ( sort keys %$typeRef )
    {
      my $newPrefix = length $prefix ? "$prefix.$key" : $key;
      collectRelationalArrayConstraints( $typeRef->{$key}, $newPrefix, $oneOfRef, $checksRef );
    }
  }
  else
  {
    my (undef, $constraints) = parseTypeSpec( $typeRef );
    push @{ $oneOfRef->{ $constraints->{'oneof'} } }, $prefix
      if( exists $constraints->{'oneof'} );
    push @$checksRef, { key => $prefix, type => $typeRef, constraints => $constraints }
      if( exists $constraints->{'lenmul'} || exists $constraints->{'lenmatch'} );
  }
}


# Enforce relational constraints after all supplied values have been assigned.
# This keeps checks independent of the order in the user's input file.
sub validateRelationalArrayLengths
{
  my ($config, $typeDef, $suppliedRef) = @_;
  my %oneOf;
  my @checks;

  collectRelationalArrayConstraints( $typeDef, '', \%oneOf, \@checks );

  foreach my $group ( sort keys %oneOf )
  {
    my @members = @{ $oneOf{$group} };
    my @supplied = grep { exists $suppliedRef->{$_} } @members;
    die "Invalid array selection for $group: got " . scalar @supplied
      . " supplied inputs, expected exactly one of " . join( ', ', @members ) . "\n"
      unless( scalar @supplied == 1 );
  }

  foreach my $check ( @checks )
  {
    my $key = $check->{'key'};
    next unless( exists $suppliedRef->{$key} );

    my $value = getConfigValue( $config, $key );
    die "Invalid array length for $key: final value is missing or not an array ($check->{'type'})\n"
      unless( ref( $value ) eq 'ARRAY' );
    my $length = scalar @$value;

    if( exists $check->{'constraints'}->{'lenmul'} )
    {
      my $rule = $check->{'constraints'}->{'lenmul'};
      $rule =~ m/^([1-9]\d*)\*(.+)$/;
      my $multiplier = $1;
      my $refKey = $2;
      my $refValue = getConfigValue( $config, $refKey );
      die "Invalid array length for $key: reference $refKey is missing or not an array ($check->{'type'})\n"
        unless( ref( $refValue ) eq 'ARRAY' );
      my $expected = $multiplier * scalar @$refValue;
      die "Invalid array length for $key: got $length values, expected $expected from lenmul=$rule ($check->{'type'})\n"
        unless( $length == $expected );
    }

    if( exists $check->{'constraints'}->{'lenmatch'} )
    {
      my $refKey = $check->{'constraints'}->{'lenmatch'};
      my $refValue = getConfigValue( $config, $refKey );
      die "Invalid array length for $key: reference $refKey is missing or not an array ($check->{'type'})\n"
        unless( ref( $refValue ) eq 'ARRAY' );
      my $expected = scalar @$refValue;
      die "Invalid array length for $key: got $length values, expected $expected to match $refKey ($check->{'type'})\n"
        unless( $length == $expected );
    }
  }
}
