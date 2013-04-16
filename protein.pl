#!/usr/bin/perl -w

# pill
foreach $dp ([1.0,0.5]) {
  my $len_cyl = sprintf("%03.1f",$$dp[0]);
  my $rad_endcap = sprintf("%03.1f",$$dp[1]);
  
  my $scriptname = "protein-p-$len_cyl-$rad_endcap.tmp.sh";
  open(SCRIPT, ">$scriptname") or die $!;
  print SCRIPT "#!/bin/sh
  ./protein p $len_cyl $rad_endcap 0.0 0.0";
  close SCRIPT;
  system("sbatch",$scriptname);
}

# cone with length len and base radius rad
foreach $dc ([1.0,1.0]) {
  my $len = sprintf("%03.1f",$$dc[0]);
  my $rad = sprintf("%03.1f",$$dc[1]);
  
  my $scriptname = "protein-c-$len-$rad.tmp.sh";
  
  open SCRIPT, ">$scriptname" or die $!;
  print SCRIPT "#!/bin/sh
  ./protein c $len $rad 0.0 0.0";
  close SCRIPT;
  system("sbatch",$scriptname);
}

# box with three side lengths
foreach $db ([0.5,1.0,0.5]) {
  my $len_x = sprintf("%03.1f",$$db[0]);
  my $len_y = sprintf("%03.1f",$$db[1]);
  my $len_z = sprintf("%03.1f",$$db[2]);
  
  my $scriptname = "protein-b-$len_x-$len_y-$len_z.tmp.sh";
  
  open SCRIPT, ">$scriptname" or die $!;
  print SCRIPT "#!/bin/sh
  ./protein b $len_x $len_y $len_z 0.0";
  close SCRIPT;
  system("sbatch",$scriptname);
  
}

# stadium shape with len of cylinder, rad along x axis, rad along y, and rad along z for end caps
foreach $dst ([0.5,1.0,0.5,0.5]) {
  my $len = sprintf("%03.1f",$$dst[0]);
  my $rad_x = sprintf("%03.1f",$$dst[1]);
  my $rad_y = sprintf("%03.1f",$$dst[2]);
  my $rad_z = sprintf("%03.1f",$$dst[3]);
  
  my $scriptname = "protein-st-$len-$rad_x-$rad_y-$rad_z.tmp.sh";
  
  open SCRIPT, ">$scriptname" or die $!;
  print SCRIPT "#!/bin/sh
  ./protein st $len $rad_x $rad_y $rad_z";
  close SCRIPT;
  system("sbatch",$scriptname);
}

# sphere with radius
foreach $dsp ([1.0]) {
  my $rad = sprintf("%03.1f",$$dsp[0]);
  
  my $scriptname = "protein-sp-$rad.tmp.sh";
  
  open SCRIPT, ">$scriptname" or die $!;
  print SCRIPT "#!/bin/sh
  ./protein b $rad 0.0 0.0 0.0";
  close SCRIPT;
  system("sbatch",$scriptname);
}

# sphere with radius
foreach $de ([0.5,1.0,0.5]) {
  my $len_x = sprintf("%03.1f",$$de[0]);
  my $len_y = sprintf("%03.1f",$$de[1]);
  my $len_z = sprintf("%03.1f",$$de[2]);
  
  
  my $scriptname = "protein-e-$len_x-$len_y-$len_z.tmp.sh";
  
  open SCRIPT, ">$scriptname" or die $!;
  print SCRIPT "#!/bin/sh
  ./protein b $len_x $len_y $len_z 0.0";
  close SCRIPT;
  system("sbatch",$scriptname);
}