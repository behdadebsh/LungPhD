#!/usr/bin/perl

# uses data point locations to place and scale an initial mesh

use strict;

my $study   = $ARGV[0];
my $subject = $ARGV[1];
my $vessel    = $ARGV[2];
my $protocol = $ARGV[3];
my $Dir     = $ARGV[4];
my $protocol_choice=$ARGV[5];

my $DirMesh = $Dir.'/Data/'.$study.'/'.$subject.'/'.$protocol.'/Vessel';

# directory containing surface data file
my $datafile = $DirMesh.'/surface_'.$vessel.'trimmed.ipdata';
my $opfile   = $DirMesh.'/SurfaceFEMesh/Initial'.$vessel.'.exnode';
my $nodefile   = $Dir.'/GeometricModels/FittingLungs/FilesSurfaceLungs/template_'.$vessel.'.exnode';

print "open data file $datafile \n";
open IPFILE, "<$datafile" or die "\033[31mError: Can't open data file\033[0m ";

my $x = 0;
my $y = 1;
my $z = 2;

my $max_x = -1e6;
my $max_y = -1e6;
my $max_z = -1e6;
my $min_x = 1e6;
my $min_y = 1e6;
my $min_z = 1e6;
my ($mean_x, $mean_y, $mean_z, $Ndata);

# get range of data and centre of mass. use to locate and scale initial mesh
my $line = <IPFILE>;
while ($line){
    if($line =~ /\d+\s*(\S*)\s*(\S*)\s*(\S*)\s*1.0/){
	my $xx = $1;
	my $yy = $2;
	my $zz = $3;
	$mean_x = $mean_x + $xx;
	$mean_y = $mean_y + $yy;
	$mean_z = $mean_z + $zz;
	if($xx < $min_x) {$min_x = $xx;}
	if($yy < $min_y) {$min_y = $yy;}
	if($zz < $min_z) {$min_z = $zz;}
	if($xx > $max_x) {$max_x = $xx;}
	if($yy > $max_y) {$max_y = $yy;}
	if($zz > $max_z) {$max_z = $zz;}
	$Ndata++;
    }
    $line = <IPFILE>;
}
close IPFILE;
$mean_x = $mean_x / $Ndata;
$mean_y = $mean_y / $Ndata;
$mean_z = $mean_z / $Ndata;

my $range_x = $max_x - $min_x;
my $range_y = $max_y - $min_y;
my $range_z = $max_z - $min_z;

###########################################################################
# CALCULATE THE CENTRE OF MASS AND RANGES FOR THE GENERIC MESH

print "open initial node file $nodefile \n";
open EXNODE, "<$nodefile" or die "\033[31mError: Can't open node file\033[0m ";
my $NumberofNodes=0;
my $max_x = -1e6;
my $max_y = -1e6;
my $max_z = -1e6;
my $min_x = 1e6;
my $min_y = 1e6;
my $min_z = 1e6;
my @xyz;
my $nversn;
my ($mesh_x, $mesh_y, $mesh_z);
my ($x_ind, $y_ind, $z_ind, $next_node);

my $line = <EXNODE>;
while ($line) {

    if ($line =~ /Versions=(\d)/){$nversn=$1;}
    if ($line =~ /Node:\s+(\d+)/) {
	$next_node = $1;
	$NumberofNodes=$NumberofNodes+1;
    
	for (my $i=1;$i<=3;$i++){
	    $line = <EXNODE>;
	    if($line =~/\s+(\S+\.\S+e\S+)\s+S*/){
		$xyz[$i] = $1; #get x, y, or z
	    }else{
		if($line =~/\s+(\S+\.\S+E\S+)\s+S*/){
		    $xyz[$i] = $1; #get x, y, or z
		}
	    }
	    for (my $nv = 2; $nv <=$nversn; $nv++){
		$line = <EXNODE>; 
	    }
	}
	
	$mesh_x=$mesh_x+$xyz[1];
	$mesh_y=$mesh_y+$xyz[2];
	$mesh_z=$mesh_z+$xyz[3];

	if($xyz[1] < $min_x) {$min_x = $xyz[1];}
	if($xyz[2] < $min_y) {$min_y = $xyz[2];}
	if($xyz[3] < $min_z) {$min_z = $xyz[3];}
	if($xyz[1] > $max_x) {$max_x = $xyz[1];}
	if($xyz[2] > $max_y) {$max_y = $xyz[2];}
	if($xyz[3] > $max_z) {$max_z = $xyz[3];}
    }
    $line = <EXNODE>;    
}

my $mesh_range_x = $max_x - $min_x;
my $mesh_range_y = $max_y - $min_y;
my $mesh_range_z = $max_z - $min_z;

close EXNODE;

my (@scale_by, @move_by);
$scale_by[1] = $range_x/$mesh_range_x; # * 0.85;
$scale_by[2] = $range_y/$mesh_range_y; # * 0.85;
$scale_by[3] = $range_z/$mesh_range_z; # * 0.85;
$mesh_x=$mesh_x/$NumberofNodes * $scale_by[1];
$mesh_y=$mesh_y/$NumberofNodes * $scale_by[2];
$mesh_z=$mesh_z/$NumberofNodes * $scale_by[3];

$move_by[1] = $mean_x - $mesh_x;
$move_by[2] = $mean_y - $mesh_y;
$move_by[3] = ($mean_z - $mesh_z); # + 20;

#print "scale by $scale_by[1], $scale_by[2], $scale_by[3]\n";
#print "move by $move_by[1], $move_by[2], $move_by[3]\n";

###########################################################################
# SCALE AND SHIFT THE INITIAL MESH

open EXNODE, "<$nodefile" or die "\033[31mError: Can't open node file\033[0m ";
open OPNODE, ">$opfile" or die "\033[31mError: Can't open node file\033[0m ";

my $line = <EXNODE>;
my $deriv = 3;

# the following assumes exnode file format as written out by cmgui
while ($line) {
#    if ($line =~ /Group name:/) { print OPNODE " Group name: initial\n"};
    if ($line =~ /Versions=(\d+)/){$nversn=$1;}
    if ($line =~ /\s*Node:\s+(\d+)/) {
	my $node = $1;
	print OPNODE "$line";
	for (1..3){
	    for my $nv (1..$nversn){
		my $line = <EXNODE>;
		chomp $line;
		my @coord = split(/ +/,$line);
		$coord[1] = $coord[1] * $scale_by[$_] + $move_by[$_];
		printf OPNODE "  %12.6e %12.6e %12.6e %12.6e \n", $coord[1],$coord[2],$coord[3],$coord[4];
	    } #nv
	}
    }else{ #write out everything else verbatim
	print OPNODE "$line";
    }
    $line = <EXNODE>;    
}
### End Loop over nodes

close EXNODE;
close OPNODE;

#########################################

print "Template mesh scaled and moved; written to $opfile\n";
system "perl FilesSurfaceLungs/ConvertExToIpnode.pl $opfile";

