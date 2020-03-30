#!/usr/bin/perl

# Perl script to control some of the fitting of a Main Pulmonary Artery. 

$cmiss_root=$ENV{CMISS_EXECUTABLE}; #where is your Cmiss?
$cmgui_root=$ENV{CMGUI_2_6_2}; #Where is your cmgui?
$root=$ENV{LUNG_ROOT};
$fiji=$ENV{FIJI};#Where is fiji or imagej executable 

my $study = 'Human_PE_Study_HRC';
my $subject = "";
my $vessel = 'MPA';
my $protocol = 'TLC';
#my $fissure = "";

for my $N (0..@ARGV-1){
    if($ARGV[$N] =~ /-h/) {
	print "Options for FITTING: \n";
	print " -s study name (HLA/HPE, default HLA)\n"; 
	print " -n name of subject (use hXXX shortcut for HLA only)\n"; 
	print " -p protocol eg (FRC/TLC, default empty)\n"; 
	die;
    }
    if($ARGV[$N] =~ /-s/) {
	$study = $ARGV[$N+1];
	if($study =~/HLA/){$study = 'Human_Lung_Atlas';}
	if($study =~/HPE/){$study = 'Human_PE_Study_HRC';}
	if($study =~/CF/){$study = 'Paediatric_CF';}
        if($study =~/IPF/){$study = 'Human_IPF';}
    }
    if($ARGV[$N] =~ /-n/) {$subject = $ARGV[$N+1];}
    if($ARGV[$N] =~ /-p/) {$protocol = $ARGV[$N+1];}
    if($ARGV[$N] =~ /-v/) {$vessel = $ARGV[$N+1];}
	#if($ARGV[$N] =~ /-f/) {$fissure = $ARGV[$N+1];}
}

my $protocol_choice;
if($protocol eq '.'){
    print "What volume would you like to initialise your fit at (FRC/TLC)\n";
    $protocol_choice = <STDIN>;
}else{
    $protocol_choice =$protocol;
}

my $subject_name = `python $root/UniversalFiles/read_props.py $root $study $subject Name`;
chomp $subject_name;
print "Subject name is: $subject_name \n";

# print out the subject information to a file
system "perl FilesSurfaceLungs/SubjectData.pl $study $subject_name $vessel $protocol $root $protocol_choice";

#################################################################
# step 1 : scale and translate the template mesh
#################################################################
# 
my $opdir   = $root.'/Data/'.$study.'/'.$subject_name.'/'.$protocol.'/Vessel/SurfaceFEMesh';
$file=$opdir.'/'.$vessel.'_fitted.exnode';
#
if(-e $file){
	print "No need for initialisation \n";

}else{
	print "No $file: initialising \n";
	system "mkdir $opdir"; # will fail if exists
	my $plfile =  $root.'/GeometricModels/FittingLungs/FilesSurfaceLungs/InitialPositionVessel.pl';
	system "perl $plfile $study $subject_name $vessel $protocol $root $protocol_choice";



#################################################################
# step 2 : first geometry fitting
#################################################################
#
	# First time. $type = 1, for fixed fissure nodes only. 
	my $comfile = $root.'/GeometricModels/FittingLungs/FilesSurfaceLungs/FitType.com';
	open(OPCOM, ">$comfile") or die "Can't open comfile\n";
	print OPCOM " \$type = 1;\n"; 
	print OPCOM " \$fix = 0;\n";
	print OPCOM " \$node_group = 0;\n";
	close OPCOM;

	my $comfile = $root.'/GeometricModels/FittingLungs/FilesSurfaceLungs/FitMesh_Vessel';
	system "$cmiss_root $comfile";

}


#################################################################
# step 3 : visualise the fitted mesh and edit nodes
#################################################################
	my $drawfile = $root.'/GeometricModels/FittingLungs/FilesSurfaceLungs/DrawVesselEdit.com';
	system "$cmgui_root $drawfile & ";
	sleep(5);
#
my $type = 1; #default
while ($type > 0){
    print "Update the fitted node locations using UpdateFittedGeometry.com\n";
    print "Edit the fitted placement of nodes\n";
    print "Use WriteFittedEdit.com to write out the edited mesh (needed for the next step)\n";
    print "Enter RETURN when ready to re-fit:\n";
    my $nothing = <STDIN>;


	print "Enter fitting type: (2) STEP 2 - INITIAL LANDMARKS; (3) STEP 3 - OPTIONAL (FIXING DEFINED NODE COORDINATES AND/OR DERIVATIVES); (4) RESET; or (0) exit: \n";

    my $option = <STDIN>;
    chomp $option;
    if($option == 0){exit;} #all done! 

   #if($option == 5) {
    if($option == 5) {
      $node_group = 0;
      print "Enter node numbers to fully fix: ";
        $node_group = <STDIN>;
        chomp $node_group;}


    my $comfile = 'FilesSurfaceLungs/FitType.com';
    open(OPCOM, ">$comfile") or die "Can't open comfile\n";
    print OPCOM " \$type = $option;\n"; 
    print OPCOM " \$fix = 0;\n";
    print OPCOM " \$node_group = \'$node_group\';\n";
    close OPCOM;

#################################################################
# step 4 : subsequent geometry fitting and viewing
#################################################################
#

    system "$cmiss_root FilesSurfaceLungs/FitMesh_Vessel";
    
}
