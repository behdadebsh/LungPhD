
read com;FilesSurfaceLungs/Subject;#define subject
read com;FilesSurfaceLungs/FitType;#define fit type


$fit_files  = $path.'/GeometricModels/FittingLungs/FilesSurfaceLungs';
$data_files = $path.'/Data/'.$study.'/'.$subject.'/'.$protocol;

#$zeroXderiv = $fit_files.'/zeroXderiv.pl';
${tot_itt}=5;#how many iterations to do


fem de param;r;$fit_files.'/3d_fitting';
fem de coor;r;$fit_files.'/versions';
fem de base;r;$fit_files.'/BiCubic_Surface_unit';

if($type == 1){ #use initial placement
   fem de node;r;$data_files.'/Vessel/SurfaceFEMesh/Initial'.$vessel;
   fem de elem;r;$fit_files.'/template_'.$vessel;
   fem export node;$data_files.'/initial' as initial offset 1000;
   fem export elem;$data_files.'/initial' as initial offset_e 1000 offset_n 1000;
}else{ #use pre-fitted and manually edited mesh
   fem de node;r;$data_files.'/Vessel/SurfaceFEMesh/FittedManual'.$vessel;
   fem de elem;r;$fit_files.'/template_'.$vessel;
}

#------------GROUP NODES AND ELEMENTS----------------------------------

#fem export node;$data_files.'/Vessel/SurfaceFEMesh/'.$vessel.'_fitted' as fitted;
#fem export elem;$data_files.'/Vessel/SurfaceFEMesh/'.$vessel.'_fitted' as fitted;
#fem quit;

fem group elem all as all_elements;

fem group node 1,2,3,4 as inlet;
fem group node 31,32,33,34 as outlet1;
fem group node 18,19,20,21 as outlet2;
fem group node 6,8..9 as bif;
fem group node in elem all_elements as all_nodes;

if($type == 5){
   fem group node $node_group as fix_everything;
}


fem update field from geometry;

  $datafile = $data_files.'/Vessel/surface_'.$vessel.'trimmed';
  fem de data;r;$datafile; #read in trimmed surface


fem de xi;c close element all_elements;


for (1..$tot_itt){

	if($type == 1){
		fem de fit;r;$fit_files.'/InitialFit'.$vessel geometry;
    }elsif($type == 2){
		fem de fit;r;$fit_files.'/SecondFit'.$vessel geometry;
	}elsif($type == 3){
	    fem de fit;r;$fit_files.'/ThirdFit'.$vessel geometry;
	}elsif($type == 4){
	    fem de fit;r;$fit_files.'/FourthFit'.$vessel geometry;
	}


	fem def mapping;r;$fit_files.'/Vessel'.$vessel;


    fem fit;
    fem update node fit;
    fem update scale_factor unit; #normalise;

    fem de xi;c close element all_elements;

    print "COMPLETED FITTING ITERATION $_ \n";
    print "data split is $DATA_SPLIT \n";

}


fem de node;w;$data_files.'/Vessel/SurfaceFEMesh/'.$vessel.'_fitted';
fem export node;$data_files.'/Vessel/SurfaceFEMesh/'.$vessel.'_fitted' as fitted;
fem export elem;$data_files.'/Vessel/SurfaceFEMesh/'.$vessel.'_fitted' as fitted;

fem quit;

