read com;FilesSurfaceLungs/Subject; #has to be in the same directory
$data_files = $path.'/Data/'.$study.'/'.$subject.'/'.$protocol;
$fit_files  = $path.'/GeometricModels/FittingLungs/FilesSurfaceLungs';

$exnode = $data_files.'/Vessel/SurfaceFEMesh/'.$vessel.'_fitted';

gfx read node $exnode; #has group name 'fitted'
gfx read elem $exnode;

gfx mod g_e fitted general clear circle_discretization 6 default_coordinate coordinates element_discretization "12*12*12" native_discretization none;
gfx mod g_e fitted lines coordinate coordinates select_on material green selected_material default_selected;
gfx mod g_e fitted node_points glyph sphere general size "4*4*4" material green;


########################################3
# Data points
#---------------

$group = 'surface_'.$vessel.'trimmed_nobase';
$datafile = $data_files.'/Vessel/'.$group;
if(-e $datafile.'.exdata'){
   gfx read data $datafile;
   gfx mod g_e $group general clear;
   gfx mod g_e $group data_points glyph point size "1*1*1" no_select material default;
   $group = 'surface_'.$vessel.'trimmed_base';
   $datafile = $data_files.'/Vessel/'.$group;
   gfx read data $datafile;
   gfx mod g_e $group general clear;
   gfx mod g_e $group data_points glyph sphere size "2*2*2" no_select material gold;
   print "Reading $group \n";
}else{
   $group = 'surface_'.$vessel;
   $datafile = $data_files.'/Vessel/'.$group.'trimmed';
   gfx read data $datafile;
   gfx mod g_e $group general clear;
   gfx mod g_e $group data_points glyph point size "1*1*1" no_select material default;
   print "Reading $group \n";

}

gfx edit scene;


open comfile;$fit_files.'/UpdateFittedGeometryVessel.com';
open comfile;$fit_files.'/WriteFittedEditVessel.com';
open comfile;$fit_files.'/edit_derivatives.com' exec 1;

