/*********************************************************************
*
*  Simple example of using Gmsh to create a simple tet4 mesh.
*  This file defines the physical groups (volumes and surfaces).
*  Note that the actual mesh creation and file exporting must either
*  be done interactively or using command-line arguments.
*
*********************************************************************/

// Include geometry definition.

Include "geometry.geo";

Delete Physicals;


/*********************************************************************
*  Define physical groups.  These are analogous to material blocks and
*  sidesets in Cubit.  I have not yet figured out how to produce the
*  equivalent of a nodeset.
*********************************************************************/

Physical Volume("Upper") = {74, 76};
Physical Volume("Lower") = {78, 80};
Physical Surface("x_neg") = {52, 50};
Physical Surface("x_pos") = {54, 56};
Physical Surface("y_neg") = {38, 40, 36, 34};
Physical Surface("y_pos") = {48, 46, 42, 44};
Physical Surface("z_neg") = {60, 58};
Physical Surface("z_pos") = {64, 62};
Physical Surface("fault") = {66, 68};
