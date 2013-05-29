/*********************************************************************
*
*  Simple example of using Gmsh to create a simple mesh analogous to
*  3d/hex8.  This file defines the geometry.
*
*********************************************************************/

Delete Model;

// Define characteristic size that is assigned to all entities.

cellSize = 1000.0;

/*********************************************************************
*  Define a set of points that will be used to define lines.
*********************************************************************/

Point( 1) = {-3000.0, -3000.0, -4000.0, cellSize};
Point( 2) = {-3000.0, -3000.0, -2000.0, cellSize};
Point( 3) = {-3000.0, -3000.0,     0.0, cellSize};
Point( 4) = {    0.0, -3000.0, -4000.0, cellSize};
Point( 5) = {    0.0, -3000.0, -2000.0, cellSize};
Point( 6) = {    0.0, -3000.0,     0.0, cellSize};
Point( 7) = { 3000.0, -3000.0, -4000.0, cellSize};
Point( 8) = { 3000.0, -3000.0, -2000.0, cellSize};
Point( 9) = { 3000.0, -3000.0,     0.0, cellSize};
Point(10) = {-3000.0,  3000.0, -4000.0, cellSize};
Point(11) = {-3000.0,  3000.0, -2000.0, cellSize};
Point(12) = {-3000.0,  3000.0,     0.0, cellSize};
Point(13) = {    0.0,  3000.0, -4000.0, cellSize};
Point(14) = {    0.0,  3000.0, -2000.0, cellSize};
Point(15) = {    0.0,  3000.0,     0.0, cellSize};
Point(16) = { 3000.0,  3000.0, -4000.0, cellSize};
Point(17) = { 3000.0,  3000.0, -2000.0, cellSize};
Point(18) = { 3000.0,  3000.0,     0.0, cellSize};

/*********************************************************************
*  Connect the points to form lines.
*********************************************************************/

// Z-lines.

Line( 1) = { 1, 2};
Line( 2) = { 2, 3};
Line( 3) = { 4, 5};
Line( 4) = { 5, 6};
Line( 5) = { 7, 8};
Line( 6) = { 8, 9};
Line( 7) = {10,11};
Line( 8) = {11,12};
Line( 9) = {13,14};
Line(10) = {14,15};
Line(11) = {16,17};
Line(12) = {17,18};

// Y-lines.

Line(13) = { 1,10};
Line(14) = { 4,13};
Line(15) = { 7,16};
Line(16) = { 2,11};
Line(17) = { 5,14};
Line(18) = { 8,17};
Line(19) = { 3,12};
Line(20) = { 6,15};
Line(21) = { 9,18};

// X-lines.

Line(22) = { 1, 4};
Line(23) = { 4, 7};
Line(24) = {10,13};
Line(25) = {13,16};
Line(26) = { 2, 5};
Line(27) = { 5, 8};
Line(28) = {11,14};
Line(29) = {14,17};
Line(30) = { 3, 6};
Line(31) = { 6, 9};
Line(32) = {12,15};
Line(33) = {15,18};

/*********************************************************************
*  Form line loops from the lines.
*********************************************************************/

Line Loop(1) = { 22,  3,-26, -1};
// Line Loop(2) = { 23,  5,-27, -3};
// Line Loop(3) = {-24,  7, 20, -9};

Plane Surface(34) = {1};
Line Loop(35) = {23, 5, -27, -3};
Plane Surface(36) = {35};
Line Loop(37) = {26, 4, -30, -2};
Plane Surface(38) = {37};
Line Loop(39) = {27, 6, -31, -4};
Plane Surface(40) = {39};
Line Loop(41) = {24, 9, -28, -7};
Plane Surface(42) = {41};
Line Loop(43) = {25, 11, -29, -9};
Plane Surface(44) = {43};
Line Loop(45) = {28, 10, -32, -8};
Plane Surface(46) = {45};
Line Loop(47) = {29, 12, -33, -10};
Plane Surface(48) = {47};
Line Loop(49) = {13, 7, -16, -1};
Plane Surface(50) = {49};
Line Loop(51) = {16, 8, -19, -2};
Plane Surface(52) = {51};
Line Loop(53) = {18, 12, -21, -6};
Plane Surface(54) = {53};
Line Loop(55) = {15, 11, -18, -5};
Plane Surface(56) = {55};
Line Loop(57) = {13, 24, -14, -22};
Plane Surface(58) = {57};
Line Loop(59) = {14, 25, -15, -23};
Plane Surface(60) = {59};
Line Loop(61) = {21, -33, -20, 31};
Plane Surface(62) = {61};
Line Loop(63) = {20, -32, -19, 30};
Plane Surface(64) = {63};
Line Loop(65) = {17, 10, -20, -4};
Plane Surface(66) = {65};
Line Loop(67) = {14, 9, -17, -3};
Plane Surface(68) = {67};
Line Loop(69) = {26, 17, -28, -16};
Plane Surface(70) = {69};
Line Loop(71) = {27, 18, -29, -17};
Plane Surface(72) = {71};
Surface Loop(73) = {64, 46, 52, 38, 70, 66};
Volume(74) = {73};
Surface Loop(75) = {62, 54, 48, 40, 66, 72};
Volume(76) = {75};
Surface Loop(77) = {70, 50, 58, 42, 34, 68};
Volume(78) = {77};
Surface Loop(79) = {72, 68, 60, 44, 56, 36};
Volume(80) = {79};
Physical Volume("Upper") = {74, 76};
Physical Volume("Lower") = {78, 80};
Physical Surface("x_neg") = {52, 50};
Physical Surface("x_pos") = {54, 56};
Physical Surface("y_neg") = {38, 40, 36, 34};
Physical Surface("y_pos") = {48, 46, 42, 44};
Physical Surface("z_neg") = {60, 58};
Physical Surface("z_pos") = {64, 62};
Physical Surface("fault") = {66, 68};
