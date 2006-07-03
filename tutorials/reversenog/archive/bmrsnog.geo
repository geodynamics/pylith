# Attempt to create geometry for SCEC BM5.
algebraic3d

solid left = plane (0., 12., -12.; -1., 0., 0.) -bc=1;
solid right = plane (24., 12., -12.; 1., 0., 0.) -bc=2;
solid front = plane (12., 0., -12.; 0., -1., 0.) -bc=3;
solid back = plane (12., 24., -12.; 0., 1., 0.) -bc=4;
solid bottom = plane (12., 12., -24.; 0., 0., -1.) -bc=5;
solid top = plane (12., 12., 0.; 0., 0., 1.) -bc=6;
solid diag = plane ( 16., 12., -12.; -1., 0., -1.) -bc=7;
solid horiz1 = plane (12., 12., -12.; 0., 0., -1.) -bc=8;
solid horiz2 = plane (12., 12., -16.; 0., 0., -1.) -bc=9;

solid mat1 = left and front and back and bottom and not horiz2 and not diag and right -maxh = 1.5;
solid mat2 = left and front and back and horiz2 and not horiz1 and not diag -maxh = 1.5;
solid mat3 = left and front and back and horiz1 and top and not diag -maxh = 1.5;
solid mat4 = right and front and back and bottom and not horiz2 and diag -maxh = 1.5;
solid mat5 = right and front and back and horiz2 and not horiz1 and diag -maxh = 1.5;
solid mat6 = right and front and back and horiz1 and top and diag -maxh = 1.5;

tlo mat1 -col=[0,0,1];
tlo mat2 -col=[0,1,0];
tlo mat3 -col=[1,0,0];
tlo mat4 -col=[0,1,1];
tlo mat5 -col=[1,1,0];
tlo mat6 -col=[1,0,1];
