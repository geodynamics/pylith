# Create geometry for vertical strike-slip fault.
algebraic3d

solid cube = orthobrick (0, 0, 0; 1, 1, 1);
solid fault = plane (0.5, 0.5, 0.5; -1., 0., 0.) -bc=1;

solid mat1 = cube and fault;
solid mat2 = cube and not fault;

tlo mat1 -col=[0,0,1];
tlo mat2 -col=[1,0,0];
