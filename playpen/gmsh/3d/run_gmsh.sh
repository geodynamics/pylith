#!/bin/bash
#
# This just creates a 3D mesh (optimized) and the output is just the
# default (gmsh) format.
gmsh mesh_tet4.geo -3 -optimize
