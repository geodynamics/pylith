c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2003  All Rights Reserved
c
c <LicenseText>
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c...  nshape.inc:  A header that defines element type information.
c
c     At present, the largest number of nodes per element is 20 (for a
c     quadratic hexahedron).  The largest number of Gauss points also
c     occurs for this element, with a value of 27 (for 3x3x3
c     integration).  The maximum number of element equations is derived
c     from nenmax, and the number of different element types is
c     presently 62 (this includes all of the infinite element
c     permutations).
c
c     This header should be included before any subroutine arguments,
c     as the parameters will generally be used to dimension the
c     arguments.
c
      integer nenmax,ngaussmax,neemax,netypes,netypesi
      parameter(nenmax=20,ngaussmax=27,neemax=3*nenmax,netypes=62,
     & netypesi=10)
c
c version
c $Id: nshape.inc,v 1.4 2004/07/09 16:25:17 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
