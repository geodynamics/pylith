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
c...  nshape.inc:  A header that defines the maximum number of nodes
c     and Gauss points per element.
c     At present, the largest number of nodes per element is 20 (for a
c     quadratic hexahedron).  The largest number of Gauss points also
c     occurs for this element, with a value of 27 (for 3x3x3
c     integration).
c     This header should be included before any subroutine arguments,
c     as the parameters will generally be used to dimension the
c     arguments.
c
      integer nenmax,ngaussmax
      parameter(nenmax=20,ngaussmax=27)
c
c version
c $Id: nshape.inc,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
