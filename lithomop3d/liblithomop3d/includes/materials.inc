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
c...  materials.inc:  A header that defines the maximum number
c     of material models and state variables.
c     At present, the maximum number of material models is 20 and the
c     maximum number of state variables is 4 (stress, strain, viscous
c     strain, and plastic strain.
c    
c     This header should be included before any subroutine arguments,
c     as the parameters will generally be used to dimension the
c     arguments.
c
      integer nmatmodmax,nstatesmax
      parameter(nmatmodmax=20,nstatesmax=4)
c
c version
c $Id: materials.inc,v 1.1 2004/06/18 13:38:12 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
