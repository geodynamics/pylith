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
c...  rtimdat_dim.inc:  A header that describes the data type and
c     dimensions of arrays for the rtimdat parameter group.  This
c     header should be included for any TECTON routine that uses
c     the group.  All 'dim' header files must be included before
c     any 'def' header files.
c
c...  arrays containing parameters
c
      real*8 rtimdat(3)
c
c...  scalar parameters associated with arrays
c
      real*8 deltp,alfap,autofac
c
c version
c $Id: rtimdat_dim.inc,v 1.2 2005/01/17 20:43:38 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
