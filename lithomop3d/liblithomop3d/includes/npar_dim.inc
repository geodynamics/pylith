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
c...  npar_dim.inc:  A header that describes the data type and
c     dimensions of arrays for the npar parameter group.  This
c     header should be included for any TECTON routine that uses
c     the group.  All 'dim' header files must be included before
c     any 'def' header files.
c
c...  arrays containing parameters
c
      integer npar(12)
c
c...  scalar parameters associated with arrays
c
      integer numelv,numat,numpr,numslp,numfn,ipstrs,ipauto,nstatesz,
     & nstatesz0,nvfamilies,numdif,intord
c
c version
c $Id: npar_dim.inc,v 1.5 2005/03/25 01:36:54 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
