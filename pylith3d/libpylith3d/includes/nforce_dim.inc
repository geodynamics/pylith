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
c...  nforce_dim.inc:  A header that describes the data type and
c     dimensions of arrays for the nforce parameter group.  This
c     header should be included for any TECTON routine that uses
c     the group.  All 'dim' header files must be included before
c     any 'def' header files.
c
c...  arrays containing parameters
c
      integer nforce(8)
c
c...  scalar parameters associated with arrays
c
      integer nextflag,ntractflag,ngravflag,nconcflag,nprestrflag,
     & nwinkflag,nwinkxflag,nprevdflag
c
c version
c $Id: nforce_dim.inc,v 1.2 2005/05/03 18:42:49 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
