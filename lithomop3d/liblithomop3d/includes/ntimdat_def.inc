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
c...  ntimdat_def.inc:  A header that associates scalars with elements
c     of an array.  In general, the array and scalars should have
c     already been defined with a 'xxx_dim.inc' file.  This include
c     file actually contains executable (assignment) statements, so
c     it should be included after all type, data, etc. statements.
c
c...  associate local names with array elements
c
      nstep=ntimdat(1)
      maxitp=ntimdat(2)
      ntdinit=ntimdat(3)
      lgdefp=ntimdat(4)
      ibbarp=ntimdat(5)
      itmaxp=ntimdat(6)
      nittot=ntimdat(7)
      nrftot=ntimdat(8)
      ndtot=ntimdat(9)
      ireform=ntimdat(10)
c
c version
c $Id: ntimdat_def.inc,v 1.2 2004/06/23 20:42:16 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
