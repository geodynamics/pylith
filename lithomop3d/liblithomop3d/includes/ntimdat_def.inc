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
      ntdinitp=ntimdat(3)
      lgdefp=ntimdat(4)
      itmaxp=ntimdat(5)
      nittot=ntimdat(6)
      nrftot=ntimdat(7)
      ndtot=ntimdat(8)
      ireform=ntimdat(9)
c
c version
c $Id: ntimdat_def.inc,v 1.4 2004/07/21 15:24:20 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
