c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2004  All Rights Reserved
c
c <LicenseText>
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c...  nforce_def.inc:  A header that associates scalars with elements
c     of an array.  In general, the array and scalars should have
c     already been defined with a 'xxx_dim.inc' file.  This include
c     file actually contains executable (assignment) statements, so
c     it should be included after all type, data, etc. statements.
c
c...  associate local names with array elements
c
      nextflag=nforce(1)
      ntractflag=nforce(2)
      ngravflag=nforce(3)
      nconcflag=nforce(4)
      nprestrflag=nforce(5)
      nwinkflag=nforce(6)
      nwinkxflag=nforce(7)
      nprevdflag=nforce(8)
c
c version
c $Id: nforce_def.inc,v 1.2 2005/05/03 18:42:49 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
