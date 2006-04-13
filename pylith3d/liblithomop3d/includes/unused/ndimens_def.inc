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
c...  ndimens_def.inc:  A header that associates scalars with elements
c     of an array.  In general, the array and scalars should have
c     already been defined with a 'xxx_dim.inc' file.  This include
c     file actually contains executable (assignment) statements, so
c     it should be included after all type, data, etc. statements.
c
c...  associate local names with array elements
c
      nsd=ndimens(1)
      ndof=ndimens(2)
      nstr=ndimens(3)
      nen=ndimens(4)
      ngem=ndimens(5)
      nddmat=ndimens(6)
      nskdim=ndimens(7)
      nsdim=ndimens(8)
      npdim=ndimens(9)
      npdir=ndimens(10)
c
c version
c $Id: ndimens_def.inc,v 1.1 2004/06/18 14:44:10 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
