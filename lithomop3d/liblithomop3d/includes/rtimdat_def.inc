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
c...  rtimdat_def.inc:  A header that associates scalars with elements
c     of an array.  In general, the array and scalars should have
c     already been defined with a 'xxx_dim.inc' file.  This include
c     file actually contains executable (assignment) statements, so
c     it should be included after all type, data, etc. statements.
c
c...  associate local names with array elements
c
      deltp=rtimdat(1)
      alfap=rtimdat(2)
      tpois=rtimdat(3)
      tyoungs=rtimdat(4)
c
c version
c $Id: rtimdat_def.inc,v 1.3 2005/02/24 00:32:40 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
