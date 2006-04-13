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
c...  nsysdat_def.inc:  A header that associates scalars with elements
c     of an array.  In general, the array and scalars should have
c     already been defined with a 'xxx_dim.inc' file.  This include
c     file actually contains executable (assignment) statements, so
c     it should be included after all type, data, etc. statements.
c
c...  associate local names with array elements
c
      numnp=nsysdat(1)
      neq=nsysdat(2)
      nnz=nsysdat(3)
      numrot=nsysdat(4)
      nprestr=nsysdat(5)
      numsn=nsysdat(6)
      numflt=nsysdat(7)
      npropsz=nsysdat(8)
      nwink=nsysdat(9)
      nwinkx=nsysdat(10)
      iskopt=nsysdat(11)
c
c version
c $Id: nsysdat_def.inc,v 1.2 2004/07/05 18:31:14 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
