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
c...  npar_def.inc:  A header that associates scalars with elements
c     of an array.  In general, the array and scalars should have
c     already been defined with a 'xxx_dim.inc' file.  This include
c     file actually contains executable (assignment) statements, so
c     it should be included after all type, data, etc. statements.
c
c...  associate local names with array elements
c
      numelv=npar(1)
      numat=npar(2)
      numpr=npar(3)
      numslp=npar(4)
      numfn=npar(5)
      ipstrs=npar(6)
      ipauto=npar(7)
      nstatesz=npar(8)
      nstatesz0=npar(9)
      nvfamilies=npar(10)
      numdif=npar(11)
      intord=npar(12)
c
c version
c $Id: npar_def.inc,v 1.5 2005/03/25 01:36:54 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
