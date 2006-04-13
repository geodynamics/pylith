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
c...  ndimens.inc:  A header that defines most of the dimensions
c     specific to 3D problems.
c
c     This header should be included before any subroutine arguments,
c     as the parameters will generally be used to dimension the
c     arguments.
c
      integer nsd,ndof,nstr,ngem,nddmat,nskdim,nsdim,npdim,npdir
      parameter(
     & nsd=3,
     & ndof=3,
     & nstr=6,
     & ngem=4,
     & nddmat=21,
     & nskdim=2,
     & nsdim=5,
     & npdim=4,
     & npdir=2)
c
c version
c $Id: ndimens.inc,v 1.2 2004/06/18 14:25:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
