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
c...  rconsts_par.inc:  A set of commonly-used constants.
c
      double precision zero,one,two,three,four,five,six,seven,eight,nine
      double precision half,third,fourth,fifth,sixth,eighth,root3,root3i
      parameter(zero   = 0.0d0,
     &          one    = 1.0d0,
     &          two    = 2.0d0,
     &          three  = 3.0d0,
     &          four   = 4.0d0,
     &          five   = 5.0d0,
     &          six    = 6.0d0,
     &          seven  = 7.0d0,
     &          eight  = 8.0d0,
     &          nine   = 9.0d0)
c
      parameter(half   = 0.5d0,
     &          third  = one/three,
     &          fourth = 0.25d0,
     &          fifth  = 0.2d0,
     &          sixth  = one/six,
     &          eighth = 0.125d0,
     &          root3  = 1.732050807568877293527446341505d0,
     &          root3i = 0.577350269189625764509148780502d0,
     &          pi     = 3.141592653589793238462643383276d0)
c
c version
c $Id: rconsts.inc,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
