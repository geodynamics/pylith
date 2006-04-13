c -*- Fortran -*-
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2005  All Rights Reserved
c
c  Copyright 2005 Rensselaer Polytechnic Institute.
c  All worldwide rights reserved.  A license to use, copy, modify and
c  distribute this software for non-commercial research purposes only
c  is hereby granted, provided that this copyright notice and
c  accompanying disclaimer is not modified or removed from the software.
c
c  DISCLAIMER:  The software is distributed "AS IS" without any express
c  or implied warranty, including but not limited to, any implied
c  warranties of merchantability or fitness for a particular purpose
c  or any warranty of non-infringement of any current or pending patent
c  rights.  The authors of the software make no representations about
c  the suitability of this software for any particular purpose.  The
c  entire risk as to the quality and performance of the software is with
c  the user.  Should the software prove defective, the user assumes the
c  cost of all necessary servicing, repair or correction.  In
c  particular, neither Rensselaer Polytechnic Institute, nor the authors
c  of the software are liable for any indirect, special, consequential,
c  or incidental damages related to the software, to the maximum extent
c  the law permits.
c
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c...  Example include file to define element interpolation information
c     in f77.  This is a temporary template, based on a linear hex
c     element with full quadrature.  It is unclear at present whether
c     it is desirable to put all of this information into an 'include'
c     file, which will require it to be compiled as a separate element
c     type.  Since all of my present routines that deal with
c     interpolation information do not require any knowledge of the
c     element characteristics, it seems a more general approach would
c     be to either tabulate (or read from a tabulated list) this
c     information at runtime.  This would prevent the problem of
c     having to define a new element type corresponding to each possible
c     integration order, for example.
c
c...  Global dimensions, defined as parameters.
c
c       nen:     Number of element nodes
c       ngauss:  Number of element Gauss points
c
      integer nen,ngauss
      parameter(nen=8,ngauss=8)
c
c...  Note that, because of f77 requirements about statement order, it
c     is easiest if the interpolation information is broken up into 2
c     files.  This first file simply defines the global dimensions.
c
c version
c $Id: f77_elem_interp_dims_example.inc,v 1.1 2005/02/28 20:48:44 willic3 Exp $

c Generated automatically by Fortran77Mill on Mon Feb 28 11:48:22 2005

c End of file 
