c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2004  All Rights Reserved
c
c  Copyright 2004 Rensselaer Polytechnic Institute.
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
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      subroutine write_element_info(numel,ngauss,ipstrs,tpois,nppts,kw,
     & idout,ofile)
c
c...subroutine to write element and prestress parameters
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer numel,ngauss,ipstrs,nppts,kw,idout
      double precision tpois
      character ofile*(*)
c
c...  included dimension and type statements
c
      include "elmlbl_dim.inc"
c
c...  included variable definitions
c
      include "elmlbl_def.inc"
c
c...  echo input to output file
c
      if(idout.gt.0) then
        open(kw,file=ofile,status="old",access="append")
        write(kw,700) elmlbl,numel,ngauss,ipstrs,tpois,nppts
        close(kw)
      end if
c
700   format(1x,///,
     &' e l e m e n t    s y s t e m   d a t a',///,5x,
     &' element type:  ',a40,//,5x,
     &' number of elements . . . . . . . . . . . . .(numel) =',i7,//,5x,
     &' number of gauss points per element . . . . (ngauss) =',i7,//,5x,
     &' prestress option. . . . . . . . . . . . . .(ipstrs) =',i5,/ ,5x,
     &'    eq.0, prestresses are read from the input file    ', / ,5x,
     &'    eq.1, prestresses computed from elastic solution  ', / ,5x,
     &'          assuming near-incompressibility             ', / ,5x,
     &' poissons ratio for prestresses. . . . . . . (tpois) =',1pe15.8,
     & /,5x,'    only used for ipstrs=1',/ ,5x,
     &' number of prestress integration points. . . (nppts) =',i5,/)
      return
      end
c
c version
c $Id: write_element_info.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
