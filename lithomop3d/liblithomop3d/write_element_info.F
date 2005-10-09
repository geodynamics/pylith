c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2005  All Rights Reserved
c
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
      subroutine write_element_info(numelv,nen,ngauss,ietypev,intord,
     & ipstrs,ipauto,tpois,tyoungs,kw,idout,ofile)
c
c...subroutine to write element and prestress parameters
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer numelv,nen,ngauss,ietypev,intord,ipstrs,ipauto,kw,idout
      double precision tpois,tyoungs
      character ofile*(*)
c
c...  included dimension and type statements
c
      include "elmlbl_dim.inc"
c
c...  local variables
c
      character intorder(3)*17
      data intorder/"             Full",
     &              "          Reduced",
     &              "Selective (B-bar)"/
c
c...  included variable definitions
c
      include "elmlbl_def.inc"
c
c...  echo input to output file
c
      if(idout.eq.izero) return
c
c...  The current setup uses autoconf to determine whether the fortran
c     compiler supports access="append", position="append", or neither.
c     If neither is supported an error is returned.
c
#if defined FORTRAN_ACCESS_APPEND
      open(kw,file=ofile,status="old",access="append")
#elif defined FORTRAN_POSITION_APPEND
      open(kw,file=ofile,status="old",position="append")
#else
#error "Cannot determine method for appending to Fortran files!"
#endif
      write(kw,700) elmlbl,numelv,nen,ngauss,ietypev,intorder(intord),
     & ipstrs,ipauto,tpois,tyoungs
      close(kw)
c
700   format(1x,///,
     &' e l e m e n t    s y s t e m   d a t a',///,5x,
     &' element type:  ',a40,//,5x,
     &' number of volume elements. . . . . . . . . (numelv) =',i7,//,5x,
     &' number of volume element nodes . . . . . .    (nen) =',i7,//,5x,
     &' number of volume element Gauss points. . . (ngauss) =',i7,//,5x,
     &' volume element type. . . . . . . . . . . .(ietypev) =',i7,//,5x,
     &' integration order . . . . . . . . . . . = ',a17,//,5x,
     &' prestress option. . . . . . . . . . . . . .(ipstrs) =',i5,/ ,5x,
     &'    eq.0, prestresses are read from the input file    ', / ,5x,
     &'    eq.1, gravitational prestresses automatically     ', / ,5x,
     &'          computed                                    ', / ,5x,
     &' prestress auto-computation option . . . . .(ipauto) =',i5,/ ,5x,
     &'    eq.0, computation uses assigned elastic properties', / ,5x,
     &'    eq.1, properties listed below are used for        ', / ,5x,
     &'          auto-computation                            ', / ,5x,
     &' poissons ratio for auto-computation . . . . (tpois) =',1pe15.8,
     & /,5x,'    only used for ipstrs=1 and ipauto=1',/,5x,
     &' youngs modulus for auto-computation . . . (tyoungs) =',1pe15.8,
     & /,5x,'    only used for ipstrs=1 and ipauto=1',/)
      return
      end
c
c version
c $Id: write_element_info.F,v 1.1 2005/08/05 19:58:06 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
