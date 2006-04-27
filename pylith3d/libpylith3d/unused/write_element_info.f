c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams
c  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
c
c  Permission is hereby granted, free of charge, to any person obtaining
c  a copy of this software and associated documentation files (the
c  "Software"), to deal in the Software without restriction, including
c  without limitation the rights to use, copy, modify, merge, publish,
c  distribute, sublicense, and/or sell copies of the Software, and to
c  permit persons to whom the Software is furnished to do so, subject to
c  the following conditions:
c
c  The above copyright notice and this permission notice shall be
c  included in all copies or substantial portions of the Software.
c
c  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
c  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
c  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
c  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
c  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
c  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
c  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
      open(kw,file=ofile,status="old",access="append")
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
c $Id: write_element_info.f,v 1.1 2005/08/05 19:58:07 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
