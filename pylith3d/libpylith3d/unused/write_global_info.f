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
      subroutine write_global_info(title,idout,idsk,numnp,
     & icode,idebug,kw,kp,ofile,pfile)
c
c...subroutine to write global control parameters
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer idout,idsk,numnp,icode,idebug,kw,kp
      character title*(*),ofile*(*),pfile*(*)
c
      if(idout.gt.izero) then
        open(kw,file=ofile,status="new",access="append")
        write(kw,500)
        write(kw,1000) title,idout,idsk,ngem,numnp,nsd,ndof,icode,idebug
        close(kw)
      end if
      if(idsk.eq.ione) then
        open(kp,file=pfile,status="new",access="append")
        write(kp,2000) title
        write(kp,3000) ngem,numnp,nsd,ndof,nstr
        close(kp)
      else if(idsk.eq.itwo) then
        open(kp,file=pfile,status="new",access="append",
     &   form="unformatted")
        write(kp) title
        write(kp) ngem,numnp,nsd,ndof,nstr
        close(kp)
      end if
c
      return
c
 500  format(//,7x,
     x'******************************************************',/,7x,
     x'*                                                    *',/,7x,
     x'*                PYLITH 0.72 OUTPUT                  *',/,7x,
     x'*                                                    *',/,7x,
     x'*        Copyright 2005 by Charles A. Williams.      *',/,7x,
     x'*          Rensselaer Polytechnic Institute.         *',/,7x,
     x'*                All rights reserved.                *',/,7x,
     x'*                                                    *',/,7x,
     x'******************************************************',//)
1000  format(1x,a70,///,
     x' c o n t r o l   i n f o r m a t i o n                 ',  //,5x,
     x' ascii file output option . . . . . . . . . . (idout) =',i5,//,
     x 5x,
     x'    eq.0,  no ascii output file                        ',   /,5x,
     x'    eq.1,  echo input                                  ',   /,5x,
     x'    eq.2,  full output                                 ',  //,5x,
     x' fedsk.dat output option. . . . . . . . . . . (idsk ) =',i5,//,
     x 5x,
     x'    eq.0,  no output                                   ',   /,5x,
     x'    eq.1,  ascii output                                ',   /,5x,
     x'    eq.2,  binary output                               ',  //,5x,
     x' geometry type  . . . . . . . . . . . . . . . (ngem ) =',i5,//,
     x 5x,
     x'    eq.0,  axisymmetric                                ',  / ,5x,
     x'    eq.1,  plane strain                                ',  / ,5x,
     x'    eq.2,  plane stress                                ',  / ,5x,
     x'    eq.3,  out-of-plane                                ',  / ,5x,
     x'    eq.4,  three-dimensional                           ',  //,5x,
     x' number of nodal points . . . . . . . . . . . (numnp) =',i7,//,
     x 5x,
     x' number of space dimensions . . . . . . . . . (nsd  ) =',i5,//,
     x 5x,
     x' number of degrees of freedom per node  . . . (ndof ) =',i5,//,
     x 5x,
     x' analysis code  . . . . . . . . . . . . . . . (icode) =',i5,/,5x,
     x'    eq.0,  set up storage, echo input data             ',  / ,5x,
     x'    eq.1,  factor stiffness, print diagonals           ',  / ,5x,
     x'    eq.2,  elastic solution only                       ',  / ,5x,
     x'    eq.3,  full viscoelastic solution                  ',  //,5x,
     x' level of debugging output . . . . . . . . . (idebug) =',i5,/,5x,
     x'    eq.0,  no special debug output                     ',  / ,5x,
     x'    eq.1,  print local and global force vectors        ',  //,5x,
     x/)
2000  format(a70)
3000  format(16i7)
c
      end
c
c version
c $Id: write_global_info.f,v 1.1 2005/08/05 19:58:08 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
