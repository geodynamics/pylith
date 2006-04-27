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
      subroutine write_sparse_info(neq,nnz,nmin,nmax,wavg,idout,kw,
     & ofile)
c
c...subroutine to write global control parameters
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer neq,nnz,nmin,nmax,idout,kw
      double precision wavg
      character ofile*(*)
c
cdebug      write(6,*) "Hello from write_sparse_info_f!"
c
      if(idout.gt.0) then
        open(kw,file=ofile,status="old",access="append")
        write(kw,3000) neq,nnz-1,nmin,nmax,wavg
        close(kw)
      end if
c
      return
c
 3000 format(///,
     &' e q u a t i o n   s y s t e m   d a t a               ',   //,
     & 5x, 
     &' number of equations  . . . . . . . . . . . . . (neq) =',i8,//,
     & 5x,
     &' number of nonzero terms in stiffness matrix. . (nnz) =',i8,//,
     & 5x,
     &' minimum number of nonzero terms per row. . . .(nmin) =',i8,//,
     & 5x, 
     &' maximum number of nonzero terms per row. . . .(nmax) =',i8,//,
     & 5x, 
     &' average number of nonzero terms per row. . . .(wavg) =',f8.2)
c
      end
c
c version
c $Id: write_sparse_info.f,v 1.1 2005/08/05 19:58:08 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
