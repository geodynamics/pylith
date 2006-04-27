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
      subroutine bsum(bextern,btraction,bgravity,bconcforce,bintern,
     & bresid,nextflag,ntractflag,ngravflag,nconcflag,neq)
c
c...subroutine to sum all contributions to the external force vector,
c   then subtract the internal force vector from the external vector
c   to obtain the residual.
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  subroutine arguments
c
      integer nextflag,ntractflag,ngravflag,nconcflag,neq
      double precision bextern(nextflag*neq),btraction(ntractflag*neq)
      double precision bgravity(ngravflag*neq),bconcforce(nconcflag*neq)
      double precision bintern(neq),bresid(neq)
c
c...  local constants
c
      double precision alpha
      data alpha/-1.0d0/
c
cdebug      integer idb
c
cdebug      write(6,*) "Hello from bsum_f!"
c
      if(nextflag.ne.izero) then
        call fill(bextern,zero,neq)
        if(ntractflag.ne.izero) call daxpy(neq,one,btraction,ione,
     &   bextern,ione)
        if(ngravflag.ne.izero) call daxpy(neq,one,bgravity,ione,
     &   bextern,ione)
        if(nconcflag.ne.izero) call daxpy(neq,one,bconcforce,ione,
     &   bextern,ione)
        call dcopy(neq,bextern,ione,bresid,ione)
cdebug        write(6,*) "bextern:"
cdebug        write(6,*) (bextern(idb),idb=1,200)
cdebug        if(ntractflag.ne.izero) then
cdebug          write(6,*) "btraction:"
cdebug          write(6,*) (btraction(idb),idb=1,200)
cdebug        end if
cdebug        if(ngravflag.ne.izero) then
cdebug          write(6,*) "bgravity:"
cdebug          write(6,*) (bgravity(idb),idb=1,200)
cdebug        end if
cdebug        if(nconcflag.ne.izero) then
cdebug          write(6,*) "bconcforce:"
cdebug          write(6,*) (bconcforce(idb),idb=1,200)
cdebug        end if
      else
        call fill(bresid,zero,neq)
      end if
      call daxpy(neq,alpha,bintern,ione,bresid,ione)
cdebug      write(6,*) "bintern:"
cdebug      write(6,*) (bintern(idb),idb=1,200)
cdebug      write(6,*) "bresid:"
cdebug      write(6,*) (bresid(idb),idb=1,200)
      return
      end
c
c version
c $Id: bsum.f,v 1.3 2005/04/08 00:30:36 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
