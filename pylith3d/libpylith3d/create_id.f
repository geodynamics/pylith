c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
c
c  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
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
      subroutine create_id(id,idx,ibond,nslip,idslp,numslp,numnp,
     & numsn,neq)
c
c...  subroutine to set up the global equation numbers.  These are
c     stored in id and idx.  Equations are eliminated for kinematic
c     BC, and new equations are added for slippery nodes.
c
c     Input information:
c
c         numslp:            Number of slippery node entries.
c         numnp:             Total number of nodes.
c         numsn:             Number of slippery nodes.
c         ibond(ndof,numnp): Integer array specifying boundary condition
c                            type for each node (zero indicates no
c                            boundary condition).
c         nslip(ndof+2,numslp):  Integer array specifying slippery
c                                nodes.
c               nslip(1,m)  = element number for slippery node entry m.
c               nslip(2,m)  = node number for slippery node entry m.
c               nslip(2+ndof,m) = slip indicator for each slippery node
c                                 dof:
c                                 -1 = negative side of fault
c                                  1 = positive side of fault
c                                  0 = no slip for this dof
c
c     Output information:
c
c         neq:             Number of global equations.
c         idslp(numsn):    Node number corresponding to each slippery
c                          node.
c         id(ndof,numnp):  Array containing the equation number for each
c                          node and degree of freedom.
c         idx(ndof,numnp): Array containing the equation numbers for
c                          each slippery node and degree of freedom.
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
      integer numslp,numnp,numsn,neq,ns
      integer id(ndof,numnp),idx(ndof,numnp),ibond(ndof,numnp)
      integer nslip(nsdim,numslp),idslp(numsn)
c
c...  included dimension and type statements
c
c
c...  local constants
c
c
c...  intrinsic functions
c
c
c...  local variables
c
      integer i,j,n,ii,nn,imode,iihist,iitype
c
c.... establish equation numbers and create id array
c
      neq=izero
      do n=1,numnp
        do i=1,ndof
          id(i,n)=izero
          imode=ibond(i,n)
          iihist=imode/10
          iitype=imode-10*iihist
          if((iitype.eq.izero).or.(iitype.eq.ithree)) then
            neq=neq+1
            id(i,n)=neq
          end if
        end do
      end do
c
c...  if there are slippery nodes, create idx array and add in new
c     equations
c
      call ifill(idx,izero,ndof*numnp)
      if(numslp.eq.izero) return
      ns=izero
      do n=1,numnp
        do j=1,numslp
          if(n.eq.nslip(2,j)) then
            ns=ns+1
            idslp(ns)=n
            do i=1,ndof
              if((nslip(2+i,j).ne.izero).and.(id(i,n).ne.izero)) then
                idx(i,n)=id(i,n)
                do ii=i,ndof
                  if(id(ii,n).ne.izero) id(ii,n)=id(ii,n)+1
                end do
                do nn=n+1,numnp
                  do ii=1,ndof
                    if(id(ii,nn).ne.izero) id(ii,nn)=id(ii,nn)+1
                  end do
                end do
                neq=neq+1
              end if
            end do
            goto 100
          end if
        end do
  100   continue
      end do
c
      return
      end
c
c version
c $Id: create_id.f,v 1.2 2005/04/14 00:57:14 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
