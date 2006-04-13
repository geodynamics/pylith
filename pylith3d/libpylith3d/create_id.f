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
