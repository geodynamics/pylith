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
      subroutine printv(a,b,id,idx,neq,numnp,iout,idout,kw)
c
c...program to print one or two vectors of length neq, along with
c   data on correlation between equation number and node
c
c     if iout = 0, print array d after factorzation a=(u)t*d*u
c        iout = 1, print global load vectors for debugging
c        iout = 2, print diagonals of stiffness matrix
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "ndimens.inc"
c
c...  subroutine arguments
c
      integer neq,numnp,iout,idout,kw
      integer id(ndof,numnp),idx(ndof,numnp)
      double precision a(neq),b(neq)
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer npage,n,i,j,iflag,idof,node
c
cdebug      write(6,*) "Hello from printv_f!"
c
      if(idout.lt.2) return
      npage=50
      do n=1,neq
c
c...find node and dof associated with equation number
c
        do i=1,ndof
          do j=1,numnp
            iflag=0
            if(id(i,j).eq.n) goto 50
            if(idx(i,j).eq.n) then
              iflag=1
              goto 50
            end if
          end do
        end do
   50   idof=i
        node=j
c
c...print out the appropriate vectors
c
        if(iout.eq.0) then
          if((n.eq.1).or.(mod(n,npage).eq.0)) write(kw,1500)
          if(iflag.eq.0) write(kw,2000) n,node,idof,a(n)
          if(iflag.eq.1) write(kw,3000) n,node,idof,a(n)
        else if(iout.eq.1) then
          if((n.eq.1).or.(mod(n,npage).eq.0)) write(kw,1000)
          if(iflag.eq.0) write(kw,2000) n,node,idof,a(n),b(n),
     &     b(n)-a(n)
          if(iflag.eq.1) write(kw,3000) n,node,idof,a(n),b(n),
     &     b(n)-a(n)
        else if(iout.eq.2) then
          if((n.eq.1).or.(mod(n,npage).eq.0)) write(kw,1700)
          if(iflag.eq.0) write(kw,2000) n,node,idof,a(n)
          if(iflag.eq.1) write(kw,3000) n,node,idof,a(n)
        end if
      end do
      return
 1000 format(//' global force vector follows'//1x,'  eq# ',' node',
     & '    dof  ',5x,' current b',7x,'applied b',5x,'difference'/)
 1500 format(//' diagonal elements of stiffness matrix:'///
     & 2x,' neq ',' node ','  dof ',7x, '    diag'//)
 1700 format(//' diagonal elements of stiffness matrix:'///
     & 2x,' neq ',' node ','  dof ',7x, '    diag'//)
 2000 format(1x,3(i7,1x ),1x,3(1x,1pe15.8))
 3000 format(1x,3(i7,'x'),1x,3(1x,1pe15.8))
      end
c
c version
c $Id: printv.f,v 1.2 2004/07/05 19:42:08 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
