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
      subroutine pcginv(alnz,b,d,bres,p,z,dprev,rmin,rmult,
     & gcurr,gprev,ja,nsiter,neq,nnz,ndtot,idout,kto,kw,used)
c
c...  subroutine to solve for the displacement increments d, using a
c     preconditioned conjugate gradient solver.  The type of
c     preconditioning is determined by the parameter nsol.  If nsol=1
c     or nsol=3, the preconditioner is simply the inverse of the
c     stiffness matrix diagonals.  If nsol=2 or nsol=4, a symmetrized
c     Gauss-Seidel preconditioner is used.  This routine also assumes
c     the stiffness matrix is stored in modified sparse row (MSR) format.
c     Note that vector b in this routine actually corresponds to vector
c     bres in the calling routine (iterate), and vice-versa.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer neq,nnz,ndtot,idout,kto,kw
      integer ja(nnz)
      double precision alnz(nnz),b(neq),d(neq),bres(neq),p(neq),z(neq)
      double precision dprev(neq)
      logical used
c
c...  included dimension and type statements
c
      include "nsiter_dim.inc"
      include "rmin_dim.inc"
      include "rmult_dim.inc"
      include "gcurr_dim.inc"
      include "gprev_dim.inc"
c
c...  defined constants
c
      include "nconsts.inc"
      include "rconsts.inc"
c
c...  intrinsic functions
c
      intrinsic min,sqrt,max
c
c...  user-defined functions
c
      double precision dnrm2,ddot
      external dnrm2,ddot
c
c...  local variables
c
cdebug      integer idb
      integer i,j,ip,l
      double precision etmp,ftmp,en,rtz0,ptap,alf,deld,rtz1,bet
      double precision sstol(3),scurr(3),sprev(3),si(3),tmp(3),acc(3)
      double precision cnvnl(3)
      logical div(3),converge
c
c...  included variable definitions
c
      include "nsiter_def.inc"
c
c...  determine convergence criteria based on convergence rate of
c     nonlinear solution
c
cdebug      write(6,*) "Hello from pcginv_f!"
cdebug      write(6,*) "alnz:",(alnz(idb),idb=1,neq)
cdebug      write(6,*) "b:",(b(idb),idb=1,neq)
cdebug      write(6,*) "d:",(d(idb),idb=1,neq)
cdebug      write(6,*) "bres:",(bres(idb),idb=1,neq)
cdebug      write(6,*) "p:",(p(idb),idb=1,neq)
cdebug      write(6,*) "z:",(z(idb),idb=1,neq)
c
      do i=1,3
        cnvnl(i)=gcurr(i)/gprev(i)
	if(cnvnl(i).eq.zero) cnvnl(i)=rmin(i)
        sstol(i)=rmult(i)*cnvnl(i)
        sstol(i)=min(sstol(i),rmin(i))
        si(i)=zero
      end do
      converge=.false.
c
c...  define initial values for first iteration
c
      if(nsol.eq.1.or.nsol.eq.3) then
        ip=1
      else
        ip=2
      end if
      call dcopy(neq,b,ione,bres,ione)
      if(.not.used) then
        call fill(d,zero,neq)
      else
cdebug2        write(6,*) "Test in pcginv!"
        call dcopy(neq,dprev,ione,d,ione)
        etmp=zero
	tmp(1)=dnrm2(neq,d,ione)
	tmp(2)=dnrm2(neq,bres,ione)
        tmp(3)=zero
        do i=1,neq
          en=d(i)*bres(i)
          tmp(3)=tmp(3)+en*en
          bres(i)=bres(i)-d(i)*alnz(i)
          do j=ja(i),ja(i+1)-1
            bres(i)=bres(i)-d(ja(j))*alnz(j)
          end do
          en=d(i)*bres(i)
          etmp=etmp+en*en
        end do
	ftmp=dnrm2(neq,bres,ione)
        etmp=sqrt(etmp)
        tmp(3)=sqrt(tmp(3))
        do i=1,3
          si(i)=tmp(i)
          scurr(i)=tmp(i)
          sprev(i)=tmp(i)
        end do
	acc(1)=one
	acc(2)=ftmp/tmp(2)
	acc(3)=etmp/tmp(3)
	l=0
	converge=.true.
	if(idout.gt.1) write(kw,800) l
	write(kto,800) l
	go to 10
      end if
      if(ip.eq.1) then
        call dcopy(neq,bres,ione,z,ione)
	call dtbsv("u","n","n",neq,izero,alnz,ione,z,ione)
        call dcopy(neq,z,ione,p,ione)
	rtz0=ddot(neq,bres,ione,z,ione)
      else if(ip.eq.2) then
        call gspre(alnz,bres,z,rtz0,ja,neq,nnz)
        call dcopy(neq,z,ione,p,ione)
      end if
c
c...  loop over number of iterations
c
      do l=1,maxcg
        ndtot=ndtot+1
        ptap=zero
c
c...  compute alfa for this iteration
c
        call dcopy(neq,p,ione,z,ione)
	call dtbmv("u","n","n",neq,izero,alnz,ione,z,ione)
        do i=1,neq
          do j=ja(i),ja(i+1)-1
            z(i)=z(i)+p(ja(j))*alnz(j)
          end do
        end do
	ptap=ddot(neq,p,ione,z,ione)
        alf=rtz0
        if(ptap.ne.zero) alf=rtz0/ptap
c
c...  compute displacements and residuals for this iteration
c
        call fill(tmp,zero,ithree)
	ftmp=dnrm2(neq,bres,ione)
        do i=1,neq
          deld=alf*p(i)
          en=deld*bres(i)
          tmp(1)=tmp(1)+deld*deld
          tmp(3)=tmp(3)+en*en
          d(i)=d(i)+deld
        end do
	call daxpy(neq,-alf,z,ione,bres,ione)
	if(ip.eq.1) then
          call dcopy(neq,bres,ione,z,ione)
	  call dtbsv("u","n","n",neq,izero,alnz,ione,z,ione)
	  rtz1=ddot(neq,bres,ione,z,ione)
	end if
	tmp(2)=dnrm2(neq,bres,ione)
        tmp(1)=sqrt(tmp(1))
        tmp(3)=sqrt(tmp(3))
        if(l.eq.1.and.(.not.used)) then
          si(1)=tmp(1)
          si(2)=max(ftmp,tmp(2))
          si(3)=tmp(3)
          do i=1,3
            scurr(i)=si(i)
            sprev(i)=si(i)
          end do
        else if(l.eq.2) then
          do i=1,3
            if(tmp(i).gt.si(i)) then
              si(i)=tmp(i)
              scurr(i)=tmp(i)
              sprev(i)=tmp(i)
            end if
          end do
        end if
        do i=1,3
          if(si(i).eq.zero) si(i)=tmp(i)
          acc(i)=zero
          if(si(i).ne.zero) acc(i)=tmp(i)/si(i)
c*          div(i)=tmp(i).gt.scurr(i).and.scurr(i).gt.sprev(i)
          div(i)=.false.
          sprev(i)=scurr(i)
          scurr(i)=tmp(i)
        end do
c
c...  return if solution has converged, otherwise compute new values for
c     p(i).
c
        if(div(1).and.div(2).and.div(3)) then
          converge=.true.
          if(idout.gt.1) write(kw,800) l
          write(kto,800) l
        else if(acc(1).le.sstol(1).and.acc(2).le.sstol(2).and.acc(3).
     &   le.sstol(3)) then
          converge=.true.
          if(idout.gt.1) write(kw,810) l
          write(kto,810) l
        end if
        if(converge) go to 10
c
c...  compute new z-vector using symmetrized Gauss-Seidel preconditioner
c
        if(ip.eq.2) call gspre(alnz,bres,z,rtz1,ja,neq,nnz)
        bet=rtz1
        if(rtz0.ne.zero) bet=rtz1/rtz0
        rtz0=rtz1
	call dscal(neq,bet,p,ione)
	call daxpy(neq,one,z,ione,p,ione)
      end do
c
c...  issue a warning message if solution didn't converge
c
      if(idout.gt.1) write(kw,820) maxcg
      write(kto,820) maxcg
10    continue
      if(nsol.gt.2) call dcopy(neq,d,ione,p,ione)
      call dcopy(neq,d,ione,z,ione)
      if(idout.gt.1) write(kw,900) (cnvnl(i),i=1,3),(scurr(i),i=1,3),
     & (si(i),i=1,3),(acc(i),i=1,3),(sstol(i),i=1,3)
      write(kto,900) (cnvnl(i),i=1,3),(scurr(i),i=1,3),
     & (si(i),i=1,3),(acc(i),i=1,3),(sstol(i),i=1,3)
cdebug      write(6,*) "From end of pcginv_f!"
cdebug      write(6,*) "alnz:",(alnz(idb),idb=1,neq)
cdebug      write(6,*) "b:",(b(idb),idb=1,neq)
cdebug      write(6,*) "d:",(d(idb),idb=1,neq)
cdebug      write(6,*) "bres:",(bres(idb),idb=1,neq)
cdebug      write(6,*) "p:",(p(idb),idb=1,neq)
cdebug      write(6,*) "z:",(z(idb),idb=1,neq)
800   format(/,
     & '     WARNING!  Apparent divergence in subiteration #',i5,'!',//)
810   format(/,
     & '     Solution converged in',i5,' subiterations',//)
820   format(/,
     & '     WARNING!  Solution failed to converge in',i5,
     & ' subiterations!',//)
900   format(29x,'Displacement',9x,'Force',14x,'Energy',/,3x,
     & '     Nonlinear         ',/,3x,
     & '      Convergence Rate ',1pe15.8,3x,1pe15.8,3x,1pe15.8,/,3x,
     & '     Final Norm        ',1pe15.8,3x,1pe15.8,3x,1pe15.8,/,3x,
     & '     Initial Norm      ',1pe15.8,3x,1pe15.8,3x,1pe15.8,/,3x,
     & '     Final/Initial     ',1pe15.8,3x,1pe15.8,3x,1pe15.8,/,3x,
     & '     Required Ratio    ',1pe15.8,3x,1pe15.8,3x,1pe15.8,//)
      return
      end
c
c version
c $Id: pcginv.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
