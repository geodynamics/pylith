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
      subroutine prints(stn,eps,deps,beta,dbeta,betb,dbetb,nstr,
     & ngauss,numel,nstep,idout,idsk,kw,kp,ivisc,iplas)
c
c...program to print stress and strain
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nstr,ngauss,numel,nstep,idout,idsk,kw,kp,ivisc,iplas
      double precision stn(nstr,ngauss,numel),eps(nstr,ngauss,numel)
      double precision deps(nstr,ngauss,numel),beta(nstr,ngauss,numel)
      double precision dbeta(nstr,ngauss,numel),betb(nstr,ngauss,numel)
      double precision dbetb(nstr,ngauss,numel)
c
c...  included dimension and type statements
c
      include "labels_dim.inc"
c
c...  intrinsic functions
c
      intrinsic mod
c
c...  local variables
c
      integer m50,n,i,l
c
c...  included variable definitions
c
      include "labels_def.inc"
c
cdebug      write(6,*) "Hello from prints_f!"
c
      m50=50/ngauss
      if(idout.gt.1.or.idsk.eq.0) then
        do n=1,numel
          if((n.eq.1.or.mod(n,m50).eq.0).and.idout.gt.1) then
            write(kw,2000) (labels(i),i=1,nstr)
            write(kw,*) ' '
          end if
          do l=1,ngauss
            if(idsk.eq.0) write(kp,1500) (stn(i,l,n),i=1,nstr)
            if(idout.gt.1) write(kw,3000) n,(stn(i,l,n),i=1,nstr)
          end do
          if(idsk.eq.0) then
            do l=1,ngauss
              write(kp,1500) (eps(i,l,n),i=1,nstr)
            end do
            if(nstep.gt.0) then
              if(ivisc.eq.1) then
                do l=1,ngauss
                  write(kp,1500) (beta(i,l,n),i=1,nstr)
                end do
              end if
              if(iplas.eq.1) then
                do l=1,ngauss
                  write(kp,1500) (betb(i,l,n),i=1,nstr)
                end do
              end if
              do l=1,ngauss
                write(kp,1500) (deps(i,l,n),i=1,nstr)
              end do
              if(ivisc.eq.1) then
                do l=1,ngauss
                  write(kp,1500) (dbeta(i,l,n),i=1,nstr)
                end do
              end if
              if(iplas.eq.1) then
                do l=1,ngauss
                  write(kp,1500) (dbetb(i,l,n),i=1,nstr)
                end do
              end if
            end if
          end if
        end do
      end if
      if(idsk.eq.1) then
        write(kp) stn
        write(kp) eps
        if(nstep.gt.0) then
          if(ivisc.eq.1) write(kp) beta
          if(iplas.eq.1) write(kp) betb
          write(kp) deps
          if(ivisc.eq.1) write(kp) dbeta
          if(iplas.eq.1) write(kp) dbetb
        end if
      end if
      return
2000  format(///1x,' s t r e s s              '///
     &' elem #',5x,6(a4,10x),a4)
3000  format(1x,i7,1x,6(1pe14.5))
1500  format(6e16.7)
      end
c
c version
c $Id: prints.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
