c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2005  All Rights Reserved
c
c  Copyright 2005 Rensselaer Polytechnic Institute.
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
      subroutine write_ucd_header(istatoutc,nstatestot,kucd,iucd)
c
c...subroutine to write the header for a UCD file containing state
c   variables
c
      include "implicit.inc"
c
c...  parameter definitions
c
      include "materials.inc"
      include "nconsts.inc"
c
c...  subroutine arguments
c
      integer nstatestot,kucd,iucd
      integer istatoutc(nstatestot)
c
c... included dimension and type statements
c
      include "labels_dim.inc"
c
c...  local constants
c
      character stress*6,strain*4,strate*9
      data stress,strain,strate/"Pascal","None","1/seconds"/
c
c...  local variables
c
      integer i,ibyte,indl,indu
      integer iout(3*nstatesmax)
      character nlabels*1024,nunits*1024
c
c...  included variable definitions
c
      include "labels_def.inc"
c
      do i=1,3*nstatesmax
        iout(i)=ione
      end do
c
      if(iucd.eq.ione) then
        write(kucd,"(100i5)") nstatestot,(iout(i),i=1,nstatestot)
c
        do i=1,nstatestot
          write(kucd,"(a11)") labels(istatoutc(i))
        end do
      else if(iucd.eq.itwo) then
        indl=ione
        indu=ione
        ibyte=ione
        do i=1,nstatestot
          nlabels(indl:indl+11)=labels(istatoutc(i))//"."
          indl=indl+12
          if(istatoutc(i).le.nstatesmax) then
            nunits(indu:indu+6)=stress//"."
            indu=indu+7
          else if(istatoutc(i).le.2*nstatesmax) then
            nunits(indu:indu+4)=strain//"."
            indu=indu+5
          else
            nunits(indu:indu+9)=strate//"."
            indu=indu+10
          end if
        end do
        nlabels(1024:1024)="0"
        nunits(1024:1024)="0"
c
c...  The g95 compiler supports stream i/o, while the gfortran compiler
c     does not.  At present, it is assumed that all f95 compilers other
c     than gfortran support this feature.  Undoubtedly, exceptions will
c     be found and the list below will need to be expanded.
c
#if defined FORTRAN_TYPE_F77
        write(kucd,rec=ibyte) nlabels,nunits
        ibyte=ibyte+2048
        write(kucd,rec=ibyte) nstatestot,(iout(i),i=1,nstatestot)
#elif defined FORTRAN_NAME_GFORTRAN
#warning "Binary UCD output not supported with gfortran!"
#else
        write(kucd,pos=ibyte) nlabels,nunits
        ibyte=ibyte+2048
        write(kucd,pos=ibyte) nstatestot,(iout(i),i=1,nstatestot)
#endif
      end if
c
      return
      end
