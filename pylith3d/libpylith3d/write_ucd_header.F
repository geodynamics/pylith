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
      subroutine write_ucd_header(istatoutc,nstatestot,kucd,iucd)
c
c...subroutine to write the header for a UCD file containing state
c   variables
c
#include "binio.h"
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
c...  external functions
c
      integer nchar
      external nchar
c
c...  local variables
c
      integer i,ibyte,indl,indu
      integer il
      integer iout(3*nstatesmax)
      character nlabels*1024,nunits*1024
      character slabel*100
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
          il=nchar(labels(istatoutc(i)))
          if(istatoutc(i).le.nstatesmax) then
            slabel=labels(istatoutc(i))(1:il)//","//stress
          else if(istatoutc(i).le.2*nstatesmax) then
            slabel=labels(istatoutc(i))(1:il)//","//strain
          else if(istatoutc(i).le.3*nstatesmax) then
            slabel=labels(istatoutc(i))(1:il)//","//strate
          end if
          write(kucd,"(a100)") slabel
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
        write(kucd,BINARY_POS=ibyte) nlabels,nunits
        ibyte=ibyte+2048
        write(kucd,BINARY_POS=ibyte) nstatestot,(iout(i),i=1,nstatestot)
      end if
c
      return
      end
