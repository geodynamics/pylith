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
c...  pntdim.inc:  A header that describes the data type and dimensions
c     of arrays containing pointers for the program TECTON.  This
c     header should be included for the main calling routines of the
c     code.
c
c...  arrays containing pointers
c
      integer npepoint(7),npfpoint(11),npgpoint(11),npitpoint(6)
      integer npopoint(4),npppoint(5),npspoint(8),npsppoint(18)
      integer nptpoint(12),npvpoint(3),ippnt(8)
c
c...  scalar pointers associated with arrays
c
      integer ipien,iplm,iplmx,iplmf,ipmat,ipinf,ipprop
      integer ipnflt,ipflt,ipdflt,ipdftn,ipnslp,ipdiff,ipdslp,ipdx,
     & ipdelx,ipdhis,ippslp
      integer ipid,ipidx,ipbdat,ipx,ipd,ipdeld,ipbond,ipskew,ipwink,
     & ipwdat,iphist
      integer ipdcur,ipdxcur,ipdeps,ipdbeta,ipdbetb,ipscur
      integer ipkout,ipsout,ipkflt,ipkx
      integer ipelno,ipside,ipphis,ippres,ippdir
      integer ipbeta,ipbetb,ipdmat,ipstn,ipeps,ipst0,ips,ipstmp
      integer iplnz,ipja,ippcg,ipzcg,ipb,ipbtot,ipbres,ipgvc1,ipgvc2,
     & ippvc,ipdprv,ipeste,ipestf,ipdiag,ipdiagl,ipneex,iplmm,iplmmi
      integer ipmstp,ipalfa,ipdelt,ipmaxi,ipmaxc,iplgdf,ipbbar,iputol,
     & ipftol,ipetol,ipitmx,ipiprint
      integer iptflt,ipwinx,ipxdat
c
c version
c $Id: pntdim.inc,v 1.1 2004/07/05 18:35:49 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
