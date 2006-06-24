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
c...  pntdef.inc:  A header that assigns pointer values to scalars from
c     arrays, as a means of avoiding common blocks.  This header should
c     be included in the main calling routines of TECTON, after all
c     dimensioning/typing statements.  This version uses equivalence
c     statements so that changes in the scalar pointers are reflected in
c     the array values.
c
      equivalence(ipien,npepoint(1))
      equivalence(iplm,npepoint(2))
      equivalence(iplmx,npepoint(3))
      equivalence(iplmf,npepoint(4))
      equivalence(ipmat,npepoint(5))
      equivalence(ipinf,npepoint(6))
      equivalence(ipprop,npepoint(7))
      equivalence(ipnflt,npfpoint(1))
      equivalence(ipflt,npfpoint(2))
      equivalence(ipdflt,npfpoint(3))
      equivalence(ipdftn,npfpoint(4))
      equivalence(ipnslp,npfpoint(5))
      equivalence(ipdiff,npfpoint(6))
      equivalence(ipdslp,npfpoint(7))
      equivalence(ipdx,npfpoint(8))
      equivalence(ipdelx,npfpoint(9))
      equivalence(ipdhis,npfpoint(10))
      equivalence(ippslp,npfpoint(11))
      equivalence(ipid,npgpoint(1))
      equivalence(ipidx,npgpoint(2))
      equivalence(ipbdat,npgpoint(3))
      equivalence(ipx,npgpoint(4))
      equivalence(ipd,npgpoint(5))
      equivalence(ipdeld,npgpoint(6))
      equivalence(ipbond,npgpoint(7))
      equivalence(ipskew,npgpoint(8))
      equivalence(ipwink,npgpoint(9))
      equivalence(ipwdat,npgpoint(10))
      equivalence(iphist,npgpoint(11))
      equivalence(ipdcur,npitpoint(1))
      equivalence(ipdxcur,npitpoint(2))
      equivalence(ipdeps,npitpoint(3))
      equivalence(ipdbeta,npitpoint(4))
      equivalence(ipdbetb,npitpoint(5))
      equivalence(ipscur,npitpoint(6))
      equivalence(ipkout,npopoint(1))
      equivalence(ipsout,npopoint(2))
      equivalence(ipkflt,npopoint(3))
      equivalence(ipkx,npopoint(4))
      equivalence(ipelno,npppoint(1))
      equivalence(ipside,npppoint(2))
      equivalence(ipphis,npppoint(3))
      equivalence(ippres,npppoint(4))
      equivalence(ippdir,npppoint(5))
      equivalence(ipbeta,npspoint(1))
      equivalence(ipbetb,npspoint(2))
      equivalence(ipdmat,npspoint(3))
      equivalence(ipstn,npspoint(4))
      equivalence(ipeps,npspoint(5))
      equivalence(ipst0,npspoint(6))
      equivalence(ips,npspoint(7))
      equivalence(ipstmp,npspoint(8))
      equivalence(iplnz,npsppoint(1))
      equivalence(ipja,npsppoint(2))
      equivalence(ippcg,npsppoint(3))
      equivalence(ipzcg,npsppoint(4))
      equivalence(ipb,npsppoint(5))
      equivalence(ipbtot,npsppoint(6))
      equivalence(ipbres,npsppoint(7))
      equivalence(ipgvc1,npsppoint(8))
      equivalence(ipgvc2,npsppoint(9))
      equivalence(ippvc,npsppoint(10))
      equivalence(ipdprv,npsppoint(11))
      equivalence(ipeste,npsppoint(12))
      equivalence(ipestf,npsppoint(13))
      equivalence(ipdiag,npsppoint(14))
      equivalence(ipdiagl,npsppoint(15))
      equivalence(ipneex,npsppoint(16))
      equivalence(iplmm,npsppoint(17))
      equivalence(iplmmi,npsppoint(18))
      equivalence(ipmstp,nptpoint(1))
      equivalence(ipalfa,nptpoint(2))
      equivalence(ipdelt,nptpoint(3))
      equivalence(ipmaxi,nptpoint(4))
      equivalence(ipmaxc,nptpoint(5))
      equivalence(iplgdf,nptpoint(6))
      equivalence(ipbbar,nptpoint(7))
      equivalence(iputol,nptpoint(8))
      equivalence(ipftol,nptpoint(9))
      equivalence(ipetol,nptpoint(10))
      equivalence(ipitmx,nptpoint(11))
      equivalence(ipiprint,nptpoint(12))
      equivalence(iptflt,npvpoint(1))
      equivalence(ipwinx,npvpoint(2))
      equivalence(ipxdat,npvpoint(3))
c
c version
c $Id: pntdef.inc,v 1.1 2004/07/05 18:35:49 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
