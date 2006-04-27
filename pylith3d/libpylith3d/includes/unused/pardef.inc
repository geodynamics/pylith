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
c...  pardef.inc:  A header that assigns parameter values to scalars
c     from arrays, as a means of avoiding common blocks.  This header
c     should be included in the main calling routines of TECTON, after
c     all dimensioning/typing statements.  This version uses
c     equivalence statements so that changes in the scalar parameters
c     are reflected in the array values.
c
      equivalence(icode,ncodat(1))
      equivalence(idebug,ncodat(2))
      equivalence(izero,nconsts(1))
      equivalence(ione,nconsts(2))
      equivalence(itwo,nconsts(3))
      equivalence(ithree,nconsts(4))
      equivalence(nsd,ndimens(1))
      equivalence(ndof,ndimens(2))
      equivalence(nstr,ndimens(3))
      equivalence(nen,ndimens(4))
      equivalence(ngem,ndimens(5))
      equivalence(ngemd,ndimens(6))
      equivalence(nddmat,ndimens(7))
      equivalence(nskdim,ndimens(8))
      equivalence(nsdim,ndimens(9))
      equivalence(npdim,ndimens(10))
      equivalence(npdir,ndimens(11))
      equivalence(numel,npar(1))
      equivalence(numat,npar(2))
      equivalence(ipoint,npar(3))
      equivalence(numpr,npar(4))
      equivalence(nsout,npar(5))
      equivalence(numslp,npar(6))
      equivalence(numfn,npar(7))
      equivalence(ngstn,npar(8))
      equivalence(ndmat,npar(9))
      equivalence(ipstrs,npar(10))
      equivalence(nppts,npar(11))
      equivalence(icontr,nprint(1))
      equivalence(idout,nprint(2))
      equivalence(idsk,nprint(3))
      equivalence(numbc,nrdimens(1))
      equivalence(nsforc,nrdimens(2))
      equivalence(nsol,nsiter(1))
      equivalence(ipcg,nsiter(2))
      equivalence(nstore,nsiter(3))
      equivalence(nstorf,nsiter(4))
      equivalence(maxcg,nsiter(5))
      equivalence(iword,nsysdat(1))
      equivalence(mtot,nsysdat(2))
      equivalence(numnp,nsysdat(3))
      equivalence(neq,nsysdat(4))
      equivalence(nnz,nsysdat(5))
      equivalence(nee,nsysdat(6))
      equivalence(nec,nsysdat(7))
      equivalence(numrot,nsysdat(8))
      equivalence(nprestr,nsysdat(9))
      equivalence(numsn,nsysdat(10))
      equivalence(numflt,nsysdat(11))
      equivalence(nprop,nsysdat(12))
      equivalence(nwink,nsysdat(13))
      equivalence(nwinkx,nsysdat(14))
      equivalence(nemax,nsysdat(15))
      equivalence(iskopt,nsysdat(16))
      equivalence(nstep,ntimdat(1))
      equivalence(maxitp,ntimdat(2))
      equivalence(maxitcp,ntimdat(3))
      equivalence(lgdefp,ntimdat(4))
      equivalence(ibbarp,ntimdat(5))
      equivalence(itmaxp,ntimdat(6))
      equivalence(nittot,ntimdat(7))
      equivalence(nrftot,ntimdat(8))
      equivalence(ndtot,ntimdat(9))
      equivalence(ivisc,ntimdat(10))
      equivalence(iplas,ntimdat(11))
      equivalence(kti,nunits(1))
      equivalence(kto,nunits(2))
      equivalence(kr,nunits(3))
      equivalence(kw,nunits(4))
      equivalence(kp,nunits(5))
      equivalence(kss,nunits(6))
      equivalence(ksf,nunits(7))
      equivalence(ncycle,nvisdat(1))
      equivalence(nintg,nvisdat(2))
      equivalence(numout,nvisdat(3))
      equivalence(nfout,nvisdat(4))
      equivalence(nxout,nvisdat(5))
      equivalence(lastep,nvisdat(6))
      equivalence(nhist,nvisdat(7))
      equivalence(reform,ltimdat(1))
      equivalence(visc,ltimdat(2))
      equivalence(plas,ltimdat(3))
      equivalence(mhist,ltimdat(4))
      equivalence(zero,rconsts(1))
      equivalence(one,rconsts(2))
      equivalence(two,rconsts(3))
      equivalence(three,rconsts(4))
      equivalence(third,rconsts(5))
      equivalence(root3,rconsts(6))
      equivalence(g,rconsts(7))
      equivalence(stol,rgiter(1))
      equivalence(dtol,rgiter(2))
      equivalence(epert,rgiter(3))
      equivalence(rmultd,rmult(1))
      equivalence(rmultf,rmult(2))
      equivalence(rmulte,rmult(3))
      equivalence(rmind,rmin(1))
      equivalence(rminf,rmin(2))
      equivalence(rmine,rmin(3))
      equivalence(deltp,rtimdat(1))
      equivalence(alfap,rtimdat(2))
      equivalence(tpois,rtimdat(3))
      equivalence(cscalx,cscal(1))
      equivalence(cscaly,cscal(2))
      equivalence(cscalz,cscal(3))
      equivalence(wscalx,wscal(1))
      equivalence(wscaly,wscal(2))
      equivalence(wscalz,wscal(3))
      equivalence(wxscalx,wxscal(1))
      equivalence(wxscaly,wxscal(2))
      equivalence(wxscalz,wxscal(3))
      equivalence(gravx,grav(1))
      equivalence(gravy,grav(2))
      equivalence(gravz,grav(3))
      equivalence(prscalxx,prscal(1))
      equivalence(prscalyy,prscal(2))
      equivalence(prscalzz,prscal(3))
      equivalence(prscalxy,prscal(4))
      equivalence(prscalxz,prscal(5))
      equivalence(prscalyz,prscal(6))
c
c version
c $Id: pardef.inc,v 1.1 2004/07/05 18:35:49 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 21:38:38 2003
c
c End of file 
