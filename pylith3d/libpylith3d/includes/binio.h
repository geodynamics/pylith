/*
 * 
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 *                               Charles A. Williams
 *                        Rensselaer Polytechnic Institute
 *                        (C) 2005 All Rights Reserved
 * 
 *  All worldwide rights reserved.  A license to use, copy, modify and
 *  distribute this software for non-commercial research purposes only
 *  is hereby granted, provided that this copyright notice and
 *  accompanying disclaimer is not modified or removed from the software.
 *
 *  DISCLAIMER:  The software is distributed "AS IS" without any express
 *  or implied warranty, including but not limited to, any implied
 *  warranties of merchantability or fitness for a particular purpose
 *  or any warranty of non-infringement of any current or pending patent
 *  rights.  The authors of the software make no representations about
 *  the suitability of this software for any particular purpose.  The
 *  entire risk as to the quality and performance of the software is with
 *  the user.  Should the software prove defective, the user assumes the
 *  cost of all necessary servicing, repair or correction.  In
 *  particular, neither Rensselaer Polytechnic Institute, nor the authors
 *  of the software are liable for any indirect, special, consequential,
 *  or incidental damages related to the software, to the maximum extent
 *  the law permits.
 * 
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#if !defined(lithomop3d_binio_h)
#define lithomop3d_binio_h


/*
 * The current setup uses autoconf to determine whether the fortran
 * compiler supports old-style f77 unformatted direct access, stream
 * I/O, or neither.  If neither is supported, binary output is turned
 * off.
*/

#if defined FORTRAN_F77_IO
#define BINARY_OPEN_SPEC access="direct",recl=1
#define BINARY_POS rec
#define BINARY_INQ_POS nextrec
#elif defined FORTRAN_STREAM_IO
#define BINARY_OPEN_SPEC access="stream"
#define BINARY_POS pos
#define BINARY_INQ_POS pos
#elif defined FORTRAN_NO_BINARY
#warning "Binary UCD output not supported for this Fortran compiler!"
#else
#error "Configuration error!"
#endif

#if defined FORTRAN_EOR
#define BINARY_ERR(s) advance='yes',eor=s,err=s
#else
#define BINARY_ERR(s) err=s
#endif

#endif

/* End of file */
