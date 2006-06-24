/*
 * 
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
 *
 *  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *
 *  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 *  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 *  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
 *  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#if !defined(pylith3d_binio_h)
#define pylith3d_binio_h


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
