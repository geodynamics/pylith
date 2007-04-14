# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
#
#  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
#
#  Permission is hereby granted, free of charge, to any person obtaining
#  a copy of this software and associated documentation files (the
#  "Software"), to deal in the Software without restriction, including
#  without limitation the rights to use, copy, modify, merge, publish,
#  distribute, sublicense, and/or sell copies of the Software, and to
#  permit persons to whom the Software is furnished to do so, subject to
#  the following conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
#  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
#  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
#  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#


cdef extern from "stdlib.h":
    void *malloc(unsigned int)
    void free(void *)


cdef class IntArray:
    
    cdef int len
    cdef int *ptr
    cdef int *_rawptr
    
    def __new__(self, unsigned int len):
        self.len = len
        self._rawptr = <int *>malloc(sizeof(int) * (len + 2))
        if self._rawptr is NULL:
            raise MemoryError("malloc failed")
        self._rawptr[0] = 0xdeadbeef
        self._rawptr[len + 1] = 0xdeadbeef
        self.ptr = self._rawptr + 1
        return

    def __dealloc__(self):
        if (self._rawptr[0] != 0xdeadbeef or
            self._rawptr[self.len + 1] != 0xdeadbeef):
            print "memory corrupted"
        free(self._rawptr)
        self._rawptr = NULL
        self.ptr = NULL
        return


cdef class DoubleArray:
    
    cdef int len
    cdef double *ptr
    cdef double *_rawptr
    cdef double guard
    
    def __new__(self, unsigned int len):
        self.len = len
        self._rawptr = <double *>malloc(sizeof(double) * (len + 2))
        if self._rawptr is NULL:
            raise MemoryError("malloc failed")
        self.guard = 3.14159
        self._rawptr[0] = self.guard
        self._rawptr[len + 1] = self.guard
        self.ptr = self._rawptr + 1
        return

    def __dealloc__(self):
        if (self._rawptr[0] != self.guard or
            self._rawptr[self.len + 1] != self.guard):
            print "memory corrupted"
        free(self._rawptr)
        self._rawptr = NULL
        self.ptr = NULL
        return


# end of file
