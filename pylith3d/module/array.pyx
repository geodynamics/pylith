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


cdef int memoryUsage


cdef class IntArray:
    
    cdef int len
    cdef int *ptr
    
    def __new__(self, l):
        cdef object init
        cdef long int size
        if isinstance(l, int):
            self.len = l
        else:
            self.len = len(l)
            init = l
        size = sizeof(int) * self.len
        self.ptr = <int *>malloc(size)
        if self.ptr is NULL:
            raise MemoryError("malloc failed")
        global memoryUsage
        memoryUsage = memoryUsage + size
        if init:
            for i from 0 <= i < self.len:
                self.ptr[i] = init[i]
        return

    def __dealloc__(self):
        free(self.ptr)
        self.ptr = NULL
        global memoryUsage
        memoryUsage = memoryUsage - (sizeof(int) * self.len)
        return

    cdef asList(self):
        return listFromIntArray(self.ptr, self.len)


cdef class DoubleArray:
    
    cdef int len
    cdef double *ptr
    
    def __new__(self, l):
        cdef object init
        cdef long int size
        if isinstance(l, int):
            self.len = l
        else:
            self.len = len(l)
            init = l
        size = sizeof(double) * self.len
        self.ptr = <double *>malloc(size)
        if self.ptr is NULL:
            raise MemoryError("malloc failed")
        global memoryUsage
        memoryUsage = memoryUsage + size
        if init:
            for i from 0 <= i < self.len:
                self.ptr[i] = init[i]
        return

    def __dealloc__(self):
        free(self.ptr)
        self.ptr = NULL
        global memoryUsage
        memoryUsage = memoryUsage - (sizeof(double) * self.len)
        return

    cdef asList(self):
        return listFromDoubleArray(self.ptr, self.len)


# utilities


cdef setIntArrayFromList(int *array, int dim, aList):
    assert(len(aList) == dim)
    for i from 0 <= i < dim:
        array[i] = aList[i]
    return


cdef listFromIntArray(int *array, int dim):
    aList = []
    for i from 0 <= i < dim:
        aList.append(array[i])
    return aList


cdef setDoubleArrayFromList(double *array, int dim, aList):
    assert(len(aList) == dim)
    for i from 0 <= i < dim:
        array[i] = aList[i]
    return


cdef listFromDoubleArray(double *array, int dim):
    aList = []
    for i from 0 <= i < dim:
        aList.append(array[i])
    return aList


# end of file
