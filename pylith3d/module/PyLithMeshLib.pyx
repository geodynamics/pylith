#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


# implementation of PyLithMeshLib


cdef public class Mesh [object PyMeshObject, type PyMesh_Type]:


    def __new__(self,
                meshInputFile,
                meshBcFile,
                interpolateMesh,
                partitioner):
        self.mesh = NULL
        self.A = NULL
        self.rhs = NULL
        self.sol = NULL
        self.meshInputFile = NULL
        self.meshBcFile = NULL
        self.interpolateMesh = 0
        self.partitioner = NULL
        return


    def __init__(self,
                 meshInputFile,
                 meshBcFile,
                 interpolateMesh,
                 partitioner):
        
        # maintain references to Python strings
        self._meshInputFile = meshInputFile   # Py_INCREF
        self._meshBcFile = meshBcFile  # Py_INCREF
        self._partitioner = partitioner  # Py_INCREF

        # fetch raw char * pointers
        self.meshInputFile = self._meshInputFile
        self.meshBcFile = self._meshBcFile
        self.partitioner = self._partitioner

        self.interpolateMesh = interpolateMesh

        _processMesh(self)

        return


    def __dealloc__(self):
        if self.A is not NULL:
            _destroyMat(self)
        if self.mesh is not NULL:
            petsc.MeshDestroy(self.mesh)
            self.mesh = NULL
        return


    cdef createMat(self):
        _createMat(self)


    cdef destroyMat(self):
        _destroyMat(self)


    cdef outputMesh(self, fileRoot):
        _outputMesh(self, fileRoot)



# C++ implementation

cdef extern from "mesh.h":
    cdef extern object _processMesh "PyLithMeshLib::Mesh::_processMesh" (Mesh self)
    cdef extern object _createMat   "PyLithMeshLib::Mesh::_createMat "  (Mesh self)
    cdef extern object _destroyMat  "PyLithMeshLib::Mesh::_destroyMat"  (Mesh self)
    cdef extern object _outputMesh  "PyLithMeshLib::Mesh::_outputMesh"  (Mesh self, char *fileRoot)


# end of file
