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


cdef extern from "const.h":
    ctypedef char const_char # Pyrex v0.9.5 doesn't support 'const'


cdef extern from "petsc.h":
    ctypedef int PetscErrorCode
    ctypedef int PetscEvent
    ctypedef int PetscCookie
    ctypedef int PetscInt
    PetscErrorCode PetscInitialize(int *, char ***, const_char *, const_char *)
    PetscErrorCode PetscFinalize()


cdef extern from "petscerror.h":
    PetscErrorCode PetscErrorMessage(int, const_char **, char **)


cdef extern from "petscvec.h":
    struct _p_Vec
    ctypedef _p_Vec *Vec


cdef extern from "petscmat.h":
    struct _p_Mat
    ctypedef _p_Mat *Mat


cdef extern from "petscmeshfwd.h": # for compilation speed
    struct _p_Mesh
    ctypedef _p_Mesh *Mesh
    PetscErrorCode MeshDestroy(Mesh)


cdef extern from "petsclog.h":
    PetscErrorCode PetscLogStageRegister(int *, char *)
    PetscErrorCode PetscLogEventRegister(PetscEvent *, char *, PetscCookie)
    PetscErrorCode PetscLogStagePush(int stage)
    PetscErrorCode PetscLogStagePop()


cdef extern from "petscsnes.h":
    enum:
        KSP_COOKIE


# end of file
