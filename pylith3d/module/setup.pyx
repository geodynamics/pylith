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


#cimport petsc

cdef extern from "string.h":
    char *strcpy(char *, char *)


gfile = None
ghelp = None


def PetscInitialize(args, file=None, help=None):
    cdef int argc
    cdef char **cargs

    argc = len(args)
    cargs = <char **>malloc((argc + 1) * sizeof(char *))

    cdef int i
    cdef char *carg
    for i from 0 <= i < argc:
        arg = args[i]
        carg = arg
        cargs[i] = <char *>malloc((len(arg)+1) * sizeof(char))
        strcpy(cargs[i], carg)
    
    cargs[argc] = NULL

    cdef char *cfile, *chelp
    cfile = NULL
    chelp = NULL
    if file is not None:
        global gfile
        gfile = file  # Py_INCREF
        cfile = gfile
    if help is not None:
        global ghelp
        ghelp = help  # Py_INCREF
        chelp = ghelp

    cdef int errnum
    cdef petsc.const_char *text
    cdef char *specific
    errnum = petsc.PetscInitialize(&argc, &cargs, cfile, chelp)
    if errnum != 0:
        petsc.PetscErrorMessage(errnum, &text, &specific)
        raise RuntimeError("PetscInitialize: error %d: %s (%s)" %
                           (errnum, <char*>text, specific))

    # Do not free 'cargs' -- PETSc saves a reference to it ('PetscGlobalArgs').

    return


def PetscFinalize():
    cdef int errnum
    cdef petsc.const_char *text
    cdef char *specific
    
    errnum = petsc.PetscFinalize()
    if errnum != 0:
        petsc.PetscErrorMessage(errnum, &text, &specific)
        raise RuntimeError("PetscFinalize: error %d: %s (%s)" %
                           (errnum, <char*>text, specific))
    
    return


cdef petsc.PetscInt PetscLogStageRegister(sname):
    cdef int errnum
    cdef petsc.const_char *text
    cdef char *specific
    cdef petsc.PetscInt stage
    
    errnum = petsc.PetscLogStageRegister(&stage, sname)
    if errnum != 0:
        petsc.PetscErrorMessage(errnum, &text, &specific)
        raise RuntimeError("PetscLogStageRegister('%s'): error %d: %s (%s)" %
                           (sname, errnum, <char*>text, specific))

    return stage


cdef petsc.PetscEvent PetscLogEventRegister(name, petsc.PetscCookie cookie):
    cdef int errnum
    cdef petsc.const_char *text
    cdef char *specific
    cdef petsc.PetscEvent event
    
    errnum = petsc.PetscLogEventRegister(&event, name, cookie)
    if errnum != 0:
        petsc.PetscErrorMessage(errnum, &text, &specific)
        raise RuntimeError("PetscLogEventRegister('%s'): error %d: %s (%s)" %
                           (name, errnum, <char*>text, specific))

    return event


cdef PetscLogStagePush(int stage):
    cdef int errnum
    cdef petsc.const_char *text
    cdef char *specific

    errnum = petsc.PetscLogStagePush(stage)
    if errnum != 0:
        petsc.PetscErrorMessage(errnum, &text, &specific)
        raise RuntimeError("PetscLogStagePush: error %d: %s (%s)" %
                           (errnum, <char*>text, specific))

    return


cdef PetscLogStagePop():
    cdef int errnum
    cdef petsc.const_char *text
    cdef char *specific
    
    errnum = petsc.PetscLogStagePop()
    if errnum != 0:
        petsc.PetscErrorMessage(errnum, &text, &specific)
        raise RuntimeError("PetscLogStagePop: error %d: %s (%s)" %
                           (errnum, <char*>text, specific))
    
    return


# end of file
