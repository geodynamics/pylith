// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/meshio/Xdmf.hh" // implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // PYLITH_ journal macros
#include "pylith/utils/Exceptions.hh" // USES Exception

#include <Python.h>

// ----------------------------------------------------------------------
// Use Python Xdmf object to write Xdmf file corresponding to HDF5 file.
void
pylith::meshio::Xdmf::write(const char* filenameH5) {
    PYLITH_METHOD_BEGIN;

    const bool alreadyInitialized = Py_IsInitialized();
    if (!alreadyInitialized) {
        Py_Initialize();
    } // if

    /* Replicate Python functionality given by
     *
     * import pylith.meshio.Xdmf
     * xdmf = pylith.meshio.Xdmf()
     * xdmf.write(filenameH5)
     */

    // Should check for NULL, decode the exception, and throw a C++ equivalent
    PyObject* mod = PyImport_ImportModule("pylith.meshio.Xdmf");
    if (!mod) {
        PYLITH_ERROR(pylith::InternalError, pylith::journal::internal,
                     "Could not import module 'pylith.meshio.Xdmf'.");
    } // if
    PyObject* cls = PyObject_GetAttrString(mod, "Xdmf");
    if (!cls) {
        PYLITH_ERROR(pylith::InternalError, pylith::journal::internal,
                     "Could not get 'Xdmf' attribute in pylith.meshio.Xdmf module.");
    } // if
    PyObject* pyXdmf = PyObject_CallFunctionObjArgs(cls, NULL);
    if (!pyXdmf) {
        PYLITH_ERROR(pylith::InternalError, pylith::journal::internal,
                     "Could not create Python Xdmf object.");
    } // if
    Py_DECREF(cls);
    Py_DECREF(mod);

    PyObject* pyWrite = PyObject_CallMethod(pyXdmf, const_cast<char*>("write"), const_cast<char*>("s"), filenameH5);
    if (!pyWrite) {
        if (PyErr_Occurred()) {
            PyErr_Clear();
        } // if
        PYLITH_ERROR(pylith::InternalError, pylith::journal::internal,
                     "Could not generate Xdmf file for HDF5 file '" << filenameH5 << "'.");
    } // if
    Py_DECREF(pyWrite);
    Py_CLEAR(pyXdmf);

    if (!alreadyInitialized) {
        Py_Finalize();
    } // if

    PYLITH_METHOD_END;
} // write


// End of file
