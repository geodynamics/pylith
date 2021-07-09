// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "Xdmf.hh" // implementation of class methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include <Python.h>

#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

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
        throw std::runtime_error("Could not import module 'pylith.meshio.Xdmf'.");
    } // if
    PyObject* cls = PyObject_GetAttrString(mod, "Xdmf");
    if (!cls) {
        throw std::runtime_error("Could not get 'Xdmf' attribute in pylith.meshio.Xdmf module.");
    } // if
    PyObject* pyXdmf = PyObject_CallFunctionObjArgs(cls, NULL);
    if (!pyXdmf) {
        throw std::runtime_error("Could not create Python Xdmf object.");
    } // if
    Py_DECREF(cls);
    Py_DECREF(mod);

    PyObject* pyWrite = PyObject_CallMethod(pyXdmf, const_cast<char*>("write"), const_cast<char*>("s"), filenameH5);
    if (!pyWrite) {
        if (PyErr_Occurred()) {
            PyErr_Clear();
        } // if
        std::ostringstream msg;
        msg << "Could not generate Xdmf file for HDF5 file '" << filenameH5 << "'.";
        throw std::runtime_error(msg.str());
    } // if
    Py_DECREF(pyWrite);
    Py_CLEAR(pyXdmf);

    if (!alreadyInitialized) {
        Py_Finalize();
    } // if

    PYLITH_METHOD_END;
} // write


// End of file
