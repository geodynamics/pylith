// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/utils/utilsfwd.hh" // forward declarations

// Version ----------------------------------------------------------
/** @brief C++ object for getting version info.
 */
class pylith::utils::DependenciesVersion { // DependenciesVersion
    friend class TestDependenciesVersion; // unit testing

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /// Default constructor.
    DependenciesVersion(void);

    /// Default destrictor.
    ~DependenciesVersion(void);

    /** Get MPI version number.
     *
     * @returns MPI version number.
     */
    static
    const char* mpiVersion(void);

    /** Get MPI implemenation info (OpenMPI, MPICH, etc).
     *
     * @returns MPI implementation info.
     */
    static
    const char* mpiImplementation(void);

    /** Get MPI standard version info.
     *
     * @returns MPI standard version info.
     */
    static
    const char* mpiStandard(void);

    /** Get NetCDF version number.
     *
     * @returns NetCDF version number.
     */
    static
    const char* netcdfVersion(void);

    /** Get HDF5 version number.
     *
     * @returns HDF5 version number.
     */
    static
    const char* hdf5Version(void);

    // PRIVATE METHODS //////////////////////////////////////////////////////
private:

    DependenciesVersion(const DependenciesVersion&); ///< Not implemented
    const DependenciesVersion& operator=(const DependenciesVersion&); ///< Not implemented

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    static const char* _mpiImplementation; ///< MPI implementation
    static const char* _mpiVersion; ///< MPI version number.
    static const char* _mpiStandard; ///< MPI standard version number.

    static const char* _netcdfVersion; ///< NetCDF version number.

    static const char* _hdf5Version; ///< HDF5 version number.

}; // DependenciesVersion

// End of file
