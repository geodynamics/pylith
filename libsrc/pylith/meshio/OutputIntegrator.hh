// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/OutputIntegrator.hh
 *
 * @brief C++ object for managing output of the solution over the domain.
 */

#if !defined(pylith_meshio_outputintegrator_hh)
#define pylith_meshio_outputintegrator_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "OutputManager.hh" // ISA OutputManager

#include "pylith/utils/array.hh" // HASA string_vector

// OutputIntegrator -----------------------------------------------------
/** @brief C++ object for managing output of the solution over the domain.
 */
class pylith::meshio::OutputIntegrator : public OutputManager {
    friend class TestOutputIntegrator;   // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    OutputIntegrator(void);

    /// Destructor
    virtual ~OutputIntegrator(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set names of information related fields to output.
     *
     * @param[in] names Array of names of fields to output.
     * @param[in] numNames Length of array.
     */
    void vertexInfoFields(const char* names[],
                          const int numNames);

    /** Verify configuration.
     *
     * @param[in] auxField Auxiliary field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& auxField) const;

    /** Write information.
     *
     * @param[in] auxField Auxiliary field.
     */
    virtual
    void writeInfo(const pylith::topology::Field& auxField);

    /** Write solution at time step.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     */
    virtual
    void writeTimeStep(const PylithReal t,
                       const PylithInt tindex,
                       const pylith::topology::Field& solution);

    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    pylith::string_vector _vertexInfoFields; ///< Names of information fields to output.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    OutputIntegrator(const OutputIntegrator&);   ///< Not implemented.
    const OutputIntegrator& operator=(const OutputIntegrator&);   ///< Not implemented

}; // OutputIntegrator

#endif // pylith_meshio_outputintegrator_hh

// End of file
