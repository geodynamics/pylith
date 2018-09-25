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
 * @brief C++ object for managing output of the solution over the domain of an integrator.
 */

#if !defined(pylith_meshio_outputintegrator_hh)
#define pylith_meshio_outputintegrator_hh

#include "meshiofwd.hh" // forward declarations

#include "OutputManager.hh" // ISA OutputManager

#include "pylith/utils/array.hh" // HASA string_vector

class pylith::meshio::OutputIntegrator : public pylith::meshio::OutputManager {
    friend class TestOutputIntegrator; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor
     *
     * @param[in] integrator Integrator to observe.
     */
    OutputIntegrator(pylith::feassemble::Integrator* const integrator);

    /// Destructor
    ~OutputIntegrator(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Verify configuration.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /// Write solution at time step.
    void _writeInfo(void);

    /** Write solution at time step.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     */
    void _writeDataStep(const PylithReal t,
                        const PylithInt tindex,
                        const pylith::topology::Field& solution);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::feassemble::Integrator* const _integrator;

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    OutputIntegrator(const OutputIntegrator&); ///< Not implemented.
    const OutputIntegrator& operator=(const OutputIntegrator&); ///< Not implemented

}; // OutputIntegrator

#endif // pylith_meshio_outputintegrator_hh

// End of file
