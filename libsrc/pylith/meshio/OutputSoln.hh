// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/OutputSoln.hh
 *
 * @brief Manager for output of solution.
 */

#if !defined(pylith_meshio_outputsoln_hh)
#define pylith_meshio_outputsoln_hh

#include "meshiofwd.hh" // forward declarations

#include "pylith/problems/ObserverSoln.hh" // ISA Observer
#include "pylith/meshio/OutputObserver.hh" // ISA OutputObserver

#include "pylith/utils/array.hh" // HASA string_vector

class pylith::meshio::OutputSoln :
    public pylith::problems::ObserverSoln,
    public pylith::meshio::OutputObserver {
    friend class TestOutputSoln; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    OutputSoln(void);

    /// Destructor
    virtual ~OutputSoln(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set names of solution subfields requested for output.
     *
     * @param[in] names Array of subfield names.
     * @param[in] numNames Length of array.
     */
    void setOutputSubfields(const char* names[],
                            const int numNames);

    /** Get names of solution subfields requested for output.
     *
     * @returns Array of subfield names.
     */
    const pylith::string_vector& getOutputSubfields(void) const;

    /** Set time scale.
     *
     * @param[in] value Time scale for dimensionalizing time.
     */
    void setTimeScale(const PylithReal value);

    /** Verify observer is compatible with solution.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Receive update from subject.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     * @param[in] infoOnly Flag is true if this update is before solution is available (e.g., after initialization).
     */
    void update(const PylithReal t,
                const PylithInt tindex,
                const pylith::topology::Field& solution);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Prepare for output.
     *
     * @param[in] mesh Finite-element mesh object.
     */
    virtual
    void _open(const pylith::topology::Mesh& mesh);

    /// Close output files.
    virtual
    void _close(void);

    /** Prepare for output at this solution step.
     *
     * @param[in] t Time associated with solution step.
     * @param[in] mesh Mesh for output.
     */
    virtual
    void _openSolnStep(const PylithReal t,
                       const pylith::topology::Mesh& mesh);

    /// Finalize output at this solution step.
    virtual
    void _closeSolnStep(void);

    /** Write output for step in solution.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     */
    virtual
    void _writeSolnStep(const PylithReal t,
                        const PylithInt tindex,
                        const pylith::topology::Field& solution);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::string_vector _subfieldNames;

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    OutputSoln(const OutputSoln&); ///< Not implemented.
    const OutputSoln& operator=(const OutputSoln&); ///< Not implemented

}; // OutputSoln

#endif // pylith_meshio_outputsoln_hh

// End of file
