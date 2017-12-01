// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/FaultCohesiveNew.hh
 *
 * @brief C++ abstract base class for a fault surface implemented with
 * cohesive elements.
 */

#if !defined(pylith_faults_faultcohesivenew_hh)
#define pylith_faults_faultcohesivenew_hh

// Include directives ---------------------------------------------------
#include "pylith/feassemble/IntegratorPointwise.hh" // ISA Integrator

#include <string> // HASA std::string

// FaultCohesive --------------------------------------------------------
/// Absract base class for fault surface implemented with cohesive cells.
class pylith::faults::FaultCohesiveNew : public pylith::feassemble::IntegratorPointwise {
    friend class TestFaultCohesiveNew; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    FaultCohesiveNew(void);

    /// Destructor.
    virtual ~FaultCohesiveNew(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set material identifier of fault.
     *
     * @param[in] value Fault identifier
     */
    void id(const int value);

    /** Get material identifier of fault.
     *
     * @returns Fault identifier
     */
    int id(void) const;

    /** Set label of group of vertices associated with fault.
     *
     * @param[in] value Label of fault
     */
    void label(const char* value);

    /** Get label of group of vertices associated with fault.
     *
     * @returns Label of fault
     */
    const char* label(void) const;

    /** Set label of group of vertices defining buried edge of fault.
     *
     * @param[in] value Label of fault
     */
    void edge(const char* value);

    /** Get label of group of vertices defining buried edge of fault.
     *
     * @returns Label of fault
     */
    const char* edge(void) const;

    /** Set up direction to discriminate among shear directions in 3-D.
     *
     * @param[in] vec Up direction unit vector.
     */
    void upDir(const double vec[3]);

    /** Adjust mesh topology for fault implementation.
     *
     * @param mesh[in] PETSc mesh.
     */
    void adjustTopology(topology::Mesh* const mesh);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& solution);

    /** Initialize integrator.
     *
     * Create fault mesh from cohesive cells and cohesive point map.
     *
     * Derived class initialize, should:
     * 1. Setup subfields in auxiliary field.
     * 2. Populate auxiliary subfields.
     * 3. Set finite-element kernels.
     *
     * @param[in] solution Solution field (layout).
     */
    virtual
    void initialize(const pylith::topology::Field& solution);


    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _faultMesh; ///< Mesh over fault surface.

    PetscIS* _cohesivePointMap; ///< Map from fault point to higher dimension point in cohesive cell.

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    int _id; ///< Identifier for cohesive cells.
    std::string _label; ///< Label for vertices associated with fault.
    std::string _edge; ///< Label for vertices along buried edges of fault.
    PylithReal _upDir[3]; ///< Up direction for unique fault orientation.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    FaultCohesiveNew(const FaultCohesiveNew&); ///< Not implemented
    const FaultCohesiveNew& operator=(const FaultCohesiveNew&); ///< Not implemented

}; // class FaultCohesiveNew

#endif // pylith_faults_faultcohesivenew_hh


// End of file
