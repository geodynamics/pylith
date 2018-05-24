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

/** @file libsrc/faults/FaultCohesive.hh
 *
 * @brief C++ abstract base class for a fault surface implemented with
 * cohesive elements.
 */

#if !defined(pylith_faults_faultcohesive_hh)
#define pylith_faults_faultcohesive_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

#include "pylith/feassemble/IntegratorPointwise.hh" // ISA Integrator

#include <string> // HASA std::string

// FaultCohesive --------------------------------------------------------
/// Absract base class for fault surface implemented with cohesive cells.
class pylith::faults::FaultCohesive : public pylith::feassemble::IntegratorPointwise {
    friend class TestFaultCohesive; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    FaultCohesive(void);

    /// Destructor.
    virtual ~FaultCohesive(void);

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

    /** Set first choice for reference direction to discriminate among tangential directions in 3-D.
     *
     * @param vec Reference direction unit vector.
     */
    void refDir1(const PylithReal vec[3]);

    /** Set second choice for reference direction to discriminate among tangential directions in 3-D.
     *
     * @param vec Reference direction unit vector.
     */
    void refDir2(const PylithReal vec[3]);

    /** Adjust mesh topology for fault implementation.
     *
     * @param mesh[in] PETSc mesh.
     */
    void adjustTopology(pylith::topology::Mesh* const mesh);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Initialize fault.
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


    // PROTECTED NETHODS //////////////////////////////////////////////////
protected:

    /** Get factory for setting up auxliary fields.
     *
     * @returns Factor for auxiliary fields.
     */
    pylith::feassemble::AuxiliaryFactory* _auxFactory(void);

    /** Setup auxiliary subfields (discretization and query fns).
     *
     * Create subfields in auxiliary fields (includes name of the field,
     * vector field type, discretization, and scale for
     * nondimensionalization) and set query functions for filling them
     * from a spatial database.
     *
     * @attention The order of the calls to subfieldAdd() must match the
     * order of the auxiliary fields in the FE kernels.
     */
    virtual
    void _auxFieldSetup(void) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _faultMesh; ///< Mesh over fault surface.
    PetscIS _cohesivePointMap; ///< Map from fault point to higher dimension point in cohesive cell.
    pylith::faults::AuxiliaryFactory* _auxFaultFactory; ///< Factory for auxiliary subfields.


    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    int _id; ///< Identifier for cohesive cells.
    std::string _label; ///< Label for vertices associated with fault.
    std::string _edge; ///< Label for vertices along buried edges of fault.
    PylithReal _refDir1[3]; ///< First choice reference direction used to compute boundary tangential directions.
    PylithReal _refDir2[3]; ///< Second choice reference direction used to compute boundary tangential directions.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    FaultCohesive(const FaultCohesive&); ///< Not implemented
    const FaultCohesive& operator=(const FaultCohesive&); ///< Not implemented

}; // class FaultCohesive

#endif // pylith_faults_faultcohesive_hh


// End of file
