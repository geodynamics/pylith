// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// See LICENSE.md.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/FaultCohesive.hh
 *
 * @brief C++ abstract base class for a fault interface implemented with
 * cohesive cells.
 */

#if !defined(pylith_faults_faultcohesive_hh)
#define pylith_faults_faultcohesive_hh

#include "faultsfwd.hh" // forward declarations

#include "pylith/problems/Physics.hh" // ISA Physics

#include "pylith/materials/materialsfwd.hh" // USES Material

#include <string> // HASA std::string

class pylith::faults::FaultCohesive : public pylith::problems::Physics {
    friend class TestFaultCohesive; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    FaultCohesive(void);

    /// Destructor.
    virtual ~FaultCohesive(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set identifier for fault cohesive cells.
     *
     * @param[in] value Fault identifier
     */
    void setInterfaceId(const int value);

    /** Get identifier for fault cohesive cells.
     *
     * @returns Fault identifier
     */
    int getInterfaceId(void) const;

    /** Set label marking surface of interface.
     *
     * @param[in] value Label of surface (from mesh generator).
     */
    void setSurfaceMarkerLabel(const char* value);

    /** Get label marking surface of interface.
     *
     * @returns Label of surface (from mesh generator).
     */
    const char* getSurfaceMarkerLabel(void) const;

    /** Set label marking buried edges of interface surface.
     *
     * @param[in] value Label of buried surface edge (from mesh generator).
     */
    void setBuriedEdgesMarkerLabel(const char* value);

    /** Get label marking buried edges of interface surface.
     *
     * @returns Label of buried surface edge (from mesh generator).
     */
    const char* getBuriedEdgesMarkerLabel(void) const;

    /** Set first choice for reference direction to discriminate among tangential directions in 3-D.
     *
     * @param vec Reference direction unit vector.
     */
    void setRefDir1(const PylithReal vec[3]);

    /** Set second choice for reference direction to discriminate among tangential directions in 3-D.
     *
     * @param vec Reference direction unit vector.
     */
    void setRefDir2(const PylithReal vec[3]);

    /** Adjust mesh topology for fault implementation.
     *
     * @param mesh[in] PETSc mesh.
     */
    void adjustTopology(pylith::topology::Mesh* const mesh);

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     * @param[in] materials Materials in problem.
     * @returns Integrator if applicable, otherwise NULL.
     */
    virtual
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution,
                                                     const std::vector<pylith::materials::Material*>& materials) = 0;

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Create single integratin patch for entire fault.
     *
     * @param[inout] integrator Integrator for fault interface.
     */
    void _createIntegrationPatch(pylith::feassemble::IntegratorInterface* integrator);

    /** Create integration patches associated with cohesive cells that have the same pairs of materials on the
     * two sides of the fault.
     *
     * @param[inout] integrator Integrator for fault interface.
     * @param[in] dmSoln PETSc DM associated with solution.
     */
    void _createIntegrationPatches(pylith::feassemble::IntegratorInterface* integrator,
                                   const PetscDM dmSoln);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    PylithReal _refDir1[3]; ///< First choice reference direction used to compute boundary tangential directions.
    PylithReal _refDir2[3]; ///< Second choice reference direction used to compute boundary tangential directions.

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    int _interfaceId; ///< Identifier for cohesive cells.
    std::string _interfaceLabel; ///< Label identifying vertices associated with fault.
    std::string _buriedEdgesLabel; ///< Label identifying vertices along buried edges of fault.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution); // Empty method

    FaultCohesive(const FaultCohesive&); ///< Not implemented
    const FaultCohesive& operator=(const FaultCohesive&); ///< Not implemented

}; // class FaultCohesive

#endif // pylith_faults_faultcohesive_hh

// End of file
