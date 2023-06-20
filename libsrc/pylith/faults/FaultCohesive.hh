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
// Copyright (c) 2010-2022 University of California, Davis
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

    /** Set name of label identifying cohesive cells.
     *
     * @param[in] value Name of label.
     */
    void setCohesiveLabelName(const char* value);

    /** Get name of label identifying cohesive cells.
     *
     * @returns Name of label.
     */
    const char* getCohesiveLabelName(void) const;

    /** Set value of label identifying cohesive cells.
     *
     * @param[in] value Value of label.
     */
    void setCohesiveLabelValue(const int value);

    /** Get value of label identifying cohesive cells.
     *
     * @returns Value of label.
     */
    int getCohesiveLabelValue(void) const;

    /** Set name of label marking surface of interface.
     *
     * @param[in] value Name of label of surface (from mesh generator).
     */
    void setSurfaceLabelName(const char* value);

    /** Get name of label marking surface of interface.
     *
     * @returns Name of label of surface (from mesh generator).
     */
    const char* getSurfaceLabelName(void) const;

    /** Set value of label marking surface of interface.
     *
     * @param[in] value Value of label of surface (from mesh generator).
     */
    void setSurfaceLabelValue(const int value);

    /** Get value of label marking surface of interface.
     *
     * @returns Value of label of surface (from mesh generator).
     */
    int getSurfaceLabelValue(void) const;

    /** Set name of label marking buried edges of interface surface.
     *
     * @param[in] value Name of label of buried surface edge (from mesh generator).
     */
    void setBuriedEdgesLabelName(const char* value);

    /** Get name of label marking buried edges of interface surface.
     *
     * @returns Name of label of buried surface edge (from mesh generator).
     */
    const char* getBuriedEdgesLabelName(void) const;

    /** Set value of label marking buried edges of interface surface.
     *
     * @param[in] value Value of label of buried surface edge (from mesh generator).
     */
    void setBuriedEdgesLabelValue(const int value);

    /** Get value of label marking buried edges of interface surface.
     *
     * @returns Value of label of buried surface edge (from mesh generator).
     */
    int getBuriedEdgesLabelValue(void) const;

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
                                                     const std::vector<pylith::materials::Material*>& materials);

    /** Create constraint and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Constraint if applicable, otherwise NULL.
     */
    std::vector<pylith::feassemble::Constraint*> createConstraints(const pylith::topology::Field& solution);

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Derived field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& domainMesh);

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

    /** Update kernel constants.
     *
     * @param[in] dt Current time step.
     */
    void _updateKernelConstants(const PylithReal dt);

    /** Set kernels for residual.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     * @param[in] materials Materials in problem.
     */
    virtual
    void _setKernelsResidual(pylith::feassemble::IntegratorInterface* integrator,
                             const pylith::topology::Field& solution,
                             const std::vector<pylith::materials::Material*>& materials) const = 0;

    /** Set kernels for Jacobian.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     * @param[in] materials Materials in problem.
     */
    virtual
    void _setKernelsJacobian(pylith::feassemble::IntegratorInterface* integrator,
                             const pylith::topology::Field& solution,
                             const std::vector<pylith::materials::Material*>& materials) const = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    PylithReal _refDir1[3]; ///< First choice reference direction used to compute boundary tangential directions.
    PylithReal _refDir2[3]; ///< Second choice reference direction used to compute boundary tangential directions.

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    inline
    static
    PetscErrorCode _zero(PetscInt dim,
                         PetscReal t,
                         const PetscReal x[],
                         PetscInt Nc,
                         PetscScalar *u,
                         void *ctx) {
        for (int c = 0; c < Nc; ++c) {
            u[c] = 0.0;
        } // for
        return 0;
    } // _zero

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    std::string _surfaceLabelName; ///< Name of label identifying points associated with fault.
    std::string _buriedEdgesLabelName; ///< Name of label identifying buried edges of fault.
    int _surfaceLabelValue; ///< Value of label identifying points associated with fault.
    int _buriedEdgesLabelValue; ///< Value of label identifying buried edges of fault.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    FaultCohesive(const FaultCohesive&); ///< Not implemented
    const FaultCohesive& operator=(const FaultCohesive&); ///< Not implemented

    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution); // Empty method

}; // class FaultCohesive

#endif // pylith_faults_faultcohesive_hh

// End of file
