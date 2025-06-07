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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/sources/MomentTensorForce.hh
 *
 * @brief C++ class for solving momenttensorforce equation.
 */

#if !defined(pylith_sources_momenttensorforce_hh)
#define pylith_sources_momenttensorforce_hh

#include "sourcesfwd.hh" // forward declarations

#include "pylith/sources/Source.hh" // ISA Source

class pylith::sources::MomentTensorForce : public pylith::sources::Source {
    friend class TestMomentTensorForce; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    MomentTensorForce(void);

    /// Destructor.
    ~MomentTensorForce(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set source time function.
     *
     * @param[in] sourceTimeFunction Source time function for momenttensor force.
     */
    void setSourceTimeFunction(pylith::sources::SourceTimeFunctionMomentTensorForce* const sourceTimeFunction);

    /** Get source time function.
     *
     * @returns Source time function for momenttensor force.
     */
    pylith::sources::SourceTimeFunctionMomentTensorForce* getSourceTimeFunction(void) const;

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Set time history database.
     *
     * @param[in] db Time history database.
     */
    void setTimeHistoryDB(spatialdata::spatialdb::TimeHistory* th);

    /** Get time history database.
     *
     * @preturns Time history database.
     */
    const spatialdata::spatialdb::TimeHistory* getTimeHistoryDB(void);

    /** Use time history term in time history expression.
     *
     * @param[in] value True if using time history term in expression.
     */
    void useTimeHistory(const bool value);

    /** Get flag associated with using time history term in time history expression.
     *
     * @returns True if using time history term in expression, false otherwise.
     */
    bool useTimeHistory(void) const;

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     *
     *  @returns Integrator if applicable, otherwise NULL.
     */
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution);

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& domainMesh);

    /** Update auxiliary subfields at beginning of time step.
     *
     * @param[out] auxiliaryField Auxiliary field.
     * @param[in] t Current time.
     */
    void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                              const double t);

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Derived field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& domainMesh);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    /** Update kernel constants.
     *
     * @param[in] dt Current time step.
     */
    void _updateKernelConstants(const PylithReal dt);

    /** Get derived factory associated with physics.
     *
     * @return Derived factory for physics object.
     */
    pylith::topology::FieldFactory* _getDerivedFactory(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    /** Set kernels for residual.
     *
     * @param[out] integrator Integrator for source.
     * @param[in] solution Solution field.
     */
    void _setKernelsResidual(pylith::feassemble::IntegratorDomain* integrator,
                             const pylith::topology::Field& solution) const;

    // /** Set kernels for Jacobian.
    //  *
    //  * @param[out] integrator Integrator for source.
    //  * @param[in] solution Solution field.
    //  */
    // void _setKernelsJacobian(pylith::feassemble::IntegratorDomain* integrator,
    //                          const pylith::topology::Field& solution) const;

    /** Set kernels for computing updated state variables in auxiliary field.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsUpdateStateVars(pylith::feassemble::IntegratorDomain* integrator,
                                    const pylith::topology::Field& solution) const;

    /** Set kernels for computing derived field.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    void _setKernelsDerivedField(pylith::feassemble::IntegratorDomain* integrator,
                                 const pylith::topology::Field& solution) const;

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::sources::SourceTimeFunctionMomentTensorForce* _sourceTimeFunction; ///< Source time function for
                                                                               ///< momenttensor force.
    pylith::sources::DerivedFactoryMomentTensorForce* _derivedFactory; ///< Factory for creating derived fields.


    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    MomentTensorForce(const MomentTensorForce&); ///< Not implemented.
    const MomentTensorForce& operator=(const MomentTensorForce&); /// Not implemented.

};

// class MomentTensorForce

#endif // pylith_sources_momenttensorforce_hh

// End of file
