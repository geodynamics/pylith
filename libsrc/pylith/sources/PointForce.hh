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

/** @file libsrc/sources/PointForce.hh
 *
 * @brief C++ class for solving pointforce equation.
 */

#if !defined(pylith_sources_pointforce_hh)
#define pylith_sources_pointforce_hh

#include "sourcesfwd.hh" // forward declarations

#include "pylith/sources/Source.hh" // ISA Source

class pylith::sources::PointForce : public pylith::sources::Source {
    friend class TestPointForce; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    PointForce(void);

    /// Destructor.
    ~PointForce(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set time history database.
     *
     * @param[in] db Time history database.
     */
    void setTimeHistoryDB(spatialdata::spatialdb::TimeHistory *th);

    /** Get time history database.
     *
     * @preturns Time history database.
     */
    const spatialdata::spatialdb::TimeHistory *getTimeHistoryDB(void);

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

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field &solution) const;

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     *
     *  @returns Integrator if applicable, otherwise NULL.
     */
    pylith::feassemble::Integrator *createIntegrator(const pylith::topology::Field &solution);

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    pylith::topology::Field *createAuxiliaryField(const pylith::topology::Field &solution,
                                                  const pylith::topology::Mesh &domainMesh);

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Derived field if applicable, otherwise NULL.
     */
    pylith::topology::Field *createDerivedField(const pylith::topology::Field &solution,
                                                const pylith::topology::Mesh &domainMesh);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory *_getAuxiliaryFactory(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    /** Set kernels for residual.
     *
     * @param[out] integrator Integrator for source.
     * @param[in] solution Solution field.
     */
    void _setKernelsResidual(pylith::feassemble::IntegratorDomain *integrator,
                             const pylith::topology::Field &solution) const;

    /** Set kernels for Jacobian.
     *
     * @param[out] integrator Integrator for source.
     * @param[in] solution Solution field.
     */
    void _setKernelsJacobian(pylith::feassemble::IntegratorDomain *integrator,
                             const pylith::topology::Field &solution) const;

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    spatialdata::spatialdb::TimeHistory *_dbTimeHistory; ///< Time history database.
    pylith::sources::AuxiliaryFactoryPointForce *_auxiliaryFactory; ///< Factory for auxiliary subfields.

    bool _useTimeHistory; ///< Use time history term.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    PointForce(const PointForce &); ///< Not implemented.
    const PointForce &operator=(const PointForce &); /// Not implemented.

};

// class PointForce

#endif // pylith_sources_pointforce_hh

// End of file
