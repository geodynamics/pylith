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

/** @file libsrc/sources/WellboreSource.hh
 *
 * @brief C++ class for solving wellboresource equation.
 */

#if !defined(pylith_sources_wellboresource_hh)
#define pylith_sources_wellboresource_hh

#include "sourcesfwd.hh" // forward declarations

#include "pylith/sources/Source.hh" // ISA Source

class pylith::sources::WellboreSource : public pylith::sources::Source
{
    friend class TestWellboreSource; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:
    /// Default constructor.
    WellboreSource(void);

    /// Destructor.
    ~WellboreSource(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

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
    pylith::sources::AuxiliaryFactoryWellboreSource *_auxiliaryFactory; ///< Factory for auxiliary subfields.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:
    WellboreSource(const WellboreSource &);                  ///< Not implemented.
    const WellboreSource &operator=(const WellboreSource &); /// Not implemented.
};

// class WellboreSource

#endif // pylith_sources_wellboresource_hh

// End of file
