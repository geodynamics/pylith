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
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/sources/RickerWavelet.hh
 *
 * @brief C++ class for Ricker source time function.
 */

#if !defined(pylith_materials_rickerwavelet_hh)
#define pylith_materials_rickerwavelet_hh

#include "sourcesfwd.hh" // forward declarations

#include "pylith/sources/SourceTimeFunctionMomentTensorForce.hh" // ISA SourceTimeFunctionMomentTensorForce

class pylith::sources::RickerWavelet : public pylith::sources::SourceTimeFunctionMomentTensorForce {
    friend class TestRickerWavelet; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    RickerWavelet(void);

    /// Destructor.
    ~RickerWavelet(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::sources::AuxiliaryFactoryMomentTensorForce* getAuxiliaryFactory(void);

    /** Add source time subfields to auxiliary field.
     *
     * @param[inout] auxiliaryField Auxiliary field.
     */
    void addAuxiliarySubfields(void);

    /** Get g1v kernel for residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return residual kernel for g1v.
     */
    PetscPointFunc getKernelg1v_explicit(const spatialdata::geocoords::CoordSys* coordsys) const;

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::sources::AuxiliaryFactorySourceTime* _auxiliaryFactory; ///< Factory for creating auxiliary subfields.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    RickerWavelet(const RickerWavelet&); ///< Not implemented.
    const RickerWavelet& operator=(const RickerWavelet&); /// Not implemented.

}; // class RickerWavelet

#endif // pylith_materials_rickerwavelet_hh

// End of file
