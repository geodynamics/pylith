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

/** @file libsrc/sources/SquareWavelet.hh
 *
 * @brief C++ class for Square source time function.
 */

#if !defined(pylith_materials_squarewavelet_hh)
#define pylith_materials_squarewavelet_hh

#include "sourcesfwd.hh" // forward declarations

#include "pylith/sources/SourceTimeFunctionMomentTensorForce.hh" // ISA SourceTimeFunctionMomentTensorForce

class pylith::sources::SquareWavelet : public pylith::sources::SourceTimeFunctionMomentTensorForce {
    friend class TestSquareWavelet; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    SquareWavelet(void);

    /// Destructor.
    ~SquareWavelet(void);

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

    SquareWavelet(const SquareWavelet&); ///< Not implemented.
    const SquareWavelet& operator=(const SquareWavelet&); /// Not implemented.

}; // class SquareWavelet

#endif // pylith_materials_squarewavelet_hh

// End of file
