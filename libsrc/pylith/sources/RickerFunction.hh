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

/** @file libsrc/sources/RickerFunction.hh
 *
 * @brief C++ class for Ricker source time function.
 */

#if !defined(pylith_materials_rickerfunction_hh)
#define pylith_materials_rickerfunction_hh

#include "sourcesfwd.hh" // forward declarations

#include "pylith/sources/SourceTimeFunctionPointForce.hh" // ISA SourceTimeFunctionPointForce

class pylith::sources::RickerFunction : public pylith::sources::SourceTimeFunctionPointForce {
    friend class TestRickerFunction; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    RickerFunction(void);

    /// Destructor.
    ~RickerFunction(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::sources::AuxiliaryFactoryPointForce* getAuxiliaryFactory(void);

    /** Add source time subfields to auxiliary field.
     *
     * @param[inout] auxiliaryField Auxiliary field.
     */
    void addAuxiliarySubfields(void);

    /** Get g1v kernel for residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual kernel for stress.
     */
    PetscPointFunc getKernelg1v_explicit(const spatialdata::geocoords::CoordSys* coordsys) const;

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::sources::AuxiliaryFactoryRickerFunction* _auxiliaryFactory; ///< Factory for creating auxiliary subfields.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    RickerFunction(const RickerFunction&); ///< Not implemented.
    const RickerFunction& operator=(const RickerFunction&); /// Not implemented.

}; // class RickerFunction

#endif // pylith_materials_rickerfunction_hh

// End of file