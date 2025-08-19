// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Discretization

#include "spatialdata/units/unitsfwd.hh" // HOLDSA Scales

class pylith::topology::FieldFactory : public pylith::utils::GenericComponent {
    friend class TestFieldFactory; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    FieldFactory(void);

    /// Destructor.
    virtual ~FieldFactory(void);

    /** Get number of subfield discretizations.
     *
     * @returns Number of subfield discretizations.
     */
    int getNumSubfields(void) const;

    /** Set discretization information for auxiliary subfield.
     *
     * @param[in] subfieldName Name of auxiliary subfield.
     * @param[in] basisOrder Polynomial order for basis.
     * @param[in] quadOrder Order of quadrature rule.
     * @param[in] dimension Dimension of points for discretization.
     * @param[in] isFaultOnly True if subfield is limited to fault degrees of freedom.
     * @param[in] cellBasis Type of basis functions to use (e.g., simplex, tensor, or default).
     * @param[in] isBasisContinuous True if basis is continuous.
     * @param[in] feSpace Finite-element space.
     */
    void setSubfieldDiscretization(const char* subfieldName,
                                   const int basisOrder,
                                   const int quadOrder,
                                   const int dimension,
                                   const bool isFaultOnly,
                                   const pylith::topology::FieldBase::CellBasis cellBasis,
                                   const pylith::topology::FieldBase::SpaceEnum feSpace,
                                   const bool isBasisContinuous);

    /** Get discretization information for subfield.
     *
     * @param[in] subfieldName Name of subfield.
     * @return Discretization information for auxiliary subfield. If
     * discretization information was not set, then use "default".
     */
    const pylith::topology::FieldBase::Discretization& getSubfieldDiscretization(const char* subfieldName) const;

    /** Initialize factory for setting up auxiliary subfields.
     *
     * @param[inout] field Auxiliary field for which subfields are to be created.
     * @param[in] scales Scales for nondimensionalization.
     * @param[in] spaceDim Spatial dimension of problem.
     * @param[in] defaultDescription Default description for new subfields.
     */
    virtual
    void initialize(pylith::topology::Field* field,
                    const spatialdata::units::Scales& scales,
                    const int spaceDim,
                    const pylith::topology::FieldBase::Description* defaultDescription=NULL);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::topology::Field* _field; ///< Auxiliary field.
    pylith::topology::FieldBase::discretizations_map _subfieldDiscretizations; ///< Discretization for each subfield.
    pylith::topology::FieldBase::Description* _defaultDescription; ///< Description for default subfield.
    spatialdata::units::Scales* _scales; ///< Scales for nondimensionalization.
    int _spaceDim; ///< Spatial dimension.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    FieldFactory(const FieldFactory &); ///< Not implemented.
    const FieldFactory& operator=(const FieldFactory&); ///< Not implemented

}; // class FieldFactory

// End of file
