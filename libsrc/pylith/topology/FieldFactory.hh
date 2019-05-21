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

/** @file libsrc/topology/FieldFactory.hh
 *
 * @brief C++ helper class for setting up subfields.
 */

#if !defined(pylith_topology_fieldfactory_hh)
#define pylith_topology_fieldfactory_hh

#include "topologyfwd.hh" // forward declarations
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Discretization

#include "spatialdata/units/unitsfwd.hh" // HOLDSA Normalizer

class pylith::topology::FieldFactory : public pylith::utils::GenericComponent {
    friend class TestFieldFactory; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    FieldFactory(void);

    /// Destructor.
    virtual ~FieldFactory(void);

    /** Set discretization information for auxiliary subfield.
     *
     * @param[in] subfieldName Name of auxiliary subfield.
     * @param[in] basisOrder Polynomial order for basis.
     * @param[in] quadOrder Order of quadrature rule.
     * @param[in] dimension Dimension of points for discretization.
     * @param[in] isBasisContinuous True if basis is continuous.
     * @param[in] feSpace Finite-element space.
     */
    void setSubfieldDiscretization(const char* subfieldName,
                                   const int basisOrder,
                                   const int quadOrder,
                                   const int dimension,
                                   const bool isBasisContinuous,
                                   const pylith::topology::FieldBase::SpaceEnum feSpace);

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
     * @param[in] normalizer Scales for nondimensionalization.
     * @param[in] spaceDim Spatial dimension of problem.
     * @param[in] defaultDescription Default description for new subfields.
     */
    virtual
    void initialize(pylith::topology::Field* field,
                    const spatialdata::units::Nondimensional& normalizer,
                    const int spaceDim,
                    const pylith::topology::FieldBase::Description* defaultDescription=NULL);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::topology::Field* _field; ///< Auxiliary field.
    pylith::topology::FieldBase::discretizations_map _subfieldDiscretizations; ///< Discretization for each subfield.
    pylith::topology::FieldBase::Description* _defaultDescription; ///< Description for default subfield.
    spatialdata::units::Nondimensional* _normalizer; ///< Scales for nondimensionalization.
    int _spaceDim; ///< Spatial dimension.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    FieldFactory(const FieldFactory &); ///< Not implemented.
    const FieldFactory& operator=(const FieldFactory&); ///< Not implemented

}; // class FieldFactory

#endif // pylith_topology_fieldfactory_hh

// End of file
