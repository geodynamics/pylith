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

/** @file libsrc/materials/Material.hh
 *
 * @brief C++ abstract base class for Material object.
 */

#if !defined(pylith_materials_material_hh)
#define pylith_materials_material_hh

// Include directives ---------------------------------------------------
#include "materialsfwd.hh" // forward declarations

#include "pylith/feassemble/IntegratorPointwise.hh" // ISA IntegratorPointwise

#include <string> // HASA std::string

// Material -------------------------------------------------------------
/** @brief C++ abstract base class for Material object.
 *
 * Interface definition for a material. A material encapsulates both
 * the rheology as well as the governing equation.
 *
 * An individual material must abide by specific rules for the
 * interface, especially the order of the fields in the solution.
 *
 * Elasticity:
 *   + displacement, [velocity, Lagrange multipliers]
 *
 * Incompressible elasticity
 *   + displacement, pressure, [velocity, Lagrange multipliers]
 */

class pylith::materials::Material : public pylith::problems::Physics {
    friend class AuxiliaryFactory; ///< Helper for setting up auxiliary fields.

    friend class TestMaterial; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param dimension Spatial dimension associated with material.
     */
    Material(const int dimension);

    /// Destructor.
    virtual ~Material(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Get spatial dimension of material.
     *
     * @returns Spatial dimension.
     */
    int getDimension(void) const;

    /** Set identifier of material.
     *
     * @param value Material identifier
     */
    void setMaterialId(const int value);

    /** Get identifier of material.
     *
     * @returns Material identifier
     */
    int getMaterialId(void) const;

    /** Set descriptive label for material.
     *
     * @param value Label of material.
     */
    void setDescriptiveLabel(const char* value);

    /** Get descruptive label of material.
     *
     * @returns Label of material
     */
    const char* getDescriptiveLabel(void) const;

    /** Create integrator and set kernels.
     *
     * @solution Solution field.
     * @returns Integrator if applicable, otherwise NULL.
     */
    virtual
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution) = 0;

    /** Create constraint and set kernels.
     *
     * @solution Solution field.
     * @returns Constraint if applicable, otherwise NULL.
     */
    virtual
    pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution) = 0;

    /** Create auxiliary field.
     *
     * @mesh Finite-element mesh associated with physics.
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    virtual
    pylith::topology::Field* createAuxiliaryField(const pylith::topology::Mesh& mesh) = 0;

    /** Create derived field.
     *
     * @mesh Finite-element mesh associated with physics.
     * @returns Derived field if applicable, otherwise NULL.
     */
    virtual
    pylith::topology::Field* createDerivedField(const pylith::topology::Mesh& mesh) = 0;

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& solution) const = 0;

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::materials::AuxiliaryFactory* _auxiliaryFactory; ///< Factory for auxiliary subfields.

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    const int _dimension; ///< Spatial dimension of material.
    int _materialId; ///< Material identifier.
    std::string _label; ///< Label of material.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    Material(const Material&); ///< Not implemented.
    const Material& operator=(const Material&); ///< Not implemented

}; // class Material

#endif // pylith_materials_material_hh

// End of file
