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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/materials/IsotropicLinearPoroelasticity.i
 *
 * Python interface to C++ IsotropicLinearPoroelasticity.
 */

namespace pylith {
    namespace materials {
        class IsotropicLinearPoroelasticity : public pylith::materials::RheologyPoroelasticity {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    IsotropicLinearPoroelasticity(void);

    /// Destructor.
    ~IsotropicLinearPoroelasticity(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Include reference stress/strain?
     *
     * @param value Flag indicating to include reference stress/strain.
     */
    void useReferenceState(const bool value);

    /** Use reference stress and strain in computation of stress and
     * strain?
     *
     * @returns True if using reference stress and strain, false otherwise.
     */
    bool useReferenceState(void) const;

    /** Include tensor permeability?
     *
     * @param value Flag indicating to include tensor permeability.
     */
    void useTensorPermeability(const bool value);

    /** Use full tensor permeability?
     *
     * @returns True if using full tensor permeability, false otherwise.
     */
    bool useTensorPermeability(void) const;

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::materials::AuxiliaryFactoryPoroelastic* getAuxiliaryFactory(void);

    /** Add rheology subfields to auxiliary field.
     *
     * @param[inout] auxiliaryField Auxiliary field.
     */
    void addAuxiliarySubfields(void);

    // ============================ Either Side ====================================
    // ---------------------------------------------------------------------------------------------------------------------
    // Get stress kernel for RHS residual, G(t,s).
    PetscPointFunc getKernelResidualStress(const spatialdata::geocoords::CoordSys* coordsys, const bool _useInertia) const;

    /** Get pressure kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for Darcy velocity.
     */
    PetscPointFunc getKernelDarcy(const spatialdata::geocoords::CoordSys* coordsys, const bool _gravityField, const bool _useInertia) const;

    // ============================= RHS ==================================== //

    // ---------------------------------------------------------------------------------------------------------------------
    // Select g0p function. Will only be used for the dynamic case.
    PetscPointFunc getKernelg0p(const spatialdata::geocoords::CoordSys* coordsys,
                                                                   const bool _useBodyForce,
                                                                   const bool _gravityField,
                                                                   const bool _useSourceDensity) const;

   // ---------------------------------------------------------------------------------------------------------------------
   /** Get pressure kernel for RHS residual, G(t,s).
   *
   * @param[in] coordsys Coordinate system.
   *
   * @return RHS residual kernel for Darcy velocity.
   */
   PetscPointFunc getKernelg1p_explicit(const spatialdata::geocoords::CoordSys* coordsys,
                                          const bool _gravityField) const;

   // ============================= LHS ==================================== //

   // ---------------------------------------------------------------------------------------------------------------------
   // Get variation in fluid content kernel for LHS residual, F(t,s,\dot{s})
   PetscPointFunc getKernelf0p_explicit(const spatialdata::geocoords::CoordSys* coordsys) const;


   // ---------------------------------------------------------------------------------------------------------------------
   // Select implicit f0p function.
   PetscPointFunc getKernelf0p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                                                  const bool _useBodyForce,
                                                                  const bool _gravityField,
                                                                  const bool _useSourceDensity) const;

  // ---------------------------------------------------------------------------------------------------------------------
  /** Get pressure kernel for LHS residual.
  *
  * @param[in] coordsys Coordinate system.
  *
  * @return RHS residual kernel for Darcy velocity.
  */
  PetscPointFunc getKernelf1p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                         const bool _gravityField) const;

  // ---------------------------------------------------------------------------------------------------------------------
  // Get poroelastic constants kernel for LHS Jacobian
  PetscPointJac getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const;

  // ---------------------------------------------------------------------------------------------------------------------
  // Get biot coefficient kernel for LHS Jacobian
  PetscPointJac getKernelJf2up(const spatialdata::geocoords::CoordSys* coordsys) const;

  // ---------------------------------------------------------------------------------------------------------------------
  // Get lambda kernel for LHS Jacobian
  PetscPointJac getKernelJf2ue(const spatialdata::geocoords::CoordSys* coordsys) const;

  // ---------------------------------------------------------------------------------------------------------------------
  // Get Specific storage kernel for LHS Jacobian F(t,s, \dot{s}).
  PetscPointJac getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const;

  // ---------------------------------------------------------------------------------------------------------------------
  // Get Darcy Conductivity kernel for LHS Jacobian
  PetscPointJac getKernelJf3pp(const spatialdata::geocoords::CoordSys* coordsys) const;

  // ---------------------------------------------------------------------------------------------------------------------
  // Get biot coefficient kernel for LHS Jacobian F(t,s, \dot{s}).
  PetscPointJac getKernelJf0pe(const spatialdata::geocoords::CoordSys* coordsys) const;

    // ============================ DERIVED FIELDS ========================== //

    /** Get stress kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing stress subfield in derived field.
     */
    PetscPointFunc getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const;

    };      // class IsotropicLinearPoroelasticity

    } // materials
} // pylith

// End of file
