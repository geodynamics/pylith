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

/** @file modulesrc/materials/RheologyPoroelasticity.i
 *
 * Python interface to C++ abstract base class RheologyPoroelasticity.
 */

namespace pylith {
    namespace materials {
        class RheologyPoroelasticity : public pylith::utils::PyreComponent {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

      /// Default constructor.
      RheologyPoroelasticity(void);

      /// Destructor.
      virtual ~RheologyPoroelasticity(void);

      /// Deallocate PETSc and local data structures.
      void deallocate(void);

      /** Get auxiliary factory associated with physics.
       *
       * @return Auxiliary factory for physics object.
       */
      virtual
      pylith::materials::AuxiliaryFactoryPoroelasticity* getAuxiliaryFactory(void) = 0;

      /// Add rheology subfields to auxiliary field.
      virtual
      void addAuxiliarySubfields(void) = 0;

      // ============================= RHS ==================================== //

      // ---------------------------------------------------------------------------------------------------------------------
      // Select g0p function. Will only be used for the dynamic case.
      virtual
      PetscPointFunc getKernelg0p(const spatialdata::geocoords::CoordSys* coordsys,
                                                                     const bool _useBodyForce,
                                                                     const bool _gravityField,
                                                                     const bool _useSourceDensity) const = 0;

     // ---------------------------------------------------------------------------------------------------------------------
     /** Get pressure kernel for RHS residual.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for Darcy velocity.
     */
     virtual
     PetscPointFunc getKernelg1p_explicit(const spatialdata::geocoords::CoordSys* coordsys,
                                            const bool _gravityField) const = 0;

     // ============================= LHS ==================================== //

     // ---------------------------------------------------------------------------------------------------------------------
     // Get variation in fluid content kernel for LHS residual, F(t,s,\dot{s})
     virtual
     PetscPointFunc getKernelf0p_explicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

     // ---------------------------------------------------------------------------------------------------------------------
     // Select implicit f0p function.
     virtual
     PetscPointFunc getKernelf0p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                                                    const bool _useBodyForce,
                                                                    const bool _gravityField,
                                                                    const bool _useSourceDensity) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    /** Get pressure kernel for LHS residual.
    *
    * @param[in] coordsys Coordinate system.
    *
    * @return LHS residual kernel for Darcy velocity.
    */
    virtual
    PetscPointFunc getKernelf1p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                         const bool _useBodyForce,
                                         const bool _gravityField,
                                         const bool _useConstantPressureSource) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get poroelastic constants kernel for LHS Jacobian
    virtual
    PetscPointJac getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get biot coefficient kernel for LHS Jacobian
    virtual
    PetscPointJac getKernelJf2up(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get lambda kernel for LHS Jacobian
    virtual
    PetscPointJac getKernelJf2ue(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get Specific storage kernel for LHS Jacobian F(t,s, \dot{s}).
    virtual
    PetscPointJac getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get Darcy Conductivity kernel for LHS Jacobian
    virtual
    PetscPointJac getKernelJf3pp(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    // Get biot coefficient kernel for LHS Jacobian F(t,s, \dot{s}).
    virtual
    PetscPointJac getKernelJf0pe(const spatialdata::geocoords::CoordSys* coordsys) const = 0;
    
    // ---------------------------------------------------------------------------------------------------------------------
    virtual
    PetscPointJac getKernelJf0ppdot(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ---------------------------------------------------------------------------------------------------------------------
    virtual
    PetscPointJac getKernelJf0pedot(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

      // ============================ DERIVED FIELDS ========================== //

      /** Get stress kernel for derived field.
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return Project kernel for computing stress subfield in derived field.
       */
      virtual
      PetscPointFunc getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

      /** Add kernels for updating state variables.
       *
       * @param[inout] kernels Array of kernels for updating state variables.
       * @param[in] coordsys Coordinate system.
       */
      virtual
      void addKernelsUpdateStateVars(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                     const spatialdata::geocoords::CoordSys* coordsys,
                                     const bool _updateFields) const;

      /** Update kernel constants.
       *
       * @param[inout] kernelConstants Array of constants used in integration kernels.
       * @param[in] dt Current time step.
       */
      virtual
      void updateKernelConstants(pylith::real_array* kernelConstants,
                                 const PylithReal dt) const;


        };

        // class RheologyPoroelasticity

    } // materials
} // pylith

// End of file
