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

      /** Get stress kernel for RHS residual, G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return RHS residual kernel for stress.
       */
      virtual
      PetscPointFunc getKernelg1u(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

      /** Get stress kernel for RHS residual, G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return RHS residual kernel for stress.
       */
      virtual
      PetscPointFunc getKernelg1v(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

      /** Get pressure kernel for RHS residual, G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return RHS residual kernel for pressure.
       */
      virtual
      PetscPointFunc getKernelg1p(const spatialdata::geocoords::CoordSys* coordsys, const bool _gravityField) const = 0;

      /** Get elastic constants kernel for RHS Jacobian G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return RHS Jacobian kernel for elastic constants.
       */
      virtual
      PetscPointJac getKernelJg3uu(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

      /** Get elastic constants kernel for RHS Jacobian G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return RHS Jacobian kernel for elastic constants.
       */
      virtual
      PetscPointJac getKernelJg3vu(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

      /** Get Biot Coefficient for RHS Jacobian G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return RHS Jacobian kernel for Biot Coefficient.
       */
      virtual
      PetscPointJac getKernelJg2up(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

      /** Get Biot Coefficient for RHS Jacobian G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return RHS Jacobian kernel for Biot Coefficient.
       */
      virtual
      PetscPointJac getKernelJg2vp(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

      /** Get lambda for RHS Jacobian G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return RHS Jacobian kernel for lambda.
       */
      virtual
      PetscPointJac getKernelJg2ue(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

      /** Get stress kernel for derived field.
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return Project kernel for computing stress subfield in derived field.
       */
      virtual
      PetscPointJac getKernelJg3pp(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

      // ============================= LHS ==================================== //

      /** Get variation in fluid content for LHS residual, F(t,s,\dot{s})
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return LHS residual kernel for variation in fluid content.
       */
      virtual
      PetscPointFunc getKernelf0p(const spatialdata::geocoords::CoordSys* coordsys, const bool _useInertia) const = 0;

      /** Get kernel for RHS Jacobian G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return RHS Jacobian kernel for tshift * 1/M (Jf0pp)
       */
      virtual
      PetscPointJac getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

      /** Get kernel for RHS Jacobian G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return LHS jacobian kernel for biot coefficient.
       */
      virtual
      PetscPointJac getKernelJf0pe(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

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
                                     const spatialdata::geocoords::CoordSys* coordsys) const;

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
