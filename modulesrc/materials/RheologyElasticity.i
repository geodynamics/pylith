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

/** @file modulesrc/materials/RheologyElasticity.i
 *
 * Python interface to C++ abstract base class RheologyElasticity.
 */

namespace pylith {
    namespace materials {
        class RheologyElasticity : public pylith::utils::PyreComponent {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            RheologyElasticity(void);

            /// Destructor.
            virtual ~RheologyElasticity(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            virtual
            pylith::materials::AuxiliaryFactoryElasticity* getAuxiliaryFactory(void) = 0;

            /// Add rheology subfields to auxiliary field.
            virtual
            void addAuxiliarySubfields(void) = 0;

            /** Get stress kernel for RHS residual, G(t,s).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return RHS residual kernel for stress.
             */
            virtual
            PetscPointFunc getKernelResidualStress(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

            /** Get elastic constants kernel for RHS Jacobian G(t,s).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return RHS Jacobian kernel for elastic constants.
             */
            virtual
            PetscPointJac getKernelJacobianElasticConstants(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

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

        // class RheologyElasticity

    } // materials
} // pylith

// End of file
