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

/** @file modulesrc/sources/SourceTimeFunctionMomentTensorForce.i
 *
 * Python interface to C++ abstract base class SourceTimeFunctionMomentTensorForce.
 */

namespace pylith {
    namespace sources {
        class SourceTimeFunctionMomentTensorForce : public pylith::utils::PyreComponent {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            SourceTimeFunctionMomentTensorForce(void);

            /// Destructor.
            virtual ~SourceTimeFunctionMomentTensorForce(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            virtual
            pylith::sources::AuxiliaryFactoryMomentTensorForce* getAuxiliaryFactory(void) = 0;

            /// Add sourcetimefunction subfields to auxiliary field.
            virtual
            void addAuxiliarySubfields(void) = 0;

            /** Get g1v kernel for residual, G(t,s).
            *
            * @param[in] coordsys Coordinate system.
            *
            * @return RHS residual kernel for stress.
            */
            virtual
            PetscPointFunc getKernelg1v_explicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;            

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

        // class SourceTimeFunctionMomentTensorForce

    } // sources
} // pylith

// End of file
