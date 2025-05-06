// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file modulesrc/materials/Material.i
 *
 * Python interface to C++ abstract base class Material.
 */

namespace pylith {
    namespace materials {
        class Material: public pylith::problems::Physics {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
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

            /** Set gravity field.
             *
             * @param g Gravity field.
             */
            void setGravityField(spatialdata::spatialdb::GravityField* const g);

            /** Create constraint and set kernels.
             *
             * @param[in] solution Solution field.
             * @returns Constraint if applicable, otherwise NULL.
             */
            virtual
            std::vector < pylith::feassemble::Constraint* > createConstraints(const pylith::topology::Field& solution);

            /** Get default PETSc solver options appropriate for material.
             *
             * @param[in] isParallel True if running in parallel, False if running in serial.
             * @param[in] hasFault True if problem has fault, False otherwise.
             * @returns PETSc solver options.
             */
            virtual
            pylith::utils::PetscOptions* getSolverDefaults(const bool isParallel,
                                                           const bool hasFault) const;

        }; // class Material

    } // materials
} // pylith

// End of file
