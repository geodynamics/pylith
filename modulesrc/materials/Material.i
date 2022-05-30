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

            /** Set descriptive label for material.
             *
             * @param value Label of material.
             */
            void setDescription(const char* value);

            /** Get descruptive label of material.
             *
             * @returns Label of material
             */
            const char* getDescription(void) const;

            /** Set name of label marking material.
             *
             * @param[in] value Name of label for material (from mesh generator).
             */
            void setLabelName(const char* value);

            /** Get name of label marking material.
             *
             * @returns Name of label for material (from mesh generator).
             */
            const char* getLabelName(void) const;

            /** Set value of label marking material.
             *
             * @param[in] value Value of label for material (from mesh generator).
             */
            void setLabelValue(const int value);

            /** Get value of label marking material.
             *
             * @returns Value of label for material (from mesh generator).
             */
            int getLabelValue(void) const;

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
