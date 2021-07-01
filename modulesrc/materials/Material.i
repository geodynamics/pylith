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
        class Material : public pylith::problems::Physics {
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

            /** Set value of label material-id used to identify material cells.
             *
             * @param value Material identifier
             */
            void setMaterialId(const int value);

            /** Get value of label material-id used to identify material cells.
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
	    pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution);

        }; // class Material

    } // materials
} // pylith

// End of file
