// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/**
 * @file modulesrc/problems/Problem.hh
 *
 * @brief Python interface to C++ Problem.
 */

namespace pylith {
    namespace problems {
        class InitialCondition : public pylith::utils::PyreComponent {
            // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Constructor
            InitialCondition(void);

            /// Destructor
            virtual ~InitialCondition(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

	    /** Set fields for initial condition.
	     *
	     * @param[in] subfields Array of names of solution subfields.
	     * @param[in] numSubfields Number of subfields.
	     */
	    %apply(const char* const* string_list, const int list_len){
		(const char* subfields[], const int numSubfields)
	    };
	    void setSubfields(const char* subfields[],
			   const int numSubfields);
	    %clear(const char* subfields[], const int numSubfields);

            /** Set solver type.
             *
             * @param[out] solution Solution field.
             * @param[in] normalizer Nondimensionalization.
             */
            virtual
            void setValues(pylith::topology::Field* solution,
                           const spatialdata::units::Nondimensional& normalizer) = 0;

        }; // InitialCondition

    } // problems
} // pylith

// End of file
