// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


namespace pylith {
    namespace meshio {
        class OutputObserver : public pylith::utils::PyreComponent {
            // PUBLIC METHODS ///////////////////////////////////////////////////////
          public:
            /// Constructor
            OutputObserver(void);

            /// Destructor
            virtual ~OutputObserver(void);

            /// Deallocate PETSc and local data structures.
            virtual void deallocate(void);

            /** Set output trigger for how often to write output.
             *
             * @param[in] trigger Output trigger.
             */
            void setTrigger(pylith::meshio::OutputTrigger* const trigger);

            /** Get trigger for how often to write otuput.
             *
             * @returns Output trigger.
             */
            const pylith::meshio::OutputTrigger* getTrigger(void) const;

            /** Set writer to write data to file.
             *
             * @param[in] datawriter Writer for data.
             */
            void setWriter(pylith::meshio::DataWriter* const writer);

            /** Set basis order for output.
             *
             * @param[in] value Basis order for output.
             */
            void setOutputBasisOrder(const int value);

            /** Set number of mesh refinement levels for output.
             *
             * @param[in] value Number of mesh refinement levels for output.
             */
            void setRefineLevels(const int value);

            /** Set time scale.
             *
             * @param[in] value Time scale for dimensionalizing time.
             */
            void setTimeScale(const PylithReal value);

        }; // Observer

    } // namespace meshio
} // namespace pylith


// End of file
