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
 * @file modulesrc/meshio/DataWriterVTK.i
 *
 * @brief Python interface to C++ DataWriterVTK object.
 */

namespace pylith {
    namespace meshio {
        class pylith::meshio::DataWriterVTK : public DataWriter {
            // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

            /// Constructor
            DataWriterVTK(void);

            /// Destructor
            ~DataWriterVTK(void);

            /** Make copy of this object.
             *
             * @returns Copy of this.
             */
            DataWriter* clone(void) const;

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set filename for VTK file.
             *
             * @param filename Name of VTK file.
             */
            void filename(const char* filename);

            /** Set time format for time stamp in name of VTK file.
             *
             * @param format C style time format for filename.
             */
            void timeFormat(const char* format);

            /** Set value used to normalize time stamp in name of VTK file.
             *
             * Time stamp is divided by this value (time in seconds).
             *
             * @param value Value (time in seconds) used to normalize time stamp in
             * filename.
             */
            void timeConstant(const PylithScalar value);

            /** Set precision of floating point values in output.
             *
             * @param value Precision for floating point values.
             */
            void precision(const int value);

            /** Prepare for writing files.
             *
             * @param mesh Finite-element mesh.
             * @param isInfo True if only writing info values.
             */
            void open(const pylith::topology::Mesh& mesh,
                      const bool isInfo);

            /// Close output files.
            void close(void);

            /** Prepare file for data at a new time step.
             *
             * @param t Time stamp for new data
             * @param mesh Finite-element mesh.
             */
            void openTimeStep(const PylithScalar t,
                              const pylith::topology::Mesh& mesh);

            /// Cleanup after writing data for a time step.
            void closeTimeStep(void);

            /** Write field over vertices to file.
             *
             * @param[in] t Time associated with field.
             * @param[in] subfield Subfield with basis order 1.
             */
            void writeVertexField(const PylithScalar t,
                                  const pylith::meshio::OutputSubfield& field);

            /** Write field over cells to file.
             *
             * @param[in] t Time associated with field.
             * @param[in] subfield Subfield with basis order 0.
             */
            void writeCellField(const PylithScalar t,
                                const pylith::meshio::OutputSubfield& subfield);

        }; // DataWriterVTK

    } // meshio
} // pylith

// End of file
