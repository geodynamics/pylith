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
 * @file modulesrc/meshio/DataWriterHDF5.i
 *
 * @brief Python interface to C++ DataWriterHDF5 object.
 */

namespace pylith {
    namespace meshio {
        class pylith::meshio::DataWriterHDF5 : public DataWriter {
            // PUBLIC METHODS /////////////////////////////////////////////////
public:

            /// Constructor
            DataWriterHDF5(void);

            /// Destructor
            ~DataWriterHDF5(void);

            /** Make copy of this object.
             *
             * @returns Copy of this.
             */
            DataWriter* clone(void) const;

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set filename for HDF5 file.
             *
             * @param filename Name of HDF5 file.
             */
            void filename(const char* filename);

            /** Generate filename for HDF5 file.
             *
             * Appends _info if only writing parameters.
             *
             * :KLUDGE: We should separate generating "info" files from the
             * DataWriter interface.
             *
             * @returns String for HDF5 filename.
             */
            std::string hdf5Filename(void) const;

            /** Open output file.
             *
             * @param mesh Finite-element mesh.
             * @param isInfo True if only writing info values.
             */
            void open(const pylith::topology::Mesh& mesh,
                      const bool isInfo);

            /// Close output files.
            void close(void);

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

            /** Write dataset with names of points to file.
             *
             * @param names Array with name for each point, e.g., station name.
             * @param mesh Finite-element mesh.
             *
             * Primarily used with OutputSolnPoints.
             */
            void writePointNames(const pylith::string_vector& names,
                                 const pylith::topology::Mesh& mesh);

        }; // DataWriterHDF5

    } // meshio
} // pylith

// End of file
