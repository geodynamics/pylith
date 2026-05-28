// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file modulesrc/topology/Distributor.i
 *
 * @brief Python interface to C++ Distributor object.
 */

namespace pylith {
    namespace topology {
        class Distributor { // Distributor
                            // PUBLIC MEMBERS /////////////////////////////////////////////////
public:

            /// Constructor
            Distributor(void);

            /// Destructor
            ~Distributor(void);

            /** Set data writer.
             *
             * @param[in] writer Data writer.
             */
            void setDataWriter(pylith::meshio::DataWriter* writer);

            /** Set edge weighting.
             *
             * @param[in] useEdgeWeighting If true, weight edges in distribution.
             */
            void setUseEdgeWeighting(const bool flag);

            /** Set partitioner.
             *
             * @param[in] partitioner Name of mesh partitioner.
             */
            void setPartitioner(const char* partitioner);

            /** Distribute mesh.
             *
             * @param[in] mesh Original mesh.
             * @param[in] faults Fault interfaces.
             */
            pylith::topology::Mesh* distribute(const pylith::topology::Mesh& mesh,
                                               const std::vector<pylith::faults::FaultCohesive*>& faults) const;

        }; // Distributor

    } // topology
} // pylith

// End of file
