// -*- swig -*-
// -*- coding: utf-8 -*-
//
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026 The PyLith authors and contributors
// SPDX-License-Identifier: MIT


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

    } // namespace topology
} // namespace pylith


// End of file
