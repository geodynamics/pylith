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
        class FieldBase {
            // PUBLIC ENUMS ///////////////////////////////////////////////////
          public:
            enum VectorFieldEnum {
                SCALAR = 0, ///< Scalar.
                VECTOR = 1, ///< Vector.
                TENSOR = 2, ///< Tensor.
                OTHER = 3, ///< Not a scalar, vector, or tensor.
                MULTI_SCALAR = 4, ///< Scalar at multiple points.
                MULTI_VECTOR = 5, ///< Vector at multiple points.
                MULTI_TENSOR = 6, ///< Tensor at multiple points.
                MULTI_OTHER = 7, ///< Not a scalar, vector, or tensor at multiple points.
            }; // VectorFieldEnum

            enum SpaceEnum {
                POLYNOMIAL_SPACE = 0, ///< Polynomial finite-element space.
                POINT_SPACE = 1, ///< Point finite-element space.
            }; // SpaceEnum

            enum CellBasis {
                SIMPLEX_BASIS = 1, ///< Simplex basis functions.
                TENSOR_BASIS = 2, ///< Tensor product basis functions.
                DEFAULT_BASIS = 10, ///< Use default for cell type.
            }; // CellBasis

            // PUBLIC TYPEDEF /////////////////////////////////////////////////
          public:
            /// Function prototype for validator functions.
            typedef const char* (*validatorfn_type)(const PylithReal);

            // PUBLIC MEMBERS /////////////////////////////////////////////////
          public:
            FieldBase(void); ///< Default constructor.
            ~FieldBase(void); ///< Default destructor.

        }; // FieldBase

    } // namespace topology
} // namespace pylith


// End of file
