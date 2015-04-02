// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/topology/Field.hh
 *
 * @brief Python interface to C++ FieldBase object.
 */

namespace pylith {
  namespace topology {

    class FieldBase
    { // FieldBase

      // PUBLIC ENUMS ///////////////////////////////////////////////////
    public :

      enum VectorFieldEnum {
	SCALAR=0, ///< Scalar.
	VECTOR=1, ///< Vector.
	TENSOR=2, ///< Tensor.
	OTHER=3, ///< Not a scalar, vector, or tensor.
	MULTI_SCALAR=4, ///< Scalar at multiple points.
	MULTI_VECTOR=5, ///< Vector at multiple points.
	MULTI_TENSOR=6, ///< Tensor at multiple points.
	MULTI_OTHER=7, ///< Not a scalar, vector, or tensor at multiple points.
      }; // VectorFieldEnum
      
      enum DomainEnum {
	VERTICES_FIELD=0, ///< FieldBase over vertices.
	CELLS_FIELD=1, ///< FieldBase over cells.
      }; // DomainEnum

      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      FieldBase(void); ///< Default constructor.
      ~FieldBase(void); ///< Default destructor.

    }; // FieldBase

  } // topology
} // pylith


// End of file
