// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/**
 * @file modulesrc/topology/FieldBase.hh
 *
 * @brief Python interface to C++ abstract base class FieldBase.
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
      }; // omainEnum
      
      // PUBLIC MEMBERS /////////////////////////////////////////////////
    public :

      /// Default constructor.
      FieldBase(void);

      /// Destructor.
      ~FieldBase(void);
      
      /** Set name of field.
       *
       * @param value Name of field.
       */
      void name(const char* value);
      
      /** Get name of field.
       *
       * @returns Name of field.
       */
      const char* name(void) const;
      
      /** Set vector field type
       *
       * @param value Type of vector field.
       */
      void vectorFieldType(const VectorFieldEnum value);
      
      /** Get vector field type
       *
       * @returns Type of vector field.
       */
      VectorFieldEnum vectorFieldType(void) const;
      
      /** Set scale for dimensionalizing field.
       *
       * @param value Scale associated with field.
       */
      void scale(const double value);
      
      /** Get scale for dimensionalizing field.
       *
       * @returns Scale associated with field.
       */
      double scale(void) const;
      
      /** Set flag indicating whether it is okay to dimensionalize field.
       *
       * @param value True if it is okay to dimensionalize field.
       */
      void addDimensionOkay(const bool value);
      
      /** Set flag indicating whether it is okay to dimensionalize field.
       *
       * @param value True if it is okay to dimensionalize field.
       */
      bool addDimensionOkay(void) const;
      
    }; // FieldBase
    
  } // topology
} // pylith


// End of file 
