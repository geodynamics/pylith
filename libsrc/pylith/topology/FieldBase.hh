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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/topology/FieldBase.hh
 *
 * @brief Basic information related to a vector field over the
 * vertices or cells of a finite-element mesh.
 */

#if !defined(pylith_topology_fieldbase_hh)
#define pylith_topology_fieldbase_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

#include "pylith/utils/array.hh" // HASA string_vector

#include <string> // USES std::string
#include <map> // USES std::map

// FieldBase ------------------------------------------------------------
/** @brief Basic information related to a vector field over the
 * vertices or cells of a finite-element mesh.
 */
class pylith::topology::FieldBase { // Field
    // PUBLIC ENUMS /////////////////////////////////////////////////////////
public:

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

    enum SpaceEnum {
        POLYNOMIAL_SPACE=0, ///< Polynomial finite-element space.
        POINT_SPACE=1, ///< Point finite-element space.
    }; // SpaceEnum

    // PUBLIC TYPEDEF ///////////////////////////////////////////////////////
public:

    /// Function prototype for validator functions.
    typedef const char* (*validatorfn_type)(const PylithReal);

    // PUBLIC STRUCTS ///////////////////////////////////////////////////////
public:

    struct Description {
        std::string label; ///< Label for field.
        std::string alias; ///< User-defined alias for field.
        VectorFieldEnum vectorFieldType; ///< Type of vector field.
        pylith::string_vector componentNames; ///< Names of components.
        size_t numComponents; ///< Number of components.
        PylithReal scale; ///< Dimension scale associated with values.
        validatorfn_type validator; ///< Validator for values in field;
        bool hasHistory; ///< Has subfields with history, i.e., state variables.
        size_t historySize; ///< Number of points in time history (currently only).

        Description(void) :
            label(""),
            alias(""),
            vectorFieldType(OTHER),
            numComponents(0),
            scale(1.0),
            validator(NULL),
            hasHistory(false),
            historySize(0) {}


        Description(const char* labelValue,
                    const char* aliasValue,
                    const pylith::string_vector& componentNamesValue,
                    const size_t numComponentsValue=0,
                    const VectorFieldEnum vectorFieldTypeValue=SCALAR,
                    const PylithReal scaleValue=1.0,
                    const validatorfn_type validatorValue=NULL,
                    bool hasHistoryValue=false,
                    const size_t historySizeValue=0) :
            label(labelValue),
            alias(aliasValue),
            vectorFieldType(vectorFieldTypeValue),
            componentNames(componentNamesValue),
            numComponents(numComponentsValue),
            scale(scaleValue),
            validator(validatorValue),
            hasHistory(hasHistoryValue),
            historySize(historySizeValue) {}


    }; // Description

    struct Discretization {
        bool tensorBasis; ///< Is this a tensor polynomial basis
        int basisOrder; ///< Order of basis functions.
        int quadOrder; ///< Order of quadrature scheme.
        int dimension; ///< Dimension of point(s) for discretization.
        bool isBasisContinuous; ///< Is basis continuous?
        SpaceEnum feSpace; ///< Finite-element space.

        Discretization(const bool tensorBasisValue=false,
                       const int basisOrderValue=1,
                       const int quadOrderValue=1,
                       const int dimensionValue=-1,
                       const bool isBasisContinuousValue=true,
                       const SpaceEnum feSpaceValue=POLYNOMIAL_SPACE) :
            tensorBasis(tensorBasisValue),
            basisOrder(basisOrderValue),
            quadOrder(quadOrderValue),
            dimension(dimensionValue),
            isBasisContinuous(isBasisContinuousValue),
            feSpace(feSpaceValue)
        {}

        bool operator==(const Discretization rhs) const
        {
          if (tensorBasis        != rhs.tensorBasis)       return false;
          if (basisOrder         != rhs.basisOrder)        return false;
          if (quadOrder          != rhs.quadOrder)         return false;
          if (dimension          != rhs.dimension)         return false;
          if (isBasisContinuous  != rhs.isBasisContinuous) return false;
          if (feSpace            != rhs.feSpace)           return false;
          return true;
        }

        bool operator<(const Discretization rhs) const
        {
          if (tensorBasis        < rhs.tensorBasis)       return true;
          if (tensorBasis       == rhs.tensorBasis) {
            if (basisOrder         < rhs.basisOrder)        return true;
            if (basisOrder        == rhs.basisOrder) {
              if (quadOrder        < rhs.quadOrder)         return true;
              if (quadOrder       == rhs.quadOrder)  {
                if (dimension      < rhs.dimension)         return true;
                if (dimension     == rhs.dimension) {
                  if (isBasisContinuous  < rhs.isBasisContinuous) return true;
                  if (isBasisContinuous == rhs.isBasisContinuous) {
                    if (feSpace          < rhs.feSpace)           return true;
                  }
                }
              }
            }
          }
          return false;
        }
    }; // Discretization

    /// Mapping from field name to discretization (intended for subfields).
    typedef std::map<std::string, pylith::topology::FieldBase::Discretization> discretizations_map;

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    FieldBase(void); ///< Default constructor.
    ~FieldBase(void); ///< Default destructor.

    /** Get string associated with vector field type.
     *
     * @param value Vector field type.
     * @returns String associated with vector field type.
     */
    static
    const char*vectorFieldString(VectorFieldEnum value);

    /** Get string associated with vector field type.
     *
     * @param value String associated with vector field type.
     * @returns Vector field type.
     */
    static
    VectorFieldEnum parseVectorFieldString(const char* value);

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    FieldBase(const FieldBase&); ///< Not implemented
    const FieldBase& operator=(const FieldBase&); ///< Not implemented

}; // FieldBase

#endif // pylith_topology_fieldbase_hh

// End of file
