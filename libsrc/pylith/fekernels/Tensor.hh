/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University at Buffalo
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2022 University of California, Davis
 *
 * See LICENSE.md for license information.
 *
 * ----------------------------------------------------------------------
 */

/** @file libsrc/fekernels/Tensor.hh
 */

#if !defined(pylith_fekernels_tensor_hh)
#define pylith_fekernels_tensor_hh

#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()
#include <iostream>

/* 2D vector: xx:0 yy:1 zz:2 xy:3
 * 2D tensor: xx:0 xy:1 yx:2 yy:3
 * 3D vector: xx:0 yy:1 zz:2 xy:3 yz:4 xz:5
 * 3D tensor: xx:0 xy:1 xz:2 yx:3 yy:4 yz:5 zx:6 zy:7 zz:7
 */

class pylith::fekernels::Tensor {
public:

    PylithReal xx;
    PylithReal yy;
    PylithReal zz;
    PylithReal xy;
    PylithReal yz;
    PylithReal xz;

    Tensor(void) :
        xx(0.0),
        yy(0.0),
        zz(0.0),
        xy(0.0),
        yz(0.0),
        xz(0.0) {}


    static const TensorOps ops2D;
    static const TensorOps ops3D;
};

class pylith::fekernels::TensorOps {
    friend class Tensor;
public:

    typedef void (*fromfn_type)(const PylithReal[],
                                Tensor*);
    typedef void (*tofn_type)(const Tensor&,
                              PylithReal[]);

    fromfn_type fromVector;
    fromfn_type fromTensor;
    tofn_type toVector;
    tofn_type toTensor;

    const int vectorSize;

    static inline
    TensorOps _create2D(void) {
        return TensorOps(_fromVector2D,
                         _fromTensor2D,
                         _toVector2D,
                         _toTensor2D,
                         4);
    }

    static inline
    TensorOps _create3D(void) {
        return TensorOps(_fromVector3D,
                         _fromTensor3D,
                         _toVector3D,
                         _toTensor3D,
                         6);
    }

    static inline
    void print(const char* name,
               const Tensor& tensor) {
        std::cout << name << ": "
                  << "xx=" << tensor.xx << ", "
                  << "yy=" << tensor.yy << ", "
                  << "zz="  << tensor.zz << ", "
                  << "xy=" << tensor.xy << ", "
                  << "yz=" << tensor.yz << ", "
                  << "xz=" << tensor.xz << std::endl;
    }

private:

    inline
    TensorOps(fromfn_type fromV,
              fromfn_type fromT,
              tofn_type toV,
              tofn_type toT,
              const int vSize) :
        fromVector(fromV),
        fromTensor(fromT),
        toVector(toV),
        toTensor(toT),
        vectorSize(vSize) {}


    static inline
    void _fromVector2D(const PylithReal vector2D[],
                       Tensor* tensor) {
        assert(vector2D);
        assert(tensor);
        tensor->xx = vector2D[0];
        tensor->yy = vector2D[1];
        tensor->zz = vector2D[2];
        tensor->xy = vector2D[3];
        tensor->yz = 0.0;
        tensor->xz = 0.0;
    }

    static inline
    void _fromTensor2D(const PylithReal tensor2D[],
                       Tensor* tensor) {
        assert(tensor);
        assert(tensor2D);
        tensor->xx = tensor2D[0];
        tensor->yy = tensor2D[3];
        tensor->zz = 0.0;
        tensor->xy = tensor2D[1];
        tensor->yz = 0.0;
        tensor->xz = 0.0;
    }

    static inline
    void _toVector2D(const Tensor& tensor,
                     PylithReal vector2D[]) {
        assert(vector2D);
        vector2D[0] = tensor.xx;
        vector2D[1] = tensor.yy;
        vector2D[2] = tensor.zz;
        vector2D[3] = tensor.xy;
    }

    static inline
    void _toTensor2D(const Tensor& tensor,
                     PylithReal tensor2D[]) {
        assert(tensor2D);
        tensor2D[0] = tensor.xx;
        tensor2D[1] = tensor.xy;
        tensor2D[2] = tensor.xy;
        tensor2D[3] = tensor.yy;
    }

    static inline
    void _fromVector3D(const PylithReal vector3D[],
                       Tensor* tensor) {
        assert(vector3D);
        assert(tensor);
        tensor->xx = vector3D[0];
        tensor->yy = vector3D[1];
        tensor->zz = vector3D[2];
        tensor->xy = vector3D[3];
        tensor->yz = vector3D[4];
        tensor->xz = vector3D[5];
    }

    static inline
    void _fromTensor3D(const PylithReal tensor3D[],
                       Tensor* tensor) {
        assert(tensor);
        assert(tensor3D);
        tensor->xx = tensor3D[0];
        tensor->yy = tensor3D[4];
        tensor->zz = tensor3D[8];
        tensor->xy = tensor3D[1];
        tensor->yz = tensor3D[5];
        tensor->xz = tensor3D[2];
    }

    static inline
    void _toVector3D(const Tensor& tensor,
                     PylithReal vector3D[]) {
        assert(vector3D);
        vector3D[0] = tensor.xx;
        vector3D[1] = tensor.yy;
        vector3D[2] = tensor.zz;
        vector3D[3] = tensor.xy;
        vector3D[4] = tensor.yz;
        vector3D[5] = tensor.xz;
    }

    static inline
    void _toTensor3D(const Tensor& tensor,
                     PylithReal tensor3D[]) {
        assert(tensor3D);
        tensor3D[0] = tensor.xx;
        tensor3D[1] = tensor.xy;
        tensor3D[2] = tensor.xz;
        tensor3D[3] = tensor.xy;
        tensor3D[4] = tensor.yy;
        tensor3D[5] = tensor.yz;
        tensor3D[6] = tensor.xz;
        tensor3D[7] = tensor.yz;
        tensor3D[8] = tensor.zz;
    }

}; // TensorOps

#endif // pylith_fekernels_tensor_hh

// End of file
