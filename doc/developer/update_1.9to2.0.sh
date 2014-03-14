#!/bin/bash

cat<<EOF > rmtemplates.sed
s/CellFilterAvgMesh/CellFilterAvg/g
s/CellFilterAvgSubMesh/CellFilterAvg/g
s/DataWriterVTKMesh/DataWriterVTK/g
s/DataWriterVTKSubMesh/DataWriterVTK/g
s/DataWriterVTKSubSubMesh/DataWriterVTK/g
s/DataWriterHDF5Mesh/DataWriterHDF5/g
s/DataWriterHDF5SubMesh/DataWriterHDF5/g
s/DataWriterHDF5SubSubMesh/DataWriterHDF5/g
s/DataWriterHDF5ExtMesh/DataWriterHDF5Ext/g
s/DataWriterHDF5ExtSubMesh/DataWriterHDF5Ext/g
s/DataWriterHDF5ExtSubSubMesh/DataWriterHDF5Ext/g
EOF

for f in `find . -name "*.cfg"`; do
  sed -f rmtemplates.sed $f > tmp && mv -f tmp $f
done
rm rmtemplates.sed
