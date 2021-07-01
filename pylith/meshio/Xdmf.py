# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

# @file pythia.pyre/meshio/Xdmf.py
##
# @brief Python class for Xdmf metadata file associated with an HDF5 file.
##
# Factory: xdmf


class Field(object):
    """Python object for data associated with vertex or cell field in HDF5 file.
    """
    DOMAIN_VERTICES = "Node"
    DOMAIN_CELLS = "Cell"

    groupToDomain = {
        "vertex_fields": DOMAIN_VERTICES,
        "cell_fields": DOMAIN_CELLS,
    }
    domainToGroup = dict((v, k) for k, v in groupToDomain.items())

    def __init__(self):
        self.name = None
        self.vectorFieldType = None
        self.data = None
        self.domain = None
        return


# Xdmf class
class Xdmf(object):
    """Python class for Xdmf metadata file associated with an HDF5 file.
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self):
        """Constructor.
        """
        self.file = None
        return

    def write(self, filenameH5, filenameXdmf=None, verbose=False):
        """Write Xdmf file corresponding to given HDF5 file.
        """
        if not filenameXdmf:
            filenameXdmf = filenameH5.replace(".h5", ".xmf")
        if verbose:
            print("Generating %s..." % filenameXdmf)

        import h5py
        import os
        if not os.path.isfile(filenameH5):
            raise IOError(
                "Cannot create Xdmf file for HDF5 file '%s'. File not found." % filenameH5)
        self.h5 = h5py.File(filenameH5, "r")

        if self._getSpaceDim() == 1:
            print("WARNING: Xdmf grids are not defined for 1-D domains.\n"
                  "Skipping creation of Xdmf file for HDF5 file '%s'." % filenameH5)
            self.h5.close()
            return

        self.file = open(filenameXdmf, "w")

        # Header
        self._openXdmf()

        # Domain
        cells = self.h5["/topology/cells"]
        vertices = self.h5["/geometry/vertices"]
        self._openDomain(cells, vertices)

        timeStamps = self._getTimeStamps()
        fields = self._getFields()

        if not timeStamps is None:
            self._openTimeCollection()
            self._writeTimeStamps(timeStamps)
            for iTime, timeStamp in enumerate(timeStamps):
                self._openTimeGrid()
                self._writeGridTopology(cells)
                self._writeGridGeometry(vertices)
                for field in fields:
                    if field.vectorFieldType in ["Tensor6", "Matrix"]:
                        numComponents = field.data.shape[-1]
                        for iComponent in range(numComponents):
                            self._writeGridFieldComponent(
                                field, iTime, iComponent)
                    else:
                        self._writeGridField(field, iTime)
                self._closeTimeGrid()
            self._closeTimeCollection()
        else:
            iTime = None
            self._openTimeGrid()
            self._writeGridTopology(cells)
            self._writeGridGeometry(vertices)
            for field in fields:
                if field.vectorFieldType in ["Tensor6", "Matrix"]:
                    numComponents = field.data.shape[-1]
                    for iComponent in range(numComponents):
                        self._writeGridFieldComponent(field, iTime, iComponent)
                else:
                    self._writeGridField(field, iTime)
            self._closeTimeGrid()

        self._closeDomain()
        self._closeXdmf()
        self._close()
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _close(self):
        self.h5.close()
        self.file.close()
        return

    def _getSpaceDim(self):
        vertices = self.h5["/geometry/vertices"]
        assert(2 == len(vertices.shape))
        return vertices.shape[1]

    def _getXdmfCellType(self, cells):
        """Get type of cell.
        """
        assert(2 == len(cells.shape))
        numCells, numCorners = cells.shape
        if "cell_dim" in cells.attrs:
            cellDim = cells.attrs["cell_dim"]
        else:
            # Use space dimension as a proxy for cell dimension.
            vertices = self.h5["/geometry/vertices"]
            assert(2 == len(vertices.shape))
            cellDim = vertices.shape[1]
        if 0 == cellDim and 1 == numCorners:
            cellType = "Polyvertex"
        elif 1 == cellDim and 2 == numCorners:
            cellType = "Polyline"
        elif 2 == cellDim and 3 == numCorners:
            cellType = "Triangle"
        elif 2 == cellDim and 4 == numCorners:
            cellType = "Quadrilateral"
        elif 3 == cellDim and 4 == numCorners:
            cellType = "Tetrahedron"
        elif 3 == cellDim and 8 == numCorners:
            cellType = "Hexahedron"
        else:
            cellType = "Unknown"
            print("WARNING: Unknown cell type with %d vertices and dimension %d." % (
                numCorners, cellDim))
        return cellType

    def _getXdmfVectorFieldType(self, vectorFieldString):
        """Get Xdmf vector field type.
        """
        vtype = "Matrix"
        if vectorFieldString.decode().lower() == "scalar":
            vtype = "Scalar"
        elif vectorFieldString.decode().lower() == "vector":
            vtype = "Vector"
        elif vectorFieldString.decode().lower() == "tensor":
            vtype = "Tensor6"
        return vtype

    def _getTimeStamps(self):
        """Get time stamps if they exist, otherwise return None.
        """
        timeStamps = None
        if "time" in self.h5:
            timeStamps = self.h5["time"][:]
        return timeStamps

    def _getFields(self):
        fields = []
        for group in ["vertex_fields", "cell_fields"]:
            if group in self.h5:
                vfields = self.h5[group]
                for name, dataset in vfields.items():
                    field = Field()
                    field.name = name
                    field.data = dataset[:]
                    field.domain = Field.groupToDomain[group]
                    if "vector_field_type" in dataset.attrs:
                        field.vectorFieldType = self._getXdmfVectorFieldType(
                            dataset.attrs["vector_field_type"])
                    else:
                        print(
                            "WARNING: Field '%s/%s' dataset missing 'vector_field_type' attribute. Guessing vector field type." % (group, name,))
                        field.vectorFieldType = "Matrix"
                        if len(dataset.shape) == 2 or len(dataset.shape) == 3:
                            numComponents = dataset.shape[-1]
                            if numComponents == 1:
                                field.vectorFieldType = "Scalar"
                            elif numComponents == 3:
                                field.vectorFieldType = "Vector"
                    fields.append(field)
        return fields

    def _openXdmf(self):
        """Write header and create Xdmf element.
        """
        import os

        filenameH5 = os.path.split(self.h5.filename)[-1]
        self.file.write(
            "<?xml version=\"1.0\" ?>\n"
            "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [\n"
            "<!ENTITY HeavyData \"%s\">\n"
            "]>\n"
            "\n"
            "<Xdmf>\n"
            % (filenameH5,)
        )
        return

    def _closeXdmf(self):
        """Close Xdmf element.
        """
        self.file.write(
            "</Xdmf>\n"
        )
        return

    def _openDomain(self, cells, vertices):
        self.file.write("  <Domain Name=\"domain\">\n")

        # Cells
        assert(2 == len(cells.shape))
        numCells, numCorners = cells.shape
        self.file.write(
            "    <DataItem Name=\"cells\" ItemType=\"Uniform\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"%d %d\">\n"
            "      &HeavyData;:/topology/cells\n"
            "    </DataItem>\n"
            % (numCells, numCorners,)
        )

        # Vertices
        assert(2 == len(vertices.shape))
        numVertices, spaceDim = vertices.shape
        if 3 == spaceDim:
            self.file.write(
                "    <DataItem Name=\"vertices\" ItemType=\"Uniform\" Format=\"HDF\" Dimensions=\"%d %d\">\n"
                "      &HeavyData;:/geometry/vertices\n"
                "    </DataItem>\n"
                % (numVertices, spaceDim),
            )
        elif 2 == spaceDim:
            # Form vector with 3 components using x and y components
            # and then a fake z-component by multiplying the
            # x-component by zero.
            self.file.write(
                "    <DataItem Name=\"vertices\" ItemType=\"Function\" Dimensions=\"%d 3\" Function=\"JOIN($0, $1, $2)\">\n" % numVertices)
            # x component
            self.file.write(
                "      <DataItem Name=\"verticesX\" ItemType=\"Hyperslab\" Type=\"HyperSlab\" Dimensions=\"%d 1\">\n"
                "        <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
                "          0 0   1 1   %d 1\n"
                "        </DataItem>\n"
                "        <DataItem Dimensions=\"%d 1\" Format=\"HDF\">\n"
                "          &HeavyData;:/geometry/vertices\n"
                "        </DataItem>\n"
                "      </DataItem>\n"
                % (numVertices, numVertices, numVertices,)
            )

            # y component
            self.file.write(
                "      <DataItem Name=\"verticesY\" ItemType=\"Hyperslab\" Type=\"HyperSlab\" Dimensions=\"%d 1\">\n"
                "        <DataItem Dimensions=\"3 2\" Format=\"XML\">\n"
                "          0 1   1 1   %d 1\n"
                "        </DataItem>\n"
                "        <DataItem Dimensions=\"%d 1\" Format=\"HDF\">\n"
                "          &HeavyData;:/geometry/vertices\n"
                "        </DataItem>\n"
                "      </DataItem>\n"
                % (numVertices, numVertices, numVertices,))

            # z component
            self.file.write(
                "      <DataItem Name=\"verticesZ\" ItemType=\"Function\" Dimensions=\"%d 1\" Function=\"0*$0\">\n"
                "        <DataItem Reference=\"XML\">\n"
                "          /Xdmf/Domain/DataItem[@Name=\"vertices\"]/DataItem[@Name=\"verticesX\"]\n"
                "        </DataItem>\n"
                "      </DataItem>\n" % (numVertices,)
            )

            self.file.write("    </DataItem>\n")
        else:
            self._close()
            raise ValueError(
                "Unexpected spatial dimension %d when writing domain vertices." % spaceDim)
        return

    def _closeDomain(self):
        """Close domain element.
        """
        self.file.write(
            "  </Domain>\n"
        )
        return

    def _openTimeCollection(self):
        """Create Grid element for collection of time grids.
        """
        self.file.write(
            "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n"
        )
        return

    def _closeTimeCollection(self):
        """Close Grid element for collection of time grids.
        """
        self.file.write(
            "    </Grid>\n"
        )
        return

    def _writeTimeStamps(self, tstamps):
        """Write time stamps.
        """
        self.file.write(
            "      <Time TimeType=\"List\">\n"
            "        <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"%d\">\n"
            "        "
            % (tstamps.shape[0],)
        )
        for t in tstamps:
            self.file.write("  %16.8e" % t)
        self.file.write(
            "\n"
            "        </DataItem>\n"
            "      </Time>\n"
        )
        return

    def _openTimeGrid(self):
        """Create Grid element for a single time step.
        """
        self.file.write(
            "      <Grid Name=\"domain\" GridType=\"Uniform\">\n"
        )
        return

    def _closeTimeGrid(self):
        """Close Grid element for a single time step.
        """
        self.file.write(
            "      </Grid>\n"
        )
        return

    def _writeGridTopology(self, cells):
        """Write topology information for current grid.
        """
        cellType = self._getXdmfCellType(cells)
        assert(2 == len(cells.shape))
        numCells = cells.shape[0]

        self.file.write(
            "        <Topology TopologyType=\"%s\" NumberOfElements=\"%d\">\n"
            "          <DataItem Reference=\"XML\">\n"
            "            /Xdmf/Domain/DataItem[@Name=\"cells\"]\n"
            "          </DataItem>\n"
            "        </Topology>\n"
            % (cellType, numCells,)
        )
        return

    def _writeGridGeometry(self, vertices):
        """Write vertices information for current grid.
        """
        self.file.write(
            "        <Geometry GeometryType=\"XYZ\">\n"
            "          <DataItem Reference=\"XML\">\n"
            "            /Xdmf/Domain/DataItem[@Name=\"vertices\"]\n"
            "          </DataItem>\n"
            "        </Geometry>\n"
        )
        return

    def _writeGridFieldComponent(self, field, iTime, iComponent):
        """Write single component of field for current time step.
        """
        if field.vectorFieldType == "Vector":
            components = ["_x", "_y", "_z"]
            componentName = field.name + components[iComponent]
        elif field.vectorFieldType == "Tensor6":
            spaceDim = self._getSpaceDim()
            if spaceDim == 2:
                components = ["_xx", "_yy", "_xy"]
            elif spaceDim == 3:
                components = ["_xx", "_yy", "_zz", "_xy", "_yz", "_xz"]
            else:
                self._close()
                raise ValueError(
                    "Unexpected spatial dimension %d for field component names." % spaceDim)
            componentName = field.name + components[iComponent]
        elif field.vectorFieldType == "Matrix":
            componentName = field.name + "_%d" % iComponent
        else:
            self._close()
            raise ValueError(
                "Unexpected vector field type '%s' for field component names." % field.vectorFieldType)
        h5Name = "/" + Field.domainToGroup[field.domain] + "/" + field.name
        if iTime is None:
            assert(2 == len(field.data.shape))
            numPoints, numComponents = field.data.shape
            numTimeSteps = 1
        else:
            assert(3 == len(field.data.shape))
            numTimeSteps, numPoints, numComponents = field.data.shape

        self.file.write(
            "        <Attribute Name=\"%(componentName)s\" Type=\"Scalar\" Center=\"%(domain)s\">\n"
            "          <DataItem ItemType=\"HyperSlab\" Dimensions=\"1 %(numPoints)d 1\" Type=\"HyperSlab\">\n"
            "            <DataItem Dimensions=\"3 3\" Format=\"XML\">\n"
            "              %(iTime)d 0 %(iComponent)d    1 1 1    1 %(numPoints)d 1\n"
            "            </DataItem>\n"
            "            <DataItem DataType=\"Float\" Precision=\"8\" Dimensions=\"%(numTimeSteps)d %(numPoints)d %(numComponents)d\" Format=\"HDF\">\n"
            "              &HeavyData;:%(h5Name)s\n"
            "            </DataItem>\n"
            "          </DataItem>\n"
            "        </Attribute>\n"
            % {"componentName": componentName,
               "domain": field.domain,
               "numPoints": numPoints,
               "iTime": iTime,
               "iComponent": iComponent,
               "numTimeSteps": numTimeSteps,
               "numComponents": numComponents,
               "h5Name": h5Name,
               }
        )
        return

    def _writeGridField(self, field, iTime):
        """Write field for current time step.
        """
        gridRef = "/Xdmf/Domain/Grid" if iTime is None else "/Xdmf/Domain/Grid/Grid[1]"
        self.file.write(
            "        <Attribute Name=\"%s\" Type=\"%s\" Center=\"%s\">\n"
            % (field.name, field.vectorFieldType, field.domain,)
        )
        h5Name = "/" + Field.domainToGroup[field.domain] + "/" + field.name
        iStep = iTime
        if iTime is None:
            iStep = 0
            if 2 == len(field.data.shape):
                numPoints, numComponents = field.data.shape
                numTimeSteps = 1
            elif 3 == len(field.data.shape):
                numTimeSteps, numPoints, numComponents = field.data.shape
            else:
                raise ValueError(
                    "Unexpected shape for dataset '%s'." % field.name)
        else:
            assert(3 == len(field.data.shape))
            numTimeSteps, numPoints, numComponents = field.data.shape

        if 2 == self._getSpaceDim() and field.vectorFieldType == "Vector":

            self.file.write(
                "          <DataItem ItemType=\"Function\" Dimensions=\"%d 3\" Function=\"JOIN($0, $1, $2)\">\n"
                % (numPoints,)
            )
            # x component
            self.file.write(
                "            <DataItem ItemType=\"HyperSlab\" Dimensions=\"%(numPoints)d 1\" Type=\"HyperSlab\">\n"
                "              <DataItem Dimensions=\"3 3\" Format=\"XML\">\n"
                "                %(iStep)d 0 0    1 1 1    1 %(numPoints)d 1\n"
                "              </DataItem>\n"
                "              <DataItem DataType=\"Float\" Precision=\"8\" Dimensions=\"%(numTimeSteps)d %(numPoints)d %(numComponents)d\" Format=\"HDF\">\n"
                "                &HeavyData;:%(h5Name)s\n"
                "              </DataItem>\n"
                "            </DataItem>\n"
                % {"numTimeSteps": numTimeSteps, "numPoints": numPoints, "iStep": iStep, "numComponents": numComponents, "h5Name": h5Name}
            )

            # y component
            self.file.write(
                "            <DataItem ItemType=\"HyperSlab\" Dimensions=\"%(numPoints)d 1\" Type=\"HyperSlab\">\n"
                "              <DataItem Dimensions=\"3 3\" Format=\"XML\">\n"
                "                %(iStep)d 0 1    1 1 1    1 %(numPoints)d 1\n"
                "              </DataItem>\n"
                "              <DataItem DataType=\"Float\" Precision=\"8\" Dimensions=\"%(numTimeSteps)d %(numPoints)d %(numComponents)d\" Format=\"HDF\">\n"
                "                &HeavyData;:%(h5Name)s\n"
                "              </DataItem>\n"
                "            </DataItem>\n"
                % {"numTimeSteps": numTimeSteps, "numPoints": numPoints, "iStep": iStep, "numComponents": numComponents, "h5Name": h5Name}
            )

            # z component
            self.file.write(
                "            <DataItem ItemType=\"Function\" Dimensions=\"%(numPoints)d 1\" Function=\"0*$0\">\n"
                "              <DataItem Reference=\"XML\">\n"
                "                %(gridRef)s/Attribute[@Name=\"%(name)s\"]/DataItem[1]/DataItem[1]\n"
                "              </DataItem>\n"
                "            </DataItem>\n"
                % {"numPoints": numPoints, "name": field.name, "gridRef": gridRef}
            )

            # close
            self.file.write(
                "          </DataItem>\n"
                "        </Attribute>\n"
            )

        else:
            self.file.write(

                "          <DataItem ItemType=\"HyperSlab\" Dimensions=\"1 %(numPoints)d %(numComponents)d\" Type=\"HyperSlab\">\n"
                "            <DataItem Dimensions=\"3 3\" Format=\"XML\">\n"
                "              %(iStep)d 0 0    1 1 1    1 %(numPoints)d %(numComponents)d\n"
                "            </DataItem>\n"
                "            <DataItem DataType=\"Float\" Precision=\"8\" Dimensions=\"%(numTimeSteps)d %(numPoints)d %(numComponents)d\" Format=\"HDF\">\n"
                "              &HeavyData;:%(h5Name)s\n"
                "            </DataItem>\n"
                "          </DataItem>\n"
                "        </Attribute>\n"
                % {"numTimeSteps": numTimeSteps, "numPoints": numPoints, "iStep": iStep, "numComponents": numComponents, "h5Name": h5Name}
            )

            return


# End of file
