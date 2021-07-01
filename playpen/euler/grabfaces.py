#!/usr/bin/env nemesis
#
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

# @file euler/grabfaces

# @brief Python application to grab a set of points specified in a UCD
# face description file along with the associated normals and write them to a
# file.

import math
import numpy

from pythia.pyre.applications.Script import Script as Application


class GrabFaces(Application):
    """Python application to grab a set of point coordinates and normals from a UCD
    face description file.
    """

    class Inventory(Application.Inventory):
        """Python object for managing GrabFaces facilities and properties.
        """

        # @class Inventory
        # Python object for managing GrabFaces facilities and properties.
        ##
        # \b Properties
        # @li \b ucd_face_file Filename of input UCD face description file.
        # @li \b fault_id_num ID number (material number) of fault to use.
        # @li \b point_output_file Filename of output set of points and normals.
        # @li \b node_values_list List specifying position of desired attributes in UCD face nodal attributes.
        # @li \b exclude_zero_normals Flag indicating whether to exclude points if the associated normal has zero magnitude.
        ##
        # \b Facilities
        # @li None

        import pythia.pyre.inventory

        ucdFaceFile = pythia.pyre.inventory.str(
            "ucd_face_file", default="test_face.inp")
        ucdFaceFile.meta['tip'] = "Filename of ucd file containing face descriptions."

        faultIDNum = pythia.pyre.inventory.int("fault_id_num", default=1)
        faultIDNum.meta['tip'] = "ID number (material number) of fault to use."

        pointOutputFile = pythia.pyre.inventory.str("point_output_file",
                                                    default="points.coordnorm")
        pointOutputFile.meta['tip'] = "Filename of output coordinates and normals."

        nodeValuesList = pythia.pyre.inventory.list(
            "node_values_list", default=[1, 2, 3])
        nodeValuesList.meta['tip'] = "Position of desired values in UCD face nodal attributes."

        excludeZeroNormals = pythia.pyre.inventory.bool("exclude_zero_normals",
                                                        default=False)
        excludeZeroNormals.meta['tip'] = "Whether to exclude points with zero normals."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="grabfaces"):
        Application.__init__(self, name)
        return

    def main(self):
        # import pdb
        # pdb.set_trace()
        self._grabFaces()
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Setup members using inventory.
        """
        Application._configure(self)
        self.ucdFaceFile = self.inventory.ucdFaceFile
        self.faultIDNum = self.inventory.faultIDNum
        self.pointOutputFile = self.inventory.pointOutputFile
        self.nodeValuesList = self.inventory.nodeValuesList
        self.excludeZeroNormals = self.inventory.excludeZeroNormals
        return

    def _grabFaces(self):
        """Reads vertex coordinates, connectivities, and vertex attributes from a
        UCD file.
        """
        f = file(self.ucdFaceFile)
        lines = f.readlines()
        fileLen = len(lines)
        firstline = lines[0].split()
        numVerts = int(firstline[0])
        numCells = int(firstline[1])
        numVertAttrs = int(firstline[2])
        vertCoords = []
        # Get vertex coordinates
        for lineCount in range(1, numVerts+1):
            data = lines[lineCount].split()
            for dim in range(1, 4):
                vertCoords.append(float(data[dim]))

        # Get cell connectivities
        faultVerts = []
        lineBegin = 1 + numVerts
        lineEnd = lineBegin + numCells
        firstCellLine = lines[lineBegin].split()
        cellType = str(firstCellLine[2])
        if cellType == "tri":
            numVertsPerCell = 3
        else:
            numVertsPerCell = 4
        for lineCount in range(lineBegin, lineEnd):
            data = lines[lineCount].split()
            cellMat = int(data[1])
            if cellMat == self.faultIDNum:
                for vert in range(3, 2 + numVertsPerCell):
                    vertNum = int(data[vert])
                    if vertNum not in faultVerts:
                        faultVerts.append(vertNum)
        faultVerts.sort()
        numFaultVerts = len(faultVerts)

        # read normals/values and write out the selected values.
        o = open(self.pointOutputFile, 'w')
        lineBegin = 2 + numVerts + numCells + numVertAttrs
        lineEnd = lineBegin + numVerts
        vertInd = 0
        ucdInd = 1
        coordCount = 0
        normals = [0.0, 0.0, 0.0]
        v0 = int(self.nodeValuesList[0])
        v1 = int(self.nodeValuesList[1])
        v2 = int(self.nodeValuesList[2])
        for lineCount in range(lineBegin, lineEnd):
            vertex = faultVerts[vertInd]

            if vertex == ucdInd:
                data = lines[lineCount].split()
                normals = [float(data[v0]), float(data[v1]), float(data[v2])]
                outputPoint = not self.excludeZeroNormals
                outputPoint = outputPoint or \
                    normals[0] != 0.0 or \
                    normals[1] != 0.0 or \
                    normals[2] != 0.0
                if outputPoint:

                    for dim in range(3):
                        o.write(' %.12e' % vertCoords[coordCount + dim])

                    for dim in range(3):
                        o.write(' %.12e' % normals[dim])

                    o.write('\n')
                vertInd += 1
            coordCount += 3

            if vertInd == numFaultVerts:
                break
            ucdInd += 1

        f.close()
        o.close()
        return


# ----------------------------------------------------------------------
if __name__ == '__main__':
    app = GrabFaces()
    app.run()

# End of file
