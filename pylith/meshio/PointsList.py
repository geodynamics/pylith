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

from pythia.pyre.components.Component import Component


def validateFilename(value):
    """Validate filename with list of points.
    """
    if 0 == len(value):
        raise ValueError("Filename for list of points not specified.")
    return value


class PointsList(Component):
    """
    Reader for a list of points from an ASCII file.

    :::{seealso}
    See [`OutputSolnPoints` Component](OutputSolnPoints.md).
    :::
    """
    DOC_CONFIG = {
        "cfg": """
            [points]
            filename = stations.txt
            comment_delimiter = #
            value_delimiter = ,
            
            coordsys = spatialdata.geocoords.CSCart
            coordsys.space_dim = 2
        """
    }

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="", validator=validateFilename)
    filename.meta['tip'] = "Filename for list of points."

    commentDelimiter = pythia.pyre.inventory.str("comment_delimiter", default="#")
    commentDelimiter.meta['tip'] = "Delimiter for comments."

    valueDelimiter = pythia.pyre.inventory.str("value_delimiter", default=None)
    valueDelimiter.meta['tip'] = "Delimiter used to separate values."

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pythia.pyre.inventory.facility("coordsys", family="coordsys", factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with points."

    def __init__(self, name="pointslist"):
        """Constructor.
        """
        Component.__init__(self, name)

    def read(self):
        """Read points from file.
        """
        import numpy

        fin = open(self.filename, "r")
        lines = fin.readlines()
        fin.close()

        npoints = 0
        ndims = None
        for line in lines:
            if line.startswith(self.commentDelimiter):
                continue

            if not self.valueDelimiter:
                fields = line.split()
            else:
                fields = line.split(self.valueDelimiter)
            if ndims:
                if len(fields) != 1 + ndims:
                    raise IOError("Error occurred while reading line '%s' in points file '%s'.\n"
                                  "Expected format: station x y [z] with %d fields. Found %d fields." %
                                  (line, self.filename, 1 + ndims, len(fields)))
            else:
                ndims = len(fields) - 1
            npoints += 1

        points = numpy.zeros((npoints, ndims), dtype=numpy.float64)
        stations = []
        ipoint = 0
        for line in lines:
            if line.startswith(self.commentDelimiter):
                continue

            fields = line.split(self.valueDelimiter)

            stations.append(fields[0].strip())
            points[ipoint,:] = list(map(float, fields[1:]))

            ipoint += 1

        return stations, points

    def _configure(self):
        """Set members based using inventory.
        """
        Component._configure(self)


# FACTORIES ////////////////////////////////////////////////////////////

def points_list():
    """Factory associated with PointsList.
    """
    return PointsList()


# End of file
