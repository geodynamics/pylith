# Coordinate systems

from spatialdata.geocoords.CSGeo import CSGeo
from spatialdata.geocoords.CSGeoProj import CSGeoProj

# Geographic lat/lon coordinates in WGS84 datum
def cs_geo():
    """Geographic lat/lon coordinates in WGS84 datum.
    """
    cs = CSGeo()
    cs.inventory.datumHoriz = "WGS84"
    cs.inventory.datumVert = "mean sea level"
    cs.inventory.spaceDim = 2
    cs._configure()
    cs.initialize()
    return cs


def cs_geo3D():
    """Geographic lat/lon/elev coordinates in WGS84 datum with mean sea
    level as vertical datum.

    """
    cs = CSGeo()
    cs.inventory.datumHoriz = "WGS84"
    cs.inventory.datumVert = "mean sea level"
    cs.inventory.spaceDim = 3
    cs._configure()
    cs.initialize()
    return cs


def cs_mesh():
    """Cartesian coordinates for mesh with Portland as the origin and mean
    sea level as vertical datum.

    """
    cs = CSGeoProj()
    cs.inventory.datumHoriz = "WGS84"
    cs.inventory.datumVert = "mean sea level"
    cs.inventory.spaceDim = 3
    cs.inventory.projector.inventory.projection = "tmerc"
    cs.inventory.projector.inventory.projOptions = "+lat_0=45.5231 +lon_0=-122.6765 +k=0.9996"
    cs.inventory.projector._configure()
    cs._configure()
    cs.initialize()
    return cs


def geoToMesh(xyz):
    """Convert coordinates from geographic coordinates to mesh coordinates.
    """
    from spatialdata.geocoords.Converter import convert
    convert(xyz, cs_mesh(), cs_geo3D())
    return
    

# End of file
