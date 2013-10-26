# -*- coding: utf-8 -*-

from astropy import coordinates as c  # not recommended for users, but used below to make things simpler
from astropy import units as u

#---- Classes for the representation of coordinates ----
c.SphericalRepresentation(c.Latitude(...), c.Longitude(...) # order doesn't matter, determine from class
c.SphericalRepresentation(c.Latitude(...), c.Longitude(...), c.Distance(...))  #also can give distance
c.SphericalRepresentation(lat=5*u.deg, lon=8*u.hour)
c.SphericalRepresentation(lat=[5, 6]*u.deg, lon=[8, 9]*u.hour)  # arrays are of course fine
c.SphericalRepresentation(lat=[5, 6]*u.deg, lon=[8, 9]*u.hour, copy=False)  # default is to make copies
c.SphericalRepresentation(lat='5rad', lon='2h6m3.3s')  # these are parsed by `Latitude` and `Longitude` constructors, so no need to implement parsing

c1 = c.SphericalRepresentation(lat=5*u.deg, lon=8*u.hour, distance=10*u.kpc)  # these could also be Angle and Distance objects instead of Quantities
c1 = c.SphericalRepresentation(c1)  #makes a copy

c.SphericalRepresentation(lat=[5, 6]*u.deg, lon=[8, 9]*u.hour, distance=[10, 11]*u.kpc)  # distance, lat, and lon must match
with raises(ValueError):
    c.SphericalRepresentation(lat=[5, 6]*u.deg, lon=[8, 9]*u.hour, distance=[10, 11, 12]*u.kpc)  # distance, lat, and lon must match
with raises(ValueError):
    c.SphericalRepresentation(lat=[5, 6]*u.deg, lon=[8, 9, 10]*u.hour, distance=[10, 11]*u.kpc)  # distance, lat, and lon must match
c2 = c.SphericalRepresentation(lat=[5, 6]*u.deg, lon=[8, 9]*u.hour, distance=10*u.kpc)  # if distance is a scalar, assume that's meant for all the angles
assert len(c2.distance) == 2

c2 = c.SphericalRepresentation('5:10:20.52 +23:23:23.5', units=(u.hourangle, u.degree))
# In the current API, `unit` is confusing because sometimes it's one object and sometimes two (issue #1421).  
# But because this is the *only* time `units` is meaningful, it's ok here
assert c2.lon.units == u.hourangle

#OPTION: add store_as keyword, as discussed in #1421
c2 = c.SphericalRepresentation('5:10:20.52 +23:23:23.5', units=(u.hourangle, u.degree), store_as=(u.radian, u.radian))
assert c2.lon.units == u.radian
#end OPTION

with raises(ValueError):
    c.SphericalRepresentation('5:10:20.52 +23:23:23.5')  # this is ambiguous so it fails

#regardless of how input, they come out as angle/distance
assert isinstance(c1.lat, c.Angle)
assert isinstance(c1.distance, c.Distance)

#but they are read-only, as representations are immutible once created
with raises(AttributeError):
    c1.lat = c.Latitude(...)

#OPTION: also support "colatitude", internally stored only as `lat`
c2 = c.SphericalRepresentation(colat=85*u.deg, lon=8*u.hour)
assert c1.lat == c2.lat
assert c2.colat.degree == 90. - c2.lat.degree
#end OPTION

#OPTION: also support "phi" and "theta"
c3 = c.SphericalRepresentation(phi=120*u.deg, theta=85*u.deg)
assert c1.lat == c3.lat
assert c1.lin == c3.lon
assert c1.phi == c3.phi
#I think this is a bad idea, because phi and theta's definition depends on your field/undergraduate classroom.
#end OPTION

c1 = c.CartesianRepresentation(randn(3, 100) * u.kpc) #first dimension must be 3
assert c1.xyz.shape[0] == 0
assert c1.unit == u.kpc
assert c1.xyz.unit == 0  # not recommended, but available because `xyz` is a quantity
assert c1.x.shape[0] == 100
assert c1.y.shape[0] == 100
assert c1.z.shape[0] == 100
c.CartesianRepresentation(x=randn(100)*u.kpc, y=randn(100)*u.kpc, z=randn(100)*u.kpc)
with raises(UnitsError):
    #units must match
    c.CartesianRepresentation(x=randn(100)*u.kpc, y=randn(100)*u.kpc, z=randn(100)*u.pc)

#OPTION: allow raw array inputs and `units` keyword
c.CartesianRepresentation(x=randn(100), y=randn(100), z=randn(100), units=u.kpc)
#end OPTION

#OPTION/extension: add more representations like cylindrical or elliptical

#---------------"Low-level" classes"
c.ICRS