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

#Future extensions: add more representations like cylindrical or elliptical



#<---------------------"Low-level"/Reference Frame classes--------------------->
#The low-level classes have a dual role: they act as specifiers of coordinate
#frames and they *may* also contain representations, in which case they are the
#actual coordinate objects themselves.


#They can always accept a representation as a first argument
icrs = c.ICRS(c.SphericalRepresentation(lat=5*u.deg, lon=8*u.hour))


#Frames that require additional data like equinoxs or obstimes get them as
#keyword parameters to the frame constructor.  Where sensible, defaults are used
fk5 = c.FK5(c.SphericalRepresentation(lat=5*u.deg, lon=8*u.hour))  # FK5 is almost always J2000 equinox
J2000 = astropy.time.Time('J2000',scale='utc')
fk5_2000 = c.FK5(c.SphericalRepresentation(lat=5*u.deg, lon=8*u.hour), equinox=J2000)
assert fk5.equinox == fk5_2000.equionx

#the frame data are immutible
J2001 = astropy.time.Time('B1950',scale='utc')
fk5.equinox = J2001  # raises AttributeError


#If no data or None is given, the class acts as a specifier of a frame, but
#without any stored data.  These frames are primarily useful for specifying
#what a coordinate should be trasnformed *into*:

newfk5 = fk5.transform_to(c.FK5(equinox=J2001))  # precesses the data point to the new equinox
assert newfk5.equinox == J2001

#If there is data, it can be accessed through representation objects
assert icrs.represent_as(c.SphericalRepresentation).lat == 5*u.deg
assert icrs.spherical.lat == 5*u.deg  # shorthand for the above
assert icrs.cartesian.z.value > 0

#Most frames have a "preferred" representation, usually the one in which they are
#typically shown, often with a special name for some of the coordinates. They can
#be accessed this way as a equivalent shorthand for the representation form
assert icrs.ra == 5*u.deg
assert fk5.dec == 8*u.hour

#the frames can also be initialized with the "preferred" names:
icrs_2 = c.ICRS(ra=8*u.hour, dec=5*u.deg, distance=1*u.kpc)
assert icrs == icrs2
fk5_2 = c.FK5(ra=8*u.hour, dec=5*u.deg, distance=1*u.kpc, equinox=J2000)


#the frames also know how to give a reasonable-looking string of themselves,
#based on the preferred coordinate system
assert str(icrs_2) == '<ICRS RA=120.000 deg, Dec=5.00000 deg, Distance=1 kpc>'
