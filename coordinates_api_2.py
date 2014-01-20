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



#<---------------------Reference Frame/"Low-level" classes--------------------->
#The low-level classes have a dual role: they act as specifiers of coordinate
#frames and they *may* also contain data as one of the representation objects,
#in which case they are the actual coordinate objects themselves.


#They can always accept a representation as a first argument
icrs = c.ICRS(c.SphericalRepresentation(lat=5*u.deg, lon=8*u.hour))


#Frames that require additional information like equinoxs or obstimes get them as
#keyword parameters to the frame constructor.  Where sensible, defaults are used
fk5 = c.FK5(c.SphericalRepresentation(lat=5*u.deg, lon=8*u.hour))  # FK5 is almost always J2000 equinox
J2000 = astropy.time.Time('J2000',scale='utc')
fk5_2000 = c.FK5(c.SphericalRepresentation(lat=5*u.deg, lon=8*u.hour), equinox=J2000)
assert fk5.equinox == fk5_2000.equionx

#the information required to specify the frame is immutible
fk5.equinox = J2001  # raises AttributeError

#OPTION: the *representation data* might not be immutible:
fk5.data = c.SphericalRepresentation(lat=6*u.deg, lon=9*u.hour)
#OR, it might be immutible:
fk5.data = ... #raises AttributeError
#end OPTION

#There is also a class-level attribute that lists the attributes needed to
#identify the frame.  These include attributes like the `equinox` above.
assert fk5.framespecattrs == ('equinox', 'obstime')
assert FK5.framespecattrs == ('equinox', 'obstime')  # defined on the *class*


#The actual position information is accessed via the representation objects
assert icrs.represent_as(c.SphericalRepresentation).lat == 5*u.deg
assert icrs.spherical.lat == 5*u.deg  # shorthand for the above
assert icrs.cartesian.z.value > 0

#Many frames have a "preferred" representation, the one in which they are
#conventionally described, often with a special name for some of the
#coordinates. E.g., most equatorial coordinate systems are spherical with RA and
#Dec. This works simply as a shorthand for the longer form above

assert icrs.ra == 5*u.deg
assert fk5.dec == 8*u.hour

assert icrs.preferred_representation == c.SphericalRepresentation

#the frames can also be initialized with the preferred names:
icrs_2 = c.ICRS(ra=8*u.hour, dec=5*u.deg, distance=1*u.kpc)
assert icrs == icrs2


#the frames also know how to give a reasonable-looking string of themselves,
#based on the preferred representation and possibly distance
assert str(icrs_2) == '<ICRS RA=120.000 deg, Dec=5.00000 deg, Distance=1 kpc>'


#<-------------------------Transformations------------------------------------->
#Transformation transform low-level classes from one frame to another.
#If no data (or `None`) is given, the class acts as a specifier of a frame, but
#without any stored data.
J2001 = astropy.time.Time('J2001',scale='utc')
fk5_J2001_frame = c.FK5(equinox=J2001)

#if they do not have data, the string instead is the frame specification
assert str(fk5_J2001_frame) == "<FK5 frame: equinox='J2000.000', obstime='B1950.000'>"

#These frames are primarily useful for specifying what a coordinate should be
# transformed *into*, as they are used by the `transform_to` method

newfk5 = fk5.transform_to(fk5_J2001_frame)  # precesses the point to the new equinox
assert newfk5.equinox == J2001

#classes can also be given to `transform_to`, which then uses the defaults for
#the frame information:
samefk5 = fk5.transform_to(c.FK5)
#`fk5` was initialized using default `obstime` and `equinox`, so:
assert samefk5.ra == fk5.ra and samefk5.dec == fk5.dec


#transforming to a new frame necessarily loses framespec information if it
#is not necessary for the new frame, so transforms are not necessarily
#round-trippable unless the frame is explicitly given:
fk5_2 =c.FK5(ra=8*u.hour, dec=5*u.deg, equinox=J2001)
ic_trans = fk5_2.transform_to(c.ICRS)
fk5_trans = fk5_2.transform_to(c.FK5)
assert fk5_2.ra == fk5_trans.ra  # AssertionError - fk5_trans is in the J2000
# equinox instead of J2001, so it does not have the same RA/Dec
fk5_trans_2 = fk5_2.transform_to(fk5_2001_frame)
assert fk5_2.ra == fk5_trans_2.ra  # Now all is fine because same equinox

#Trying to tansforming a frame with no data is of course an error:
c.FK5(equinox=J2001).transform_to(c.ICRS)  # ValueError

#To actually define a new transformation, the same basic scheme as in the
#0.2/0.3 coordinates framework is used - a graph of transform functions
#connecting various coordinate classes together.  The main changes are:
# 1) The transform functions now get the frame object they are trasnforming the
#    current data into.
# 2) Frames with additional information need to have a way to transform between
#    objects of the same class, but with different framespecinfo values

#An example transform function:
@c.dynamic_transform_matrix(SomeNewSystem, FK5)
def new_to_fk5(newobj, fk5frame):
    ot = newobj.obstime
    eq = fk5frame.equinox
    # ... build a *cartesian* transform matrix using `eq` that trasnforms from
    # the `newobj` frame as observed at `ot` to FK5 an equinox `eq`
    return matrix



#<---------------------------"High-level" class-------------------------------->
