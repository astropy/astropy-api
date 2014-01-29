# -*- coding: utf-8 -*-

from astropy import coordinates as coords
from astropy import units as u

#<-----------------Classes for representation of coordinate data--------------->
#These classes inherit from a common base class and internally contain Quantity
#objects, which are arrays (although they may act as scalars, like numpy's
#length-0  "arrays")

#If `Latitude` or `Logitude` objects are given, order doesn't matter
coords.SphericalRepresentation(coords.Latitude(...), coords.Longitude(...)

#similar if `Distance` object is given
coords.SphericalRepresentation(coords.Latitude(...), coords.Longitude(...), coords.Distance(...))
coords.SphericalRepresentation(lat=5*u.deg, lon=8*u.hour)

#Arrays of any of the inputs are fine
coords.SphericalRepresentation(lat=[5, 6]*u.deg, lon=[8, 9]*u.hour)

#Default is to copy arrays, but optionally, it can be a reference
coords.SphericalRepresentation(lat=[5, 6]*u.deg, lon=[8, 9]*u.hour, copy=False)

#strings are parsed by `Latitude` and `Longitude` constructors, so no need to
#implement parsing in the Representation classes
coords.SphericalRepresentation(lat='5rad', lon='2h6m3.3s')

#Or, you can give `Quantity`s with keywords, and they will be internally
#converted to Angle/Distance
c1 = coords.SphericalRepresentation(lat=5*u.deg, lon=8*u.hour, distance=10*u.kpc)

# Can always just accept another representation, and that means just copy it
c1 = coords.SphericalRepresentation(c1)

 # distance, lat, and lon must match
coords.SphericalRepresentation(lat=[5, 6]*u.deg, lon=[8, 9]*u.hour, distance=[10, 11]*u.kpc)
with raises(ValueError):
    coords.SphericalRepresentation(lat=[5, 6]*u.deg, lon=[8, 9]*u.hour, distance=[10, 11, 12]*u.kpc)
with raises(ValueError):
    coords.SphericalRepresentation(lat=[5, 6]*u.deg, lon=[8, 9, 10]*u.hour, distance=[10, 11]*u.kpc)
#if the inputs are not the same, if possible they will be broadcast following
#numpy's standard broadcasting rules.
c2 = coords.SphericalRepresentation(lat=[5, 6]*u.deg, lon=[8, 9]*u.hour, distance=10*u.kpc)
assert len(c2.distance) == 2
#when they can't be broadcast, it is a ValueError (same as Numpy)
c2 = coords.SphericalRepresentation(lat=[5, 6]*u.deg, lon=[8, 9, 10]*u.hour) #raises ValueError

#OPTION: have the representation objects handle strings and lists of strings.
#if this option is *not* accepted, instead it will be in the "high-level"
#class discussed at the bottom of this document
c2 = coords.SphericalRepresentation('5:10:20.52 +23:23:23.5', units=(u.hourangle, u.degree))
# In the current API, `unit` is confusing because sometimes it's one object
# and sometimes two (issue #1421). But because this is the *only* time `units`
# is meaningful, it's ok here
assert c2.lon.units == u.hourangle

with raises(ValueError):
    coords.SphericalRepresentation('5:10:20.52 +23:23:23.5')  # this is ambiguous so it fails

#a `store_as` option is important for some cases - see astropy/astropy#1421 for justification
c2 = coords.SphericalRepresentation('5:10:20.52 +23:23:23.5', units=(u.hourangle, u.degree), store_as=(u.radian, u.radian))
assert c2.lon.units == u.radian

#end OPTION


#regardless of how input, the `lat` and `lon` come out as angle/distance
assert isinstance(c1.lat, coords.Angle)
assert isinstance(c1.distance, coords.Distance)

#but they are read-only, as representations are immutible once created
with raises(AttributeError):
    c1.lat = coords.Latitude(...)
#Note that it is still possible to modify the array in-place, but this is not
#sanctioned by the API, as this would prevent things like caching.
c1.lat[:] = [0] * u.deg  # possible, but NOT SUPPORTED

#OPTION: also support "colatitude", internally stored only as `lat`
c2 = coords.SphericalRepresentation(colat=85*u.deg, lon=8*u.hour)
assert c1.lat == c2.lat
assert c2.colat.degree == 90. - c2.lat.degree
#end OPTION

#OPTION: also support "phi" and "theta"
c3 = coords.SphericalRepresentation(phi=120*u.deg, theta=85*u.deg)
assert c1.lat == c3.lat
assert c1.lin == c3.lon
assert c1.phi == c3.phi
#I think this is a bad idea, because phi and theta's definition depends on your
#subfield/undergraduate classroom.
#end OPTION
#OPTION (suggested by @taldcroft): instead of the above, include subclasses that
#have the appropriate names
c3 = coords.PhysicistSphericalRepresentation(phi=120*u.deg, theta=85*u.deg, rho=3*u.kpc)
#end OPTION

#first dimension must be length-3 if a lone `Quantity` is passed in.
c1 = coords.CartesianRepresentation(randn(3, 100) * u.kpc)
assert c1.xyz.shape[0] == 0
assert c1.unit == u.kpc
#using `c1.xyz.unit` is the same as `c1.unit`, so it's not necessary to do this,
#but it illustrates that `xyz` is a full-fledged `Quantity`.
assert c1.xyz.unit == u.kpc
assert c1.x.shape[0] == 100
assert c1.y.shape[0] == 100
assert c1.z.shape[0] == 100
# can also give each as separate keywords
coords.CartesianRepresentation(x=randn(100)*u.kpc, y=randn(100)*u.kpc, z=randn(100)*u.kpc)
#if the units don't match but are all distances, they will automatically be
#converted to match `x`
xarr, yarr, zarr = randn(3, 100)
c1 = coords.CartesianRepresentation(x=xarr*u.kpc, y=yarr*u.kpc, z=zarr*u.kpc)
c2 = coords.CartesianRepresentation(x=xarr*u.kpc, y=yarr*u.kpc, z=zarr*u.pc)
assert c2.unit == u.kpc
assert c1.z.kpc / 1000 == c2.z.kpc
# although if they are not `Distance` compatible, it's an error
c2 = coords.CartesianRepresentation(x=randn(100)*u.kpc, y=randn(100)*u.kpc, z=randn(100)*u.deg) # raises UnitsError

#OPTION: allow raw array inputs and `units` keyword
coords.CartesianRepresentation(x=randn(100), y=randn(100), z=randn(100), units=u.kpc)
#end OPTION

#representations convert into other representations via the `represent_as` method
srep = coords.SphericalRepresentation(lat=0*u.deg, lon=90*u.deg, distance=1*u.pc)
crep = srep.represent_as(coords.CartesianRepresentation)
assert crep.x.value == 0 and crep.y.value == 1 and crep.z.value == 0
#The functions that actually do the conversion are defined via methods on the
#representation classes. This may later be expanded into a full registerable
#transform graph like the coordinate frames, but initially it will be a simpler
#method system


#Future extensions: add additional useful representations like cylindrical or elliptical


#<---------------------Reference Frame/"Low-level" classes--------------------->
#The low-level classes have a dual role: they act as specifiers of coordinate
#frames and they *may* also contain data as one of the representation objects,
#in which case they are the actual coordinate objects themselves.


#They can always accept a representation as a first argument
icrs = coords.ICRS(coords.SphericalRepresentation(lat=5*u.deg, lon=8*u.hour))


#Frames that require additional information like equinoxs or obstimes get them
#as keyword parameters to the frame constructor.  Where sensible, defaults are
#used. E.g., FK5 is almost always J2000 equinox
fk5 = coords.FK5(coords.SphericalRepresentation(lat=5*u.deg, lon=8*u.hour))
J2000 = astropy.time.Time('J2000',scale='utc')
fk5_2000 = coords.FK5(coords.SphericalRepresentation(lat=5*u.deg, lon=8*u.hour), equinox=J2000)
assert fk5.equinox == fk5_2000.equionx

#the information required to specify the frame is immutible
fk5.equinox = J2001  # raises AttributeError

#OPTION: the *representation data* might not be immutible:
fk5.data = coords.SphericalRepresentation(lat=6*u.deg, lon=9*u.hour)
#OR, it might be immutible:
fk5.data = ... #raises AttributeError
#end OPTION

#There is also a class-level attribute that lists the attributes needed to
#identify the frame.  These include attributes like the `equinox` above.
assert FK5.framespecattrs == ('equinox', 'obstime')  # defined on the *class*
assert fk5.framespecattrs == ('equinox', 'obstime')  # and hence also in the objects


#The actual position information is accessed via the representation objects
assert icrs.represent_as(coords.SphericalRepresentation).lat == 5*u.deg
assert icrs.spherical.lat == 5*u.deg  # shorthand for the above
assert icrs.cartesian.z.value > 0

#Many frames have a "preferred" representation, the one in which they are
#conventionally described, often with a special name for some of the
#coordinates. E.g., most equatorial coordinate systems are spherical with RA and
#Decoords. This works simply as a shorthand for the longer form above

assert icrs.ra == 5*u.deg
assert fk5.dec == 8*u.hour

assert icrs.preferred_representation == coords.SphericalRepresentation

#low-level classes can also be initialized with the preferred names:
icrs_2 = coords.ICRS(ra=8*u.hour, dec=5*u.deg, distance=1*u.kpc)
assert icrs == icrs2

#and these are taken as the default if keywords are not given:
icrs_nokwarg = coords.ICRS(8*u.hour, 5*u.deg, distance=1*u.kpc)
assert icrs_nokwarg.ra == icrs_2.ra and icrs_nokwarg.dec == icrs_2.dec

#they also are capable of computing on-sky or 3d separations from each other,
#which will be a direct port of the existing methods:
coo1 = coords.ICRS(ra=0*u.hour, dec=0*u.deg)
coo2 = coords.ICRS(ra=0*u.hour, dec=1*u.deg)
# `separation` is the on-sky separation
assert coo1.separation(coo2).degree == 1.0

#while `separation_3d` includes the 3D distance information
coo3 = coords.ICRS(ra=0*u.hour, dec=0*u.deg, distance=1*u.kpc)
coo4  = coords.ICRS(ra=0*u.hour, dec=0*u.deg, distance=2*u.kpc)
assert coo3.separation_3d(coo4).kpc == 1.0

#The next example raises an error because `coo1` and `coo2` don't have distances
assert coo1.separation_3d(coo2).kpc == 1.0  # raises ValueError


#the frames also know how to give a reasonable-looking string of themselves,
#based on the preferred representation and possibly distance
assert str(icrs_2) == '<ICRS RA=120.000 deg, Dec=5.00000 deg, Distance=1 kpc>'


#<-------------------------Transformations------------------------------------->
#Transformation transform low-level classes from one frame to another.
#If no data (or `None`) is given, the class acts as a specifier of a frame, but
#without any stored data.
J2001 = astropy.time.Time('J2001',scale='utc')
fk5_J2001_frame = coords.FK5(equinox=J2001)

#if they do not have data, the string instead is the frame specification
assert str(fk5_J2001_frame) == "<FK5 frame: equinox='J2000.000', obstime='B1950.000'>"

#These frames are primarily useful for specifying what a coordinate should be
#transformed *into*, as they are used by the `transform_to` method
# E.g., this snippet precesses the point to the new equinox
newfk5 = fk5.transform_to(fk5_J2001_frame)
assert newfk5.equinox == J2001

#classes can also be given to `transform_to`, which then uses the defaults for
#the frame information:
samefk5 = fk5.transform_to(coords.FK5)
#`fk5` was initialized using default `obstime` and `equinox`, so:
assert samefk5.ra == fk5.ra and samefk5.dec == fk5.dec


#transforming to a new frame necessarily loses framespec information if that
#information is not applicable to the new frame.  This means transforms are not
#always round-trippable:
fk5_2 =coords.FK5(ra=8*u.hour, dec=5*u.deg, equinox=J2001)
ic_trans = fk5_2.transform_to(coords.ICRS)

#`ic_trans` does not have an `equinox`, so now when we transform back to FK5,
#it's a *different* RA and Dec
fk5_trans = ic_trans.transform_to(coords.FK5)
assert fk5_2.ra == fk5_trans.ra  # raises AssertionError

# But if you explicitly give the right equinox, all is fine
fk5_trans_2 = fk5_2.transform_to(coords.FK5(equinox=J2001))
assert fk5_2.ra == fk5_trans_2.ra

#Trying to tansforming a frame with no data is of course an error:
coords.FK5(equinox=J2001).transform_to(coords.ICRS)  # raises ValueError


#To actually define a new transformation, the same scheme as in the
#0.2/0.3 coordinates framework can be re-used - a graph of transform functions
#connecting various coordinate classes together.  The main changes are:
# 1) The transform functions now get the frame object they are trasnforming the
#    current data into.
# 2) Frames with additional information need to have a way to transform between
#    objects of the same class, but with different framespecinfo values

#An example transform function:
@coords.dynamic_transform_matrix(SomeNewSystem, FK5)
def new_to_fk5(newobj, fk5frame):
    ot = newobj.obstime
    eq = fk5frame.equinox
    # ... build a *cartesian* transform matrix using `eq` that trasnforms from
    # the `newobj` frame as observed at `ot` to FK5 an equinox `eq`
    return matrix

#other options for transform functions include one that simply returns the new
#coordinate object, and one that returns a cartesian matrix but does *not*
#require `newobj` or `fk5frame` - this allows more optimization of the transform


#<---------------------------"High-level" class-------------------------------->
#The "high-level" class is intended to wrap the lower-level classes in such a
#way that they can be round-tripped, as well as providing a variety of convenience
#functionality.  This document is not intended to show *all* of the possible
#high-level functionality, rather how the high-level classes are initialized and
#interact with the low-level classes

#this creates an object that contains an `ICRS` low-level class, initialized
#identically to the first ICRS example further up.
sc = coords.SkyCoordinate(coords.SphericalRepresentation(lat=5*u.deg, lon=8*u.hour, distance=1*u.kpc), system='icrs')
#Other representations and `system` keywords delegate to the appropriate
#low-level class. The already-existing registry for user-defined coordinates
#will be used by `SkyCoordinate` to figure out what various the `system`
#keyword actually means.

#they can also be initialized using the preferred representation names
sc = coords.SkyCoordinate(ra=8*u.hour, dec=5*u.deg, system='icrs')
sc = coords.SkyCoordinate(l=120*u.deg, b=5*u.deg, system='galactic')

#High-level classes can also be initialized directly from low-level objects
sc = coords.SkyCoordinate(coords.ICRS(ra=8*u.hour, dec=5*u.deg))
#The next example raises an error because the high-level class must always
#have position data
sc = coords.SkyCoordinate(coords.FK5(equinox=J2001))  # raises ValueError

#similarly, the low-level object can always be accessed
assert str(scoords.frame) == '<ICRS RA=120.000 deg, Dec=5.00000 deg>'

#Should (eventually) support a variety of possible complex string formats
#OPTION: if this is implemented on the low-level class, this would also delegate
sc = coords.SkyCoordinate('8h00m00s +5d00m00.0s', system='icrs')
#In the next example, the unit is only needed b/c units are ambiguous.  In
#general, we *never* accept ambiguity
sc = coords.SkyCoordinate('8:00:00 +5:00:00.0', unit=(u.hour, u.deg), system='icrs')
#The next one would yield length-2 array coordinates, because of the comma
sc = coords.SkyCoordinate(['8h 5d', '2°5\'12.3" 0.3rad'], system='icrs')
#It should also interpret common designation styles as a coordinate
sc = coords.SkyCoordinate('SDSS J123456.89-012345.6', system='icrs')

#the string representation can be inherited from the low-level class.
assert str(sc) == '<SkyCoordinate (ICRS) RA=120.000 deg, Dec=5.00000 deg>'
#but it should also be possible to provide formats for outputting to strings,
#similar to `Time`.  This can be added right away or at a later date.

#transformation is done the same as for low-level classes, which it delegates to
scfk5_j2001 = scoords.transform_to(coords.FK5(equinox=J2001))

#the key difference is that the high-level class remembers frame information
#necessary for round-tripping, unlike the low-level classes:
sc1 = coords.SkyCoordinate(ra=8*u.hour, dec=5*u.deg, equinox=J2001, system='fk5')
sc2 = sc1.transform_to(coords.ICRS)
#The next assertion succeeds, but it doesn't mean anything for ICRS, as ICRS
#isn't defined in terms of an equinox
assert sc2.equinox == J2001
#But it *is* necessary once we transform to FK5
sc3 = sc2.transform_to(coords.FK5)
assert sc3.equinox == J2001
assert sc1.ra == sc3.ra
#Note that this did *not* work in the low-level class example shown above,
#because the ICRS low-level class does not store `equinox`.


#`SkyCoordinate` will also include the attribute-style access that is in the
#v0.2/0.3 coordinate objects.  This will *not* be in the low-level classes
sc = coords.SkyCoordinate(ra=8*u.hour, dec=5*u.deg, system='icrs')
scgal = scoords.galactic
assert str(scgal) == '<SkyCoordinate (Galactic) l=216.31707 deg, b=17.51990 deg>'


#the existing `from_name` and `match_to_catalog_*` methods will be moved to the
#high-level class as convinience functionality.

m31icrs = SkyCoordinate.from_name('M31', system='icrs')
assert str(m31icrs) == '<SkyCoordinate (ICRS) RA=10.68471 deg, Dec=41.26875 deg>'

cat1 = SkyCoordinate(ra=[...]*u.hr, dec=[...]*u.deg, distance=[...]*u,kpc)
cat2 = SkyCoordinate(ra=[...]*u.hr, dec=[...]*u.deg, distance=[...]*u,kpc
idx2, sep2d, dist3d = cat1.match_to_catalog_sky(cat2)
idx2, sep2d, dist3d = cat1.match_to_catalog_3d(cat2)


#additional convinience functionality for the future should be added as methods
#on `SkyCoordinate`, *not* the low-level classes.
