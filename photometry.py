# About
# =====
#
# Current photometry functionality is being developed in the photutils
# affiliated package, with the aim of ultimately merging into
# astropy.photometry, so the present API document presents the final
# user-level API once the code is merged into the core package. In
# addition, the present document gives no indications of how we deal
# with things under the hood - we are only concerned with the user-level
# API.
#
# Authors
# =======
#
# * Thomas Robitaille
#
# Requirements
# ============
#
# The main requirements of the API are:
#
# * Provide an interface for users to programmatically perform aperture and
#   PSF photometry.
#
# * Allow users to easily define custom apertures and PSFs for use with the
#   photometry routines.
#
# Usage
# =====
#
# Apertures
# ---------
#
# Apertures are objects that have a specific geometry, and are assumed
# to be centered on the origin. Values for sizes are assumed to be in pixels
# by default.

from astropy.photometry import CircularAperture
a1 = CircularAperture(1.)

from astropy.photometry import EllipticalAnnulusAperture
a2 = EllipticalAnnulusAperture(1.5, 2.0, 1.0, np.pi / 3.)

# The Aperture class can be sub-classed for user-defined shapes:

a3 = CustomAperture()

# The details of the attributes to set and methods to overload is not
# described here, since it is unimportant for the user.
#
# If WCS information is available in the data, then apertures can be specified
# in arcseconds, arcminutes, or degrees:

from astropy import units as u
a4 = CircularAperture(1. * u.arcsec)
a5 = CircularAperture(2.3 * u.arcmin)
a6 = CircularAperture(0.002 * u.degree)

# PSFs
# ----
#
# PSFs are 2-d maps of source shapes, and can be defined either
# analytically or numerically on a grid. Some common PSF shapes are
# provided for convenience.
#
# Define a gaussian PSF with sigma=2.

from astropy.photometry import GaussianPSF
p1 = GaussianPSF(sigma=2.)

# Define a gaussian with an arbitrary function:

from astropy.photometry import PSF

def my_custom_psf(x):
    return np.sin(x) / x

p2 = PSF(my_custom_psf)

# which, in this case, can also be written as a lambda function:

p2 = PSF(lambda x: np.sin(x) / x)

# PSFs can also be specified numerically:

p3 = PSF(psf_array, sampling=5)

# where the PSF array should be centered on the origin, and `sampling`
# indicates the oversampling factor of the PSF (defaults to 1).
#
# The PSF class can be sub-classed for specific instruments/missions,
# and for example to read PSFs from files (by default, the PSF class
# does not accept any files, since there is no single format for
# PSFs):

p4 = SpitzerPSF('irac_8.0_psf.fits')

# The details of the attributes to set and methods to overload is not
# described here, since it is unimportant for the user.
#
# Photometry
# ----------
#
# We need to ensure that the API for aperture and PSF photometry is as
# similar as possible to make it easier for users. The basic API for the
# aperture and PSF photometry is therefore of the form:
#
# results = <type>_photometry(<data>, <coordinates>, <aperture or PSF>)
#
# The simplest case is to do photometry on a 2-d Numpy array using pixel coordinates:

from astropy.io import fits
image = fits.getdata('image.fits')

x = np.array([1., 4., 2., 7.])
y = np.array([3., 7., 2., 3.])

from astropy.photometry import GaussianPSF
psf = GaussianPSF(sigma=2.)

from astropy.photometry import psf_photometry
results_psf = psf_photometry(image, (x, y), psf)

from astropy.photometry import CircularAperture
ap = CircularAperture(radius=3.)

from astropy.photometry import aperture_photometry
results_ap = aperture_photometry(image, (x, y), ap)

# The object returned is an astropy.table.Table object, which for all
# intents and purposes, behaves like a Numpy structured array, but also
# includes meta-data and allows easy I/O. The rationale for using a
# structured-like array is that any number of fields can be returned,
# and new fields can be added in future without breaking backward
# compatibility. The rationale for using an astropy.table.Table object
# is because it adds table-manipulation methods that users normally
# expect of structured arrays, but which structured arrays do not have.

# Make a plot of the flux vs. error

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.loglog(results['flux'], results['error'])

# Write out the photometry results to different files

results.write('photometry.txt', format='daophot')  # DAOphot table

results.write('photometry.xml', format='vo')  # VO table

# For convenience, any object initializing the PSF class can be
# passed directly to psf_photometry:

results = psf_photometry(data, galactic, lambda x: np.sin(x) / x)

# Integration with NDData
# -----------------------
#
# While the photometry can be performed on data in the form of a Numpy
# array, one can also pass an astropy.nddata.NDData object (or
# anything initializing an NDData object):

from astropy.nddata import NDData
data = NDData(image, mask=mask)
results = psf_photometry(data, (x, y), psf)

# which allows masks and uncertainties to be set on the NDData object
# (the behavior of the mask when doing photometry will need to be well
# documented) Not all types of uncertainties will be supported by
# psf_photometry and aperture_photometry, so if a type of uncertainty is
# passed which is not supported, an exception will be raised:

data = NDData(image, uncertainty=CustomUncertainty(unc_image))
results = psf_photometry(data, (x, y), psf)

# UnsupportedUncertaintyError: cannot perform PSF photometry on data
# with uncertainties of type CustomUncertainty.
#
# Once photutils is integrated into astropy.photometry, one could
# envisage also allowing these function to be accessible as a method on
# NDData, with exactly the same behavior:

results_psf = image.psf_photometry((x, y), psf)
results_ap = image.aperture_photometry((x, y), ap)

# Integration with astropy.wcs and astropy.coordinates
# ----------------------------------------------------

# If and NDData object with WCS is passed to psf_photometry, then
# non-pixel coordinates can be used:

galactic = GalacticCoordinates([1.,2.,3.], [4.,5.,6.])
results = psf_photometry(data, galactic, psf)

# The behavior of the rotation of the PSF will need to be well
# documented in cases where non-pixel coordinates are used, and the
# PSF is non-axisymmetric.
#
# Syntax
# ------
#
# Wherever possible, arguments will be named arguments so that code
# can more easily understandable, e.g.:

results = psf_photometry(data=image, coords=galactic,
                         psf=lambda x: np.sin(x) / x)

# Multiple apertures
# ------------------

# In some cases, one might want to do photometry for multiple
# apertures for a given source, or different apertures for each source
# (or both). For the simplest case where we want to do photometry for
# different apertures (but the same apertures for all sources), we can just
# envisage passing a list of apertures:

results = aperture_photometry(image, (x, y), [ap1, ap2, ap3])

# and similiarly to do PSF photometry with multiple PSFs. The
# resulting table will then have vector rather than scalar columns for
# columns that depend on the apertures (flux, error, etc.).
#
# For the case where we want a different aperture for each source, if each
# aperture is a different shape, then the user can just call
# aperture_photometry for each object, but if all the apertures are the same
# type (e.g. circular), then we can make things easier by initializing the
# aperture object with a 1-d sequence of properties instead of scalars:

aps1 = CircularAperture(radius=[1.,2.,3.])
results = aperture_photometry(image, (x, y), aps1)

# One can then combine this with the multiple apertures above, by
# specifying lists of array-based apertures:

aps1 = CircularAperture(radius=[1.,2.,3.])
aps2 = CircularAperture(radius=[1.5,2.5,3.5])
results = aperture_photometry(image, (x, y), [aps1, aps2])

# which as before, returns vector columns in the table.
