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
# * Thomas Robitaille (@astrofrog)
# * Kyle Barbary (@kbarbary)
# * Rene Breton (@bretonr)
# * David Shupe (@stargaser)
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

# and similarly for other aperture types.

# PSFs
# ----
#
# PSFs are 2-d maps of source shapes, and can be defined either
# analytically or numerically on a grid. Some common PSF shapes are
# provided for convenience.
#
# Define a gaussian PSF with sigma=2.
p1 = GaussianPSF(sigma=2.)

# Define a Moffat PSF
p2 = MoffatPSF(alpha=5., beta=0.1)

# Define a Lorentzian PSF
p3 = LorentzianPSF(gamma=2.)

# Define a gaussian with an arbitrary function:

from astropy.photometry import PSF

def my_custom_psf(x):
    return np.sin(x) / x

p4 = AnalyticalPSF(my_custom_psf)

# which, in this case, can also be written as a lambda function:

p5 = AnalyticalPSF(lambda x: np.sin(x) / x)

# PSFs can also be specified numerically:

p6 = DiscretePSF(psf_array, sampling=5)

# where the PSF array should be centered on the origin, and `sampling`
# indicates the oversampling factor of the PSF (defaults to 1).
#
# The PSF class can be sub-classed for specific instruments/missions,
# and for example to read PSFs from files (by default, the PSF class
# does not accept any files, since there is no single format for
# PSFs):

p7 = SpitzerPSF('irac_8.0_psf.fits')

# The details of the attributes to set and methods to overload is not
# described here, since it is unimportant for the user.
#
# For analytical PSFs, we need to worry about the truncation radius of the
# PSF. We parametrize this by an argument ``truncation`` which should be set
# to the radius at which the PSF can be truncated, e.g.:

p8 = GaussianPSF(sigma=3., truncation=20.)

# Building a PSF from a dataset should be made possible, using:

p_new = create_psf(image, (x, y), truncation=10., sampling=2, mode='median')

# the above will create a 40x40 DiscretePSF by median-combining the individual
# PSFs of stars at the locations specified by ``(x, y)``.
#
# We should allow for position-dependent PSFs and PSFs with distortions.
#
# We may also want to support PRFs (Point Response Functions) which are
# slightly different - see
#
# http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/calibrationfiles/psfprf/
#
# for more details. These are equivalent to over-sampled PSF, but are actually
# easier to construct and use in fitting because the oversampling/downsampling
# steps are not needed. Provided the PRF class has the same public interface
# as the PSF class, it can be used interchangeably.
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

# PSF-fitting algorithm(s)
# ------------------------

# By default, the PSF photometry function would attempt a simultaneous fit of
# the PSFs to the positions specified, but there are cases where a sequential
# approach might work well and be much faster, so we could have a ``mode``
# argument to specify how to operate:

# By default, fit PSFs simultaneously at all positions

results1 = psf_photometry(data, (x, y), psf)

# or

results1 = psf_photometry(data, (x, y), psf, mode='simultaneous')

# Fit sources sequentially
results2 = psf_photometry(data, (x, y), psf, mode='sequential')

# Centroiding/tuning coordinates
# ------------------------------
#
# We should provide a function to improve the accuracy of coordinates based on
# centroiding algorithms, inside a set of apertures. This could be done with:

coords_new = centroid(data, (x, y), <parameters>)

# where <parameters> are arguments that define how the centroiding is done.
#
# In the case of PSF photometry, the coordinates could be tuned during the
# fitting itself:

results = psf_photometry(data, (x, y), psf, tune_coordinates=True, tune_limit=7.)

# where ``tune_limit`` would be the maximum radius the coordinates can change
# by.

# Making residual images
# ----------------------

# The psf_photometry routine would take an argument that if specified, returns
# an array (or NDData object) containing the residual image.

# By default, only return table of results
results = psf_photometry(data, (x, y), psf)

# If requested, return residuals
results, residual = psf_photometry(data, (x, y), psf, residual=True)

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
# Photometry without NDData objects
# ---------------------------------
#
# Since it would be nice to be able to use the photometry functions without
# using the NDData class in all cases, we will allow ``wcs``, ``uncertainty``,
# and ``mask`` to be passed to the photometry routines, e.g.:

results = psf_photometry(image, (x, y), psf, wcs=WCS(...),
                         uncertainty=StdDevUncertainty(...))

# In addition, ``gain`` can also be passed, which is a scalar value to convert
# the data units into the number of photons. This is used to add an additional
# error term based on Poisson statistics.
#
# If an NDData object is passed as input, then these additional arguments
# cannot be specified.
#
# Syntax
# ------
#
# In order to ensure that code is more easily understandable, we should
# recommend that users always explicilty state the argument names:

results = psf_photometry(data=image, coords=galactic,
                         psf=lambda x: np.sin(x) / x)

# Multiple apertures/PSFs
# -----------------------

# In some cases, one might want to do photometry for multiple apertures or
# PSFs for a given source. In this case, one can simply pass a list of
# apertures or PSF objects:

results_ap = aperture_photometry(image, (x, y), [ap1, ap2, ap3])

results_psf = psf_photometry(image, (x, y), [psf1, psf2, psf3])

# We can also provide a single function that allows users to do either, or
# both aperture and PSF photometry:

results = photometry(image, (x, y), [ap1, psf1, ap2, psf2])

# Custom statistic for aperture photometry
# ----------------------------------------

# We may want to have `aperture_photometry` accept a keyword that specifies
# what statistic to calculate on pixels in the aperture:

aperture_photometry(image, (x, y), ap, statistic='sum')  # sum pixels in aperture (default)
aperture_photometry(image, (x, y), ap, statistic='mean')  # average of pixels in aperture
aperture_photometry(image, (x, y), ap, statistic='median')  # median of pixels in aperture

# In the above cases, the behavior for partial pixels should be well
# documented, especially in the case of 'median'. Users could define their own
# functions, which would need to have a 'weights' argument, so that the
# function can deal with partial pixel coverage. For example, one could
# implement a custom median function that uses only pixels with more than 50%
# overlap with the aperture:

def my_median(a, weights):
    """Median of pixels more than halfway in aperture"""
    return np.median(a[weights > 0.5])

aperture_photometry(image, (x, y), ap, statistic=my_median)