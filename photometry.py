# Regions
# -------

# PSF
# ---

# Define a gaussian PSF with sigma=2.

from photutils.psf import GaussianPSF

p1 = GaussianPSF(sigma=2.)

# Define a gaussian with an arbitrary function:

from photutils.psf import PSF

def my_custom_psf(x):
    return np.sin(x) / x

p2 = PSF(my_custom_psf)

# ... or as a lambda function:

p2 = PSF(lambda x: np.sin(x) / x)

# The PSF class can be sub-classed for specific instruments/missions,
# and for example to read PSFs from files (by default, the PSF class
# does not accept any files, since there is no single format for
# PSFs):

p3 = SpitzerPSF('irac_8.0_psf.fits')

# The details of the attributes to set and methods to overload is not
# described here, since it is unimportant for the user.

# PSF Photometry
# --------------

# To perform the actual photoemtry, the simplest is to pass a Numpy
# array with pixel coordinates at which to do the photometry, as well as
# the PSF itself:

from astropy.io import fits
image = fits.getdata('image.fits')

x = np.array([1., 4., 2., 7.])
y = np.array([3., 7., 2., 3.])

from photutils.psf import GaussianPSF
psf = GaussianPSF(sigma=2.)

from photutils.psf import psf_photometry
results = psf_photometry(image, (x, y), psf)

# The object returned is an astropy.table.Table object, which for all
# intents and purposes, behaves like a Numpy structured array, but also
# includes meta-data and allows easy I/O. The rationale for using a
# structured-like array is that any number of fields can be returned,
# and new fields can be added in future without breaking backward
# compatibility.

# Make a plot of the flux vs error

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

# While the photometry can be performed on data in the form of a Numpy
# array, one can also pass an astropy.nddata.NDData object (or
# anything initializing an NDData object) with uncertainties:

from astropy.nddata import NDData
data = NDData(image, mask=mask)
results = psf_photometry(data, (x, y), psf)

# which allows masks and uncertainties to be set on the NDData object
# (the behavior of the mask when doing photometry will need to be well
# documented) Not all types of uncertainties will be supported by
# psf_photometry, so if a type of uncertainty is passed which is not
# supported, an exception will be raised:

data = NDData(image, uncertainty=CustomUncertainty(unc_image))
results = psf_photometry(data, (x, y), psf)

# UnsupportedUncertaintyError: cannot perform PSF photometry on data
# with uncertainties of type CustomUncertainty.

# Once photutils is integrated into astropy.photometry, one could
# envisage also allowing this function to be accessible as a method on
# NDData, with exactly the same behavior:

results = image.psf_photometry((x, y), psf)

# Integration with astropy.wcs and astropy.coordinates
# ----------------------------------------------------

# If and NDData object with WCS is passed to psf_photometry, then
# non-pixel coordinates can be used:

galactic = GalacticCoordinates([1.,2.,3.], [4.,5.,6.])
results = psf_photometry(data, galactic, psf)

# The behavior of the rotation of the PSF will need to be well
# documented in cases where non-pixel coordinates are used, and the
# PSF is non-axisymmetric.

# Syntax
# ------

# The arguments will be named arguments so that code can more easily
# understandable:

results = psf_photometry(data=image, coords=galactic,
                         psf=lambda x: np.sin(x) / x)

# but this can only be done once photutils is part of the core astropy
# package.