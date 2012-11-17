# About
# =====
#
# A number of packages, such as pywcsgrid2 and APLpy, have developed ways of
# making plots in world coordinate systems (WCS) in Matplotlib, and since
# this is a common issue in Astronomy, we are interested in combining all the
# common functionality for inclusion in the Astropy core. The present document
# describes the proposed API.
#
# The aim is not to provide the simplest possible user-friendly API to the
# user, but to design an API consistent with the matplotlib API upon which
# user-friendly tools can be developed.
#
# Authors
# =======
#
# * Thomas Robitaille
# * Jae-Joon Lee
#
# Requirements
# ============
#
# The main requirements of the API are:
#
# * Provide a WCSAxes class/projection which, given a WCS projection object,
#   will automatically set up the infrastracture to show the relevant ticks,
#   labels, and grid.
#
# * Make it easy for users to plot and overlay data in world coordinates, and
#   in different coordinate systems.
#
# Usage
# =====
#
# Initialization
# --------------
#
# Object-oriented interface (similar for plt.figure() and true OO). This
# uses a WCSAxes class which takes the same argument as Axes, but with an
# additional arguement for the WCS object.

import matplotlib.pyplot as plt
from astropy.wcs import WCS, WCSAxes

fig = plt.figure()
ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=WCS('image.fits'))
fig.add_axes(ax)

# We can also provide pyplot-style initialization:

from astropy import wcs
wcs.subplot(1, 1, 1, wcs=wcs.WCS('image.fits'))
wcs.axes([0.1, 0.1, 0.8, 0.8], wcs=wcs.WCS('image.fits'))

# Images
# ------
#
# Plotting images as bitmaps or contours should be done via the usual
# matplotlib methods (which can accept e.g. Numpy arrays or PIL Image objects
# for RGB images)

ax.imshow(fits.getdata('image.fits'))
ax.imshow(Image.open('rgb_image.png'))
ax.contour(fits.getdata('scuba.fits'))
ax.contourf(fits.getdata('scuba.fits'))

# Since WCS information technically doesn't require any information about
# the image size, there would be no such thing as a 'mismatch' of dimensions
# compared to the WCS. In the case where we want to overplot an image with a
# different WCS projection, we could do e.g.:

ax[other_wcs].contour(scuba_image)

# Multi-dimensional WCS objects
# -----------------------------

# The WCS object can have more than two dimensions, but since we are
# ultimately plotting two-dimensional data, we have to select which dimensions
# are used for the x and y pixel coordinates in the axes, and we also have to
# select which slices are selected for the dimensions that are not being
# shown. WCSAxes will also need information about the pixel coordinates in the
# other dimensions (i.e. the 'slice' position). We can try and encapsulate
# this in a single 'slice' argument which should be specified as a list, and
# where each element is either an integer (indicating the slice position) or a
# string containing 'x' or 'y', indicating the dimension that should be used
# for the x and y pixel coordinates.
#
# The user could set the slices to use when instantiating WCSAxes:

ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=WCS('cube.fits'),
             slice=['x', 3, 'y'])

# or after instantiation:

ax.set_slice(['x', 3, 'y'])

# or when calling imshow/contour/other methods that accept 2-d data:

ax.imshow(image, slice=['x', 3, 'y'])

# In the above examples, the first and third dimension will be used as the x
# and y pixel coordinates in the Axes, and the second pixel coordinate will be
# set to 4 (could also be set to a floating-point value).
#
# Info: we tried using separate dimensions= and slices= argument in APLpy (to
# separate the selection of dimensions from the slices) but since the two are
# closely related, an approach such as the above is simpler.
#
# Allowing the slice to be changed on the fly will make it easy to develop
# for example a utility to slice through cubes without having to
# re-instantiate the WCSAxes object for each slice.
#
# Coordinate system
# -----------------
#
# To set the coordinate system to use for the labels/grid:

ax.set_coordinate_system('fk5')
ax.set_coordinate_system('gal')

# We should have the following built-in coordinate systems:
#
# * 'fk4' or 'b1950': B1950 equatorial coordinates
# * 'fk5' or 'j2000': J2000 equatorial coordinates
# * 'gal' or 'galactic': Galactic coordinates
# * 'ecl' or 'ecliptic': Ecliptic coordinates
# * 'sgal' or 'supergalactic': Super-Galactic coordinates
#
# In addition, we should make it easy to make plots offset from given
# coordinates, which we could do with a 'offset' argument which takes a tuple
# of coordinates:

ax.set_coordinate_system('gal', offset=(30., 1.2))
ax.set_coordinate_system('fk5', offset=(272.3, 44.3))

# Finally, we should make it easy for users to specify their own coordinate
# system.

# Coordinate grid
# ---------------

# To overlay a grid, the normal matplotlib method can be used:

ax.grid()

# and to modify the visual properties of the grid, users can capture a
# handle to the grid object:

g = ax.grid()
g.set_alpha(0.5)
g.set_color('red')

# Axis labels
# -----------
#
# The format of the labels can be specified with:

ax.set_xticklabels_format('x.xxx')  # decimal, non-angle coordinates, 3 decimal places

ax.set_xticklabels_format('d.ddddd')  # decimal degrees, 5 decimal places

ax.set_xticklabels_format('dd:mm')  # sexagesimal, 1' precision
ax.set_xticklabels_format('dd:mm:ss')  # sexagesimal, 1" precision
ax.set_xticklabels_format('dd:mm:ss.ss')  # sexagesimal, 0.01" precision

ax.set_xticklabels_format('hh:mm')  # sexagesimal (hours)
ax.set_xticklabels_format('hh:mm:ss')  # sexagesimal (hours)
ax.set_xticklabels_format('hh:mm:ss.ss')  # sexagesimal (hours),

# and similarly for set_yticklabels_format.
#
# Question: how do we set the format for e.g. an offset in arcseconds?
# 's.ssss' is ambiguous because it could be arcseconds, or seconds of an hour.
#
# Question: at the moment, there is nothing in this document about what axes
# users want to show which ticks and coordinates on. For example, they might
# want to make a plot which shows the RA on both the top and right axes, and
# the Dec on the left and bottom axes, or combinations like that. I'm not sure
# what the best way is to make this possible/easy. An example of that is the
# top left plot in:
#
#       http://www.atnf.csiro.au/people/mcalabre/WCS/PGSBOX/index.html
#
# Tick/label spacing
# ------------------
#
# The spacing of ticks/tick labels should default to something sensible, but
# users often want to be able to specify the spacing of the ticks. The tick
# positions could be set manually:

ax.set_xticks([242.2, 242.3, 242.4])

# but it would also be nice to allow the spacing to be set:

ax.set_xticks_spacing(0.1)

# Question: allowing users to specify their own tick position or spacing,
# and the label format raises the issue of what happens when the value cannot
# be well represented by the format (this confused a lot of users in APLpy
# initially). For example, if one specifies the format to be 'dd:mm:ss' and a
# tick spacing of 0.001 degrees then the ticks will go:
#
# 13:42:00
# 13:42:01
# 13:42:03
# 13:42:04
# 13:42:06
#
# and so on, which confuses users. In APLpy, we now raise an exception if
# the tick spacing is not compatible with the labelling, but is this the best
# way to go? This also means that if a tick spacing is not specified, the
# minimum tick spacing should be set by the minimum one representable by the
# format. For example, if the format is set to 'dd:mm', then the default tick
# spacing cannot end up being smaller than 1' (and if it does it gets reset to
# 1').
#


# Alternative API for Labels

# For sky coordinate systems with high curvature (for example, near the
# pole), the coupling between x- and y-axis to actual coordinates
# become less well-defined.
#
# For example :  http://www.atnf.csiro.au/people/mcalabre/WCS/PGSBOX/pgsbox4.gif
#
# Therefore, the meaning of the method like "set_xticks" become ambiguous.
#
# One option is to use a number 1 & 2 to specify the sky coordinates.

ax.set_ticklabel1_format("hms") # 1st coordinate (e.g., R.A.)
ax.set_ticklabel2_format("dms") # 2nd coordinate (e.d., Dec)

#
ax.set_xaxis_coordinate(2) # x-axis show tick and ticklabels for 2nd
                           # sky coordinate.
ax.set_yaxis_coordinate(1) # x-axis show tick and ticklabels for 1st
                           # sky coordinate.

ax.set_yaxis_coordinate(1, which="right")
# only the right y-axis show tick and ticklabels for the 1st sky coordinate.

# to adjust tick locations
ax.set_tick1_params(spacing=2, unit="m") # spacing in 2 minutes

ax.set_tick1_params(nbins=4) # (approximately) 4 ticks

ax.set_tick1_params(locs=[242.2, 242.3, 242.4])

ax.set_tick1_params(locs=[242.2, 242.3, 242.4],
                    labels=["242.2", "242.3", "242.4"])
# ticks in specified locations with specified labels

# These keyword parameters will be allowed for set_ticklabel[12]_format.
ax.set_ticklabel1_format("hms", nbins=4)




# Patches/shapes/lines
# --------------------
#
# To overlay arbitrary matplotlib patches in pixel coordinates:

from matplotlib.patches import Rectangle
r = Rectangle((43., 44.), 23., 11.)
ax['pixel'].add_patch(r)

# To overlay arbitrary matplotlib patches in world coordinates:

from matplotlib.patches import Rectangle
r = Rectangle((272., -44.), 1., 0.8)
ax['fk5'].add_patch(r)

# And similarly:

ax['gal'].add_collection(c)
ax['fk4'].add_line(l)

# The same applies to normal plotting methods such as scatter/plot/etc:

ax['gal'].scatter(l, b)

# Question: If the coordinate system is not specified, does it default to
# pixel coordinates or the default world coordinate system of the WCS?
#
# Multiple coordinate systems
# ---------------------------
#
# Question: users might want to show multiple coordinate systems at the same
# time (e.g. Galactic and Equatorial). They could simply create two WCSAxes
# instances and make one empty and transparent, but the axes will not be
# linked. Should:

ax['gal'].grid()

# work and plot a Galactic grid? If so, then should:

ax['gal'].set_xtick_spacing(0.1)

# and so on work? Or should we have a since method that spawns a new linked WCSAxes object, e.g.:

ax_gal = ax.get_axes('gal')

# or something similar?
#
# Other issues
# ------------
#
# * It would be useful to provide conveniences for overlaying beams and
#   compasses for example, but this does not need to be included in the
#   original API.
