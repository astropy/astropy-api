About
=====

A number of packages, such as pywcsgrid2, APLpy, and kapteyn, have developed
ways of making plots in world coordinate systems (WCS) in Matplotlib, and
since this is a common issue in Astronomy, we are interested in combining all
the common functionality for inclusion in the Astropy core. The present
document describes the proposed API.

The aim is not to provide the simplest possible user-friendly API to the
user, but to design an API consistent with the matplotlib API upon which
user-friendly tools can be developed.

Authors
=======

* Thomas Robitaille (@astrofrog)
* Jae-Joon Lee (@leejjoon)

Requirements
============

The main requirements of the API are:

* Provide a ``WCSAxes`` class/projection which, given a WCS projection object,
  will automatically set up the infrastructure to show the relevant ticks,
  labels, and grid.

* Make it easy for users to plot and overlay data in world coordinates, and
  in different coordinate systems.

* Support n-dimensional datasets

Usage
=====

Initialization
--------------

The object-oriented interface (similar for ``plt.figure()`` and the true OO
API) uses a ``WCSAxes`` class which takes the same argument as Axes, but with an
additional argument for the WCS object.

    import matplotlib.pyplot as plt
    from astropy.wcs import WCS, WCSAxes

    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=WCS('image.fits'))
    fig.add_axes(ax)

We can also provide pyplot-style initialization:

    from astropy import wcs
    wcs.subplot(1, 1, 1, wcs=wcs.WCS('image.fits'))
    wcs.axes([0.1, 0.1, 0.8, 0.8], wcs=wcs.WCS('image.fits'))

Coordinate systems
------------------

As a convention, and to avoid any assumptions, if methods are called as:

    ax.method()

then the method applies to pixel coordinates. To plot in world coordinates,
one should use:

    ax['world'].method()

One can also explicitly specify pixel coordinates:

    ax['pixel'].method()

If the two world coordinates corresponding to the pixel axes define sky
coordinates, we can also use:

    ax[coordinate_system].method()

We should have the following built-in sky coordinate systems:

* 'fk4' or 'b1950': B1950 equatorial coordinates
* 'fk5' or 'j2000': J2000 equatorial coordinates
* 'gal' or 'galactic': Galactic coordinates
* 'ecl' or 'ecliptic': Ecliptic coordinates
* 'sgal' or 'supergalactic': Super-Galactic coordinates

It would also be nice to allow users to define their own sky and non-sky
coordinate transformation.

Images
------

Plotting images as bitmaps or contours should be done via the usual
matplotlib methods (which can accept e.g. Numpy arrays or PIL Image objects
for RGB images)

    ax.imshow(fits.getdata('image.fits'))
    ax.imshow(Image.open('rgb_image.png'))
    ax.contour(fits.getdata('scuba.fits'))
    ax.contourf(fits.getdata('scuba.fits'))

Since WCS information technically doesn't require any information about
the image size, there would be no such thing as a 'mismatch' of dimensions
compared to the WCS.

In the case where we want to overplot an image with a different WCS
projection, we could do e.g.:

    ax[other_wcs].contour(scuba_image)

Note that ``ax[other_wcs].imshow()`` would not work (due to matplotlib
limitations).

Multi-dimensional WCS objects
-----------------------------

The WCS object can have more than two dimensions, but since we are ultimately
plotting two-dimensional data, we have to select which dimensions are used for
the x and y pixel coordinates in the axes, and we also have to select which
slices are selected for the dimensions that are not being shown. We can try
and encapsulate this in a single 'slice' argument which should be specified as
a list, and where each element is either an integer (indicating the slice
position) or a string containing 'x' or 'y', indicating the dimension that
should be used for the x and y pixel coordinates.

The user should set the slices to use when instantiating WCSAxes:

    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=WCS('cube.fits'),
                 slice=['x', 3, 'y'])

In the above example, the first and third dimension will be used as the x
and y pixel coordinates in the Axes, and the second pixel coordinate will be
set to 4 (could also be set to a floating-point value).

Info: we tried using separate dimensions= and slices= argument in APLpy (to
separate the selection of dimensions from the slices) but since the two are
closely related, an approach such as the above is simpler.

Allowing the slice to be changed on the fly will make it easy to develop
for example a utility to slice through cubes without having to
re-instantiate the ``WCSAxes`` object for each slice. To do this, we can use a
VariableSlice() object:

    v_slice = VariableSlice()

    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=WCS('cube.fits'),
                 slice=['x', v_slice, 'y'])

Subsequent calls to e.g.

    v_slice.set(3)

would then force the axes to refresh.

Coordinate grid
---------------

To overlay a grid, the normal matplotlib method can be used:

    ax['fk5'].grid()

and to modify the visual properties of the grid, users can capture a
handle to the grid object:

    g = ax['fk5'].grid()
    g.set_alpha(0.5)
    g.set_color('red')

Ticks and tick label properties
-------------------------------

While for many images, the coordinate axes are aligned with the pixel axes,
this is not always the case, especially in coordinate systems with high
curvature, where the coupling between x- and y-axis to actual coordinates
become less well-defined.

Therefore rather than referring to ``x`` and ``y`` ticks (as in APLpy), we
propose a new API, which allows very flexible customization of plots:

    # Access the first WCS coordinate
    ra = ax.coords[0]

    # Show ticks for this coordinate on the bottom and right axes
    ra.set_tick_location('br')

    # Show tick labels only for the bottom axis
    ra.set_ticklabel_location('b')

    # Set the format of the tick labels
    ra.set_ticklabel_format('hh:mm:ss')

    # Set the spacing of the ticks
    from astropy import units as u
    ra.set_tick_spacing(2. * u.arcmin)

and so on - we can also consider allowing aliases to certain coordinates, for
example:

    ra = ax.coords['RA--']

or even:

    ra = ax.coords['ra']

This interface can also be shown to plot only the grid lines for a specific
coordinate:

    ax.coords[2].grid()

This could be useful for example to show contours of constant wavelength in an
IFU cube with distortions.

Tick label format
-----------------

The format of the tick labels can be specified with:

    ra.set_ticklabel_format('x.xxx')  # decimal, non-angle coordinates,
                                      # 3 decimal places

    ra.set_ticklabel_format('d.ddddd')  # decimal degrees, 5 decimal places

    ra.set_ticklabel_format('dd:mm')  # sexagesimal, 1' precision
    ra.set_ticklabel_format('dd:mm:ss')  # sexagesimal, 1" precision
    ra.set_ticklabel_format('dd:mm:ss.ss')  # sexagesimal, 0.01" precision

    ra.set_ticklabel_format('hh:mm')  # sexagesimal (hours)
    ra.set_ticklabel_format('hh:mm:ss')  # sexagesimal (hours)
    ra.set_ticklabel_format('hh:mm:ss.ss')  # sexagesimal (hours)

Tick/label spacing
------------------

The spacing of ticks/tick labels should default to something sensible, but
users often want to be able to specify the spacing of the ticks. The tick
positions could be set manually:

    ra.set_ticks([242.2, 242.3, 242.4])

but users can also set the spacing explicitly:

    ra.set_ticks(0.1)

or can set the approximate number of ticks:

    ra.set_ticks(number=4)

In the case of angles, it is often convenient to specify the spacing as a
fraction of a degree - users should then use the units framework:

    from astropy import units as u
    ra.set_ticks(5. * u.arcmin)

This is to avoid roundoff errors.

Note: allowing users to specify their own tick position or spacing,
and the label format raises the issue of what happens when the value cannot
be well represented by the format (this confused a lot of users in APLpy
initially). For example, if one specifies the format to be 'dd:mm:ss' and a
tick spacing of 0.001 degrees then the ticks will go:

    13:42:00
    13:42:01
    13:42:03
    13:42:04
    13:42:06

and so on, which confuses users. In APLpy, we now raise an exception if the
tick spacing is not compatible with the labeling, but it might be best to
simply clip to the nearest value and emit an INFO message. If a tick spacing
is not specified, the minimum tick spacing should be set by the minimum one
representable by the format. For example, if the format is set to 'dd:mm',
then the default tick spacing cannot end up being smaller than 1' (and if it
does it gets reset to 1').

Patches/shapes/lines
--------------------

To overlay arbitrary matplotlib patches in pixel coordinates:

    from matplotlib.patches import Rectangle
    r = Rectangle((43., 44.), 23., 11.)
    ax['pixel'].add_patch(r)

To overlay arbitrary matplotlib patches in world coordinates:

    from matplotlib.patches import Rectangle
    r = Rectangle((272., -44.), 1., 0.8)
    ax['fk5'].add_patch(r)

And similarly:

    ax['gal'].add_collection(c)
    ax['fk4'].add_line(l)

The same applies to normal plotting methods such as scatter/plot/etc:

    ax['gal'].scatter(l, b)

Multiple coordinate systems
---------------------------

As shown above, the API could be set up to allow users to plot objects in
different coordinate systems. In general, using

    ax[coordinate_system]

has the effect of creating a new ``WCSAxes`` that is linked to the original one
but contains the additional coordinate transformation. By default, the axes
and labels are not shown, but the user could enable these, e.g.:

    glon = ax['gal'].coords[0]
    glon.show()

then adjust the parameters accordingly. By default, As a convenience, users
can simply use

    ax.switch_default_label_system('gal')

to hide all other labels

Offset systems
--------------

In order to switch to showing offset rather than absolute coordinates, one can
use

    ax['world'].enable_offset_mode(12.14, 15.55)

To revert to absolute coordinates:

    ax['world'].disable_offset_mode()

Examples
========

Example 1
---------

![](https://github.com/astrofrog/astropy-api/blob/wcs-plotting-proposal/wcs/pgsbox4.png?raw=true)

    import matplotlib.pyplot as plt
    from astropy.wcs import WCS, WCSAxes

    # Initialize figure
    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=WCS('image.fits'))
    fig.add_axes(ax)

    # Set format
    ax.set_tick_format(0, 'hh')
    ax.set_tick_format(1, 'dd:mm:ss')

    # Set tick label location
    ax.set_tick_location(0, 'lb')  left and bottom
    ax.set_tick_location(1, 'tr')  top and right

    # Draw grid
    ax.grid()

    # Save image
    fig.savefig('example1.png')

Example 2
---------

![](https://github.com/astrofrog/astropy-api/blob/wcs-plotting-proposal/wcs/pgsbox5.png?raw=true)

    import matplotlib.pyplot as plt
    from astropy.wcs import WCS, WCSAxes

    # Initialize figure
    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=WCS('image.fits'))
    fig.add_axes(ax)

    # Set format for the RA tick labels
    ax.set_tick_format(0, 'ddd')
    ax.set_tick_format(1, 'ddd')

    # Set tick label location
    ax.set_tick_location(0, 'bt')
    ax.set_tick_location(1, 'lr')

    # Draw grid for coordinates separately
    ax.grid(0, color='red')
    ax.grid(1, color='blue')

    # Set title
    ax.set_title("WCS conic equal area projection", color='turquoise')

    # Save image
    fig.savefig('example2.png')

Example 3
---------

![](https://github.com/astrofrog/astropy-api/blob/wcs-plotting-proposal/wcs/pgsbox6.png?raw=true)

    import matplotlib.pyplot as plt
    from astropy.wcs import WCS, WCSAxes

    # Initialize figure
    fig = plt.figure()
    ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=WCS('image.fits'))
    fig.add_axes(ax)

    # Set format for the RA tick labels
    ax.set_tick_format(0, 'hh')
    ax.set_tick_format(1, 'ddd')

    # Set tick label location
    ax.set_tick_location(0, 'btl')
    ax.set_tick_location(1, 'blr')

    # Draw grid for coordinates separately, and pass functions to determine
    the color rather than fixed colors.
    ax.grid(0, color=color_funtion_ha)
    ax.grid(1, color=color_funtion_dec)

    # Set title
    ax.set_title("WCS polyconic projection")

    # Save image
    fig.savefig('example3.png')

Example 4
---------

![](https://github.com/astrofrog/astropy-api/blob/wcs-plotting-proposal/wcs/pgsbox7.png?raw=true)

    import matplotlib.pyplot as plt
    from astropy.wcs import WCS, WCSAxes

    # Initialize figure
    fig = plt.figure()
    ax_gal = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=WCS('image.fits'))
    fig.add_axes(ax_gal)

    # Get an Axes with ecliptic coordinates
    ax_ecl = ax.get_axes('ecliptic')

    # Set format for the RA tick labels
    ax_gal.set_tick_format(0, 'ddd')
    ax_gal.set_tick_format(1, 'ddd')
    ax_gal.set_tick_color(0, 'green')
    ax_gal.set_tick_color(1, 'green')

    # Set tick label location
    ax_gal.set_tick_location(0, 't')
    ax_gal.set_tick_location(1, 'r')

    # Draw grid
    ax_gal.grid(color='green')

    # Set ticks for Galactic axes (defaults of x=longitude and y=latitude are
    fine here, so don't change)
    ax_ecl.set_tick_color(0, 'orange)
    ax_ecl.set_tick_color(1, 'blue)

    # Plot grid separately
    ax_ecl.grid(0, color='orange')
    ax_ecl.grid(1, color='blue')

    # Set title
    ax.set_title("WCS plate caree projection")

    # Save image
    fig.savefig('example4.png')
