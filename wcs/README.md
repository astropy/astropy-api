Example 1
=========

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
    ax.set_tick_location(0, 'lb')  # left and bottom
    ax.set_tick_location(1, 'tr')  # top and right

    # Draw grid
    ax.grid()

    # Save image
    fig.savefig('example1.png')

Example 2
=========

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
=========

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
    # the color rather than fixed colors.
    ax.grid(0, color=color_funtion_ha)
    ax.grid(1, color=color_funtion_dec)

    # Set title
    ax.set_title("WCS polyconic projection")

    # Save image
    fig.savefig('example3.png')

Example 4
=========

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
    # fine here, so don't change)
    ax_ecl.set_tick_color(0, 'orange)
    ax_ecl.set_tick_color(1, 'blue)

    # Plot grid separately
    ax_ecl.grid(0, color='orange')
    ax_ecl.grid(1, color='blue')

    # Set title
    ax.set_title("WCS plate caree projection")

    # Save image
    fig.savefig('example4.png')
