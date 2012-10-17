#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import astropy.units as u

""" The Quantity class will represent a number + unit + uncertainty """

# -------------------
# Creating quantities
# -------------------

# One method for creating a quantity is through operations with Unit objects:
q = 11.42 * u.meter # returns a Quantity object
q = 11.42 / u.meter
q = 182.234 + u.meter # raises an error

# You could also use the Quantity constructor, but 'unit' must be supplied:
q1 = u.Quantity(11.412, unit=u.meter)
q2 = u.Quantity(21.52, "cm")
q3 = u.Quantity(11.412) # raises an exception because no unit specified

# A Quantity object may also have an uncertainty associated with them, but this will be 
#   implemented at a later date, probably in a subclass of Quantity. The usage might look
#   something like this:
q1 = u.Quantity(14.6, u.meter, uncertainty=some_probability_distribution)
q2 = u.Quantity(14.6, u.meter, uncertainty=(0.1, 0.2))
q1 = u.Quantity(11.41, u.meter, uncertainty=0.01) 
q2 = u.Quantity(14.6, u.meter, uncertainty=0.2)

# We will also support a number of astropy Uncertainty classes, but the names have not been
#   decided upon. That will look like:
q1 = u.Quantity(14.6, u.meter, uncertainty=StandardDeviationUncertainty(0.5))
q2 = u.Quantity(14.6, u.meter, uncertainty=VarianceUncertainty(0.25))
# etc.

# ----------
# Operations
# ----------

# Each of these operations would return a new Quantity object
q1 = u.Quantity(11.41, u.meter)
q2 = u.Quantity(15.37, u.meter)
new_quantity = q1 + q2
new_quantity = q1 - q2
new_quantity = q1 * q2
new_quantity = q1 / q2
new_quantity = -q1

# Operations will *adopt the units from the object on the left* if the units are equivalent
#   for both objects involved, and may raise an exception otherwise:
q1 = u.Quantity(11.41, u.meter)
q2 = u.Quantity(15.37, u.kilometer)
new_quantity = q1 + q2 # has unit u.meter!

q1 = u.Quantity(11.41, u.meter)
q2 = u.Quantity(15.37, u.second)
q1 / q2 # valid, returns an object in meters per second
q1 + q2 # raises an exception

# For objects with uncertainties, we can't propagate errors unless the user uses one of the astropy
#   Uncertainty classes, which know how to do their own error propagation (but are still in development).

# In the future, something like this will work (names may change):
q1 = u.Quantity(11.41, u.meter, uncertainty=StandardDeviationUncertainty(0.35))
q2 = u.Quantity(15.37, u.meter, uncertainty=StandardDeviationUncertainty(0.11))
new_quantity = q1 + q2 # error propagation done for the user

# The Units package will have a native type that represents "no unit." Operations with such dimensionless
#   objects will look like this (replace unit="" with whatever is decided for Unit package):
Quantity(15.1234, unit=u.kilogram) * Quantity(0.75, unit="") # should work
Quantity(15.1234, unit=u.kilogram) / Quantity(0.75, unit="") # should work

Quantity(15.1234, unit=u.kilogram) + Quantity(0.75, unit="") # should NOT work
Quantity(15.1234, unit=u.kilogram) - Quantity(0.75, unit="") # should NOT work

# ----------------
# Converting units
# ----------------

# to get a new quantity object in some new, equivalent units, you would use the to() method:
q = 11.412*u.meter
q.to(u.kilometer) # this returns a new quantity object with the supplied units
q.to(u.zettastokes) # Fails because not equivalent

# you may also just want the *value* of the quantity in some units. To do this, use the value attribute:
q1 = Quantity(194.118, unit=u.kilometer)
q1.value # returns 194.118
q1.to(u.centimeter).value # returns the *value* (just a number, not an object) of this quantity in the specified units

# --------
# Equality
# --------
Quantity(1000, unit=u.m) == Quantity(1, units=u.km) # returns True
Quantity(1, units=u.m) == Quantity(1, units=u.km) # returns False
Quantity(1, units=u.m) == Quantity(1, units=u.s) # raises IncompatibleUnitError()

# ----------
# Displaying
# ----------
# Printing the quantity should display the number and the units when converting to a string
str(quantity) # returns '11.412 m'
repr(quantity) # returns '<Quantity Value: 11.412 Units: meter>' or something