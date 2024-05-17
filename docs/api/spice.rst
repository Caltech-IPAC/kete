spice
=====

This is a thread-safe, read only re-implementation of a SPICE kernel interpreter.
Outputs of this exactly match the common cSPICE interpreter, but can be easily
used among an arbitrary number of cores. SPICE kernel files are loaded directly
into RAM.

.. note::

   This does not use cSPICE, or any original SPICE library code.

   cSPICE is difficult to use in a thread safe manner which is limiting when
   performing orbit calculations on millions of objects.

Data files which are automatically included:

DE440 - A SPICE file containing the planets within a few hundred years.

BSP files are also automatically included for the 5 largest asteroids, which are
used for numerical integrations when the correct flags are set.

PCK Files which enable coordinate transformations between Earths surface and the
common inertial frames.

.. automodule:: neospy.spice
   :members:
   :inherited-members:
