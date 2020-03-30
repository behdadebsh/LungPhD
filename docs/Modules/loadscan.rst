============
load_scan.py
============

This script is used to read the stack of DICOM images. CT images are loaded into Python via Pydicom library. There is more documentation on how to install this package on this `link <https://pydicom.github.io/pydicom/stable/getting_started.html#installing-pydicom>`_. You also need to have `NumPy <http://www.numpy.org/>`_ Python package which you can also add it to your python packages using `pip <https://pypi.org/project/pip/>`_.

This script is written using two functions:
- load_scan
- get_pixels_hu
The one named "load_scan" requires a path to the directory that the DICOM images are stored and the second function, "get_pixels_hu", converts everything to hounsfield unit (HU).
