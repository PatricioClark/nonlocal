# Cythonizes the C library
# Run by the makefile

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

examples_extension = Extension(
    name="spherical4",
    sources=["spherical.pyx"],
    libraries=["spherical4"],
    library_dirs=["lib"],
    include_dirs=["lib"]
)
setup(
    name="spherical4",
    ext_modules=cythonize([examples_extension])
)
