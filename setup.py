# Cythonizes the C library
# Run by the makefile

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

examples_extension = Extension(
    name="spherical",
    sources=["spherical.pyx"],
    libraries=["spherical"],
    library_dirs=["lib"],
    include_dirs=["lib"]
)
setup(
    name="spherical",
    ext_modules=cythonize([examples_extension])
)