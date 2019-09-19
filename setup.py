# Cythonizes the C library
# Run by the makefile

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

num = 2
examples_extension = Extension(
    name="spherical{}".format(num),
    sources=["spherical.pyx"],
    libraries=["spherical"],
    library_dirs=["lib{}".format(num)],
    include_dirs=["lib{}".format(num)]
)
setup(
    name="spherical{}".format(num),
    ext_modules=cythonize([examples_extension])
)
