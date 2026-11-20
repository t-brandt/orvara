"""
Build script for orvara's Cython extension.

All project metadata lives in pyproject.toml; this file exists only because the
``orvara.orbit`` extension module has to be cythonised and compiled.
"""
import os
from setuptools import setup, Extension

# Cython is a build-time requirement declared in pyproject.toml, so it is always
# importable here when pip/uv builds the project.
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = ['-O3']
if os.name == 'nt':                       # MSVC does not understand -O3
    EXTRA_COMPILE_ARGS = ['/O2']

extensions = [
    Extension(
        'orvara.orbit',
        sources=[os.path.join('orvara', 'orbit.pyx')],
        extra_compile_args=EXTRA_COMPILE_ARGS,
    )
]

setup(
    # Only language_level is set globally.  orbit.pyx applies
    # @cython.boundscheck(False) / @cython.wraparound(False) per function,
    # where the author established they are safe -- turning them on (or
    # cdivision) for the whole module changes the fitted results.
    ext_modules=cythonize(
        extensions,
        language_level='3',
        compiler_directives={'language_level': '3'},
    ),
)
