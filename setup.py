from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_module = Extension("orbit",
                       ['orbit.pyx'],
                       extra_compile_args=['-O3'],
)

setup(
    name='orbit',
    cmdclass={'build_ext': build_ext},
    ext_modules=[ext_module],
    version='0.1.0',
    python_requires='>=3.5',
    install_requires=['numpy>=1.13', 'htof', 'emcee', 'Cython', 'pandas', 'astropy'],
)

