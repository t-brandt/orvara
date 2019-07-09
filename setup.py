from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

source_files = ['orbit3d/orbit.pyx']
modules = [Extension("orbit3d.orbit", source_files, extra_compile_args=['-O3'])]

setup(name='orbit3d',
      ext_modules=cythonize(modules),
      version='0.1.0',
      python_requires='>=3.5',
      install_requires=['numpy>=1.13', 'htof', 'emcee',
                        'Cython', 'pandas', 'astropy', 'pytest'])
