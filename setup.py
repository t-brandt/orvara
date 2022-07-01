from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

source_files = ['orvara/orbit.pyx']
modules = [Extension("orvara.orbit", source_files, extra_compile_args=['-O3'])]

setup(name='orvara',
      ext_modules=cythonize(modules),
      version='1.1.2',
      python_requires='>=3.5',
      package_dir={'orvara': 'orvara'},
      install_requires=['numpy>=1.13', 'htof>=1.1.5', 'emcee', 'ptemcee',
                        'Cython', 'pandas', 'astropy', 'pytest', 'mock'],
      entry_points={'console_scripts': ['fit_orbit=orvara.main:run', 'plot_orbit=orvara.main_plotting:run']})
