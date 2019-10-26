from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

source_files = ['orbit3d/orbit.pyx']
modules = [Extension("orbit3d.orbit", source_files, extra_compile_args=['-O3'])]

setup(name='orbit3d',
      ext_modules=cythonize(modules),
      version='0.1.2',
      python_requires='>=3.5',
      package_dir={'orbit3d': 'orbit3d'},
      package_data={'orbit3d': ['data/*.fits']},
      install_requires=['numpy>=1.13', 'htof>=0.1.3', 'emcee',
                        'Cython', 'pandas', 'astropy', 'pytest'],
      entry_points={'console_scripts': ['fit_orbit=orbit3d.main:run' ,'plot_orbit=orbit3d.main_plotting:run' ]})
