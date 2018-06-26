from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

ext_modules = [
    Extension("illstack.cyprof"
              ,["illstack/cyprof.pyx"],
              include_dirs=[np.get_include()]
    ),]

setup(name='illstack',
      ext_modules=cythonize(ext_modules),
      include_dirs=[np.get_include()],
      version='0.1',
      description='for converting halo catalogs to maps for a given hod',
      url='http://github.com/marcelo-alvarez/illstack',
      author='Marcelo Alvarez and Nick Battaglia',
      author_email='marcelo.alvarez@berkeley.edu',
      license='GNU General Public License v3.0',
      packages=['illstack'],
      zip_safe=False)

