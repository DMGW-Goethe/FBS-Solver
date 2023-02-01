from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import os
import numpy as np

pyfbs_file = "pyfbs.pyx"

root_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
include_folder = os.path.join(root_folder, "include")
pyfbs_folder = os.path.join(root_folder, "pyfbs")

liblist = ["fbs"]


pyfbs_ext = Extension("pyfbs_cython",
                        [os.path.join(pyfbs_folder, pyfbs_file)],
                        include_dirs=[np.get_include(), include_folder],
                        language="c++",
                        libraries=liblist,
                        library_dirs=[root_folder],
                        extra_compile_args = ["-O3", "-ffast-math", "-fopenmp" ],
                        extra_link_args=['-lgomp']
                    )

import six
pyfbs_ext.cython_directives = {'language_level': "3" if six.PY3 else "2"}

setup(
    name='pyfbs',
    cmdclass={'build_ext':build_ext},
    ext_modules=[pyfbs_ext],
)


