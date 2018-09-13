#!/usr/bin/env python

"""
setup.py file for splashpipe
"""

from distutils.core import setup, Extension
import numpy

splashback_module = Extension('_splashback',
                          sources=['src/splashback.i', 'src/splashback.cpp', 'src/cosmology.cpp', 'src/haloes.cpp', 'src/powerspectrum.cpp', 'src/gauleg.cpp', 'src/kdtree2.cpp'],
                          swig_opts=["-keyword","-c++"],
                          libraries=['m','gsl','gslcblas', 'boost_system'],
                          include_dirs=[numpy.get_include()],
                          )

setup (name        = 'splashpipe',
      version     = '0.01',
      author      = "Surhud More",
      url         = "http://member.ipmu.jp/surhud.more/research",
      author_email= "surhud.more@ipmu.jp",
      description = """Splashback radius pipeline tools""",
      ext_modules = [splashback_module],
      license     = ['GPL'],
      py_modules  = ["splashback"],
      package_dir = { '':'src'},
      )
