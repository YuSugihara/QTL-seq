#!/usr/bin/env python

from distutils.core import setup
from qtlseq.__init__ import __version__

setup(name='qtlseq',
      version='{}'.format(__version__),
      description='QTL-seq: pipeline to identify causative mutations responsible for a phenotype',
      author='Yu Sugihara',
      author_email='yu57th@gmail.com',
      url='https://github.com/YuSugihara/QTL-seq',
      license='GPL',
      packages=['qtlseq'],
      entry_points={'console_scripts': [
            'qtlseq = qtlseq.qtlseq:main',
            'qtlplot = qtlseq.qtlplot:main',
            ]
        }
    )
