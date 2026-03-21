#!/usr/bin/env python

from setuptools import setup
from qtlseq.__init__ import __version__


setup(
    name='qtlseq',
    version='{}'.format(__version__),
    description='QTL-seq: pipeline to identify causative mutations responsible for a phenotype',
    author='Yu Sugihara',
    author_email='sugihara.yu.85s@kyoto-u.jp',
    url='https://github.com/YuSugihara/QTL-seq',
    license='GPL',
    packages=['qtlseq'],
    python_requires='>=3.8',
    install_requires=[
        'matplotlib',
        'numpy',
        'pandas',
    ],
    entry_points={
        'console_scripts': [
            'qtlseq = qtlseq.qtlseq:main',
            'qtlplot = qtlseq.qtlplot:main',
        ]
    },
)
