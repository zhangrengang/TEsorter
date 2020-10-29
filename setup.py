#!/usr/bin/env python
from TEsorter.version import __version__

from setuptools import setup, find_packages

with open('README.md') as f:
    long_description = f.read()


setup(
    name='TEsorter',
    version=__version__,
    description='Lineage-level classification of transposable elements using conserved protein domains',
    url='https://github.com/zhangrengang/TEsorter',
    download_url='https://github.com/zhangrengang/TEsorter/archive/v1.2.6.tar.gz',
    author='Zhang, Ren-Gang and Wang, Zhao-Xuan, Jacques Dainat',
    license='GPL-3.0',

    python_requires='>=3',
    packages=find_packages(),
    include_package_data=True,
    scripts=['TEsorter/modules/get_record.py',
            'TEsorter/modules/concatenate_domains.py',
            'TEsorter/modules/RepeatMasker.py',
            'TEsorter/modules/LTR_retriever.py',],
    entry_points={
        'console_scripts': ['TEsorter = TEsorter.app:main',
        'TEsorter-test = TEsorter.test.test_app:main',
        ],
    }
)
