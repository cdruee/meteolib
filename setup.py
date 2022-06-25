from setuptools import setup, find_packages


setup(
    name='meteolib',
    version='0.0.1',
    license='EUPL-1.1',
    author="Clemens Druee",
    author_email='druee@uni-trier.de',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    url='https://github.com/cdruee/meteolib',
    keywords='meteorology',
    install_requires=[
          'numpy', 'pandas'
      ],
    long_description='''
    This module conatins standard equations, constants and conversions
    for general use in meteorology. Values and methods herin are are either
    adopted or recommended by the World meteorological Organization (WMO)
    and similar bodies or backed by peer-reviewed literature and fully citable.

    Note: This module is in a stub state until publication of the corresponding paper
    '''
)