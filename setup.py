from setuptools import setup, find_packages


setup(
    name='meteolib',
    version='0.0.0',
    license='EUPL-1.1',
    author="Clemens Drüe",
    author_email='druee@uni-trier.de',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    url='https://github.com/cdruee/meteolib',
    keywords='meteorology',
    install_requires=[
          'numpy', 'pandas'
      ],

)