from setuptools import setup

setup(name='scopebin',
      version='0.1',
      description='A package for viewing and analyzing microscope images in Python.',
      url='https://github.com/omidalam/scopebin',
      author='Omid Gholamalamdari',
      author_email='omid.gholamalamdari@mcgill.ca',
      license='GPL-3.0 license',
      packages=['scopebin'],
      install_requires=[
          'scikit-image',
          'numpy',
          'scipy',
          'matplotlib',
          'python-bioformats',
          'reportlab',
          'PyPDF2',
          'pandas',
      ],      
      zip_safe=False)
