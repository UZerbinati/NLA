from setuptools import setup

setup(name='NLA',
      version='0.0.1',
      description='Numeric Linear Algebra Library',
      url='http://uzerbinati.github.io',
      author='Umberto Zerbinati',
      author_email='uzerbinati@me.com',
      license='MIT',
      packages=['NLA'],
      install_requires=[
          'numpy',
          'mpmath',
          'matplotlib.pyplot'
      ],
      zip_safe=False)
