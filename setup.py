from setuptools import setup, find_packages

setup(name='EAGLE',
      version='0.5.0',
      packages=find_packages(),
      package_data={'EAGLE': ['configs/*', ]},
      install_requires=[
            'wget',
            'pyaml',
            'pandas'
      ])
