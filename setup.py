from setuptools import setup, find_packages

setup(name='EAGLE',
      version='0.2.1',
      packages=find_packages(),
      install_requires=[
            'wget',
            'pyaml'
      ])
