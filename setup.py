from setuptools import setup, find_packages

setup(name='EAGLE',
      version='0.6.0',
      packages=find_packages(),
      package_data={'EAGLE': ['configs/*', ]},
      install_requires=[
            'wget',
            'pyaml',
            'pandas'
      ],
      entry_points={
            'console_scripts': [
                  "EAGLEdb.prepare_ncbi_summary = EAGLEdb.prepare_ncbi_summary:prepare_summary_table",
            ]
      })
