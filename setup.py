from setuptools import setup, find_packages

setup(name='EAGLE',
      version='0.8.0',
      packages=find_packages(),
      package_data={'EAGLE': ['configs/*', ],
                    'EAGLEdb': ['org_tables/*', ]},
      install_requires=[
            'wget',
            'pyaml',
            'pandas',
            'redis'
      ],
      entry_points={
            'console_scripts': [
                  "EAGLEdb.prepare_ncbi_summary = EAGLEdb.prepare_ncbi_summary:prepare_summary_table",
                  "EAGLEdb = EAGLEdb.__main__:main",
            ]
      })
