from setuptools import setup, find_packages

setup(name='EAGLE',
      version='0.7.3',
      packages=find_packages(),
      package_data={'EAGLE': ['configs/*', ],
                    'EAGLEdb': ['org_tables/*', ]},
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
