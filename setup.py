from setuptools import setup, find_packages

setup(name='EAGLE',
      version='0.8.4',
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
                  "EAGLEdb.get_analyzed_bacteria = EAGLEdb.get_analyzed_bacteria:are_bacteria_analyzed",
                  "EAGLEdb = EAGLEdb.__main__:main",
            ]
      })
