from setuptools import setup, find_packages

setup(name='EAGLE',
      version='0.9.10',
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
                  "EAGLEdb.prepare_ncbi_summary = EAGLEdb.files_utils:prepare_summary_table",
                  "EAGLEdb.get_analyzed_bacteria = EAGLEdb.files_utils:are_bacteria_analyzed",
                  "EAGLEdb.join_bacteria_lists = EAGLEdb.files_utils:join_bacteria_list_files",
                  "EAGLEdb = EAGLEdb.__main__:main",
            ]
      })
