from setuptools import setup, find_packages

setup(name='EAGLE',
      version='0.13.2',
      packages=find_packages(exclude=find_packages("Research").__add__(["EAGLE.tests", "EAGLEdb.tests", "Research"])),
      package_data={'EAGLE': ['configs/*', 'tests/*'],
                    'EAGLEdb': ['org_tables/*', 'tests/*']},
      install_requires=[
            'wget >= 3.2',
            'pyaml >= 3.12',
            'numpy >= 1.14.3',
            'pandas == 0.22.0',
            'scipy >= 1.1.0',
            'biopython >= 1.72',
            'redis >= 2.10.6',
      ],
      entry_points={
            'console_scripts': [
                  "EAGLEdb.prepare_ncbi_summary = EAGLEdb.files_utils:prepare_summary_table",
                  "EAGLEdb.get_analyzed_bacteria = EAGLEdb.files_utils:are_bacteria_analyzed",
                  "EAGLEdb.join_bacteria_lists = EAGLEdb.files_utils:join_bacteria_list_files",
                  "EAGLEdb = EAGLEdb.__main__:main",
                  "EAGLE = EAGLE.__main__:main"
            ]
      })
