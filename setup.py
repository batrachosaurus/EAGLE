from setuptools import setup, find_packages

setup(name='EAGLE',
      version='0.21.11',
      author="Denis Moshensky",
      author_email="loven-doo@fbb.msu.ru", 
      description="Essential and Advantageous Genes Location Explorer",
      url="https://github.com/loven-doo/EAGLE",
      classifiers=[
          "Programming Language :: Python :: 2.7",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
      ],
      packages=find_packages(exclude=find_packages("Research").__add__(["eagle.tests", "eagledb.tests", "Research"])),
      package_data={'eagle': ['configs/*', 'tests/*'],
                    'eagledb': ['org_tables/*', 'tests/*']},
      install_requires=[
          'wget >= 3.2',
          'pyaml >= 3.12',
          'numpy >= 1.14.3',
          'pandas == 0.22.0',
          'matplotlib >= 2.2.3',
          'scipy >= 1.1.0',
          'biopython >= 1.72',
          'DendroPy >= 4.4.0',
          'redis >= 2.10.6',
          'psutil >= 5.6.1',
      ],
      entry_points={
          'console_scripts': [
              "eagle_db.prepare_ncbi_summary = eagledb.files_utils:prepare_summary_table",
              "eagle_db.get_analyzed_bacteria = eagledb.files_utils:are_bacteria_analyzed",
              "eagle_db.join_bacteria_lists = eagledb.files_utils:join_bacteria_list_files",
              "eagle_db = eagledb.__main__:main",
              "eagle = eagle.__main__:main",
              "eagle.explore_orfs = eagle.orfs_explorer:explore_orfs_cmd",
              "eagle.classify_orfs = eagle.orfs_classifier:classify_orfs_cmd",
          ]
      })
