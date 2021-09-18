from setuptools import setup, find_packages

setup(name='EAGLE',
      version='0.0.1',
      author="Denis Moshensky",
      author_email="loven-doo@fbb.msu.ru", 
      description="Essential and Advantageous Genes Location Explorer",
      url="https://github.com/loven-doo/EAGLE",
      classifiers=[
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
      ],
      packages=find_packages(exclude=find_packages("Research") +
                                     ["eagle.tests", "eaglib.tests", "eagledb.tests", "Research"]),
      package_data={'eagle': ['tests/*'],
                    'eaglib': ['tests/*'],
                    'eaglib.logging': ['eagle_log_conf.yml'],
                    'eagledb': ['org_tables/*', 'tests/*']},
      install_requires=[
          'wget >= 3.2',
          'pyaml >= 3.12',
          'numpy >= 1.14.3',
          'pandas >= 0.22.0',
          'matplotlib >= 2.2.3',
          'scipy >= 1.1.0',
          'biopython >= 1.72',
          'DendroPy >= 4.4.0',
          'psutil >= 5.6.1',
          'jsondler >= 0.0.4',
          'Deprecated >= 1.2.9',
      ],
      entry_points={
          'console_scripts': [
              "eagle_db.prepare_ncbi_summary = eagledb.files_utils:prepare_summary_table",
              "eagle_db.get_analyzed_bacteria = eagledb.files_utils:are_bacteria_analyzed",
              "eagle_db.join_bacteria_lists = eagledb.files_utils:join_bacteria_list_files",
              "eagle-db.create = eagledb.creation:create_cmd",
              "eagle-db.bacteria.create-ncbi = eagledb.bacteria.scenarios:create_ncbi_cmd",
              "eagle = eagle.__main__:main",
              "eagle.explore_orfs = eagle.orfs_explorer:explore_orfs_cmd",
              "eagle.classify_orfs = eagle.orfs_classifier:classify_orfs_cmd",
              "eagle.btax_name = eagle.btax_scanner:get_btax_name_cmd",
          ]
      })
