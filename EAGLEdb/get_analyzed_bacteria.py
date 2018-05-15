import os
import sys
import json


def are_bacteria_analyzed(bacteria_list_path=None, analyzed_bacteria_path=None):
    if not bacteria_list_path and not analyzed_bacteria_path:
        try:
            bacteria_list_path = sys.argv[1]
            analyzed_bacteria_path = sys.argv[2]
        except IndexError:
            print("Number of arguments must be 2: 1 - bacteria list json path; 2 - analyzed bacteria json path")
    bacteria_list_f = open(bacteria_list_path)
    bacteria_list = json.load(bacteria_list_f)
    bacteria_list_f.close()
    analyzed_bacteria = dict()
    for bacterium in bacteria_list:
        if os.path.exists(bacterium[u'16S_rRNA_file']):
            analyzed_bacteria[bacterium[u'strain'].replace("_", " ")] = True
    analyzed_bacteria_f = open(analyzed_bacteria_path, 'w')
    json.dump(analyzed_bacteria, analyzed_bacteria_f)
    analyzed_bacteria_f.close()
