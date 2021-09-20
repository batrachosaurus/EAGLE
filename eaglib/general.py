# TODO: reduce this module
# This code can have only standard Python imports
from collections import OrderedDict


def filter_list_(in_list):  # This can be reduced with 'list(filter(lambda li: li.strip(), in_list))'
    filtered_list = list()
    for elm_ in in_list:
        elm = None
        elm = elm_.strip()
        if elm:
            filtered_list.append(elm)
    return filtered_list


def revert_dict_(in_dict):  # This can be reduced with '{v: k for k, v in in_dict}'
    out_dict = OrderedDict()
    for key in in_dict:
        out_dict[in_dict[key]] = key
    return out_dict
