import random
import re
import string


def get_un_fix(un_num, fix_len):
    un_codes = ["_", '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E']
    # 'N' - undefined (num duplicates is bigger than len(un_codes))
    if fix_len == 1:
        try:
            return un_codes[un_num]
        except IndexError:
            return 'N'
    elif fix_len == 0:
        return ""
    elif un_num < len(un_codes):
        return un_codes[0] + get_un_fix(un_num, fix_len - 1)
    else:
        filled_rank = len(un_codes)**(fix_len-1)
        return un_codes[un_num//filled_rank] + get_un_fix(un_num % filled_rank, fix_len - 1)


def generate_random_string(l=10):
    return "".join(random.choice(string.ascii_letters + string.digits) for i in range(l))


def fullmatch_regexp_list(pattern, target_list):
    return list(map(lambda x: re.fullmatch(pattern, x), target_list))
