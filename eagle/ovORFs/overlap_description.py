# coordinates are ALWAYS in forward strand: c[0] < c[1]

def detect_frame_parallel(orf_c, ovorf_c, orf_ori):
    if (orf_c[0] - ovorf_c[0]) % 3 == 0:
        return +1
    if orf_ori == "+" or orf_ori > 0:
        if (orf_c[0] - (ovorf_c[0] - 1)) % 3 == 0:
            return +2
        if (orf_c[0] - (ovorf_c[0] - 2)) % 3 == 0:
            return +3
    else:
        if (orf_c[0] - (ovorf_c[0] + 1)) % 3 == 0:
            return +2
        if (orf_c[0] - (ovorf_c[0] + 2)) % 3 == 0:
            return +3


def detect_frame_antiparallel(orf_c, ovorf_c, orf_ori):
    if (orf_c[0] - ovorf_c[0]) % 3 == 0:
        return -1
    if orf_ori == "+" or orf_ori > 0:
        if (orf_c[0] - (ovorf_c[0] + 1)) % 3 == 0:
            return -2
        if (orf_c[0] - (ovorf_c[0] + 2)) % 3 == 0:
            return -3
    else:
        if (orf_c[0] - (ovorf_c[0] - 1)) % 3 == 0:
            return -2
        if (orf_c[0] - (ovorf_c[0] - 2)) % 3 == 0:
            return -3


def get_overlap_coords(orf_c, ovorf_c):
    return max(orf_c[0], ovorf_c[0]), min(orf_c[1], ovorf_c[1])


def get_overlap_length(orf_c, ovorf_c):
    overlap_coords = get_overlap_coords(orf_c=orf_c, ovorf_c=ovorf_c)
    return overlap_coords[1] - overlap_coords[0] + 1
