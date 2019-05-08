def detect_frame_parallel(orf_c, ovorf_c, orf_ori):
    if int(float(orf_c[0] - ovorf_c[0]) % 3.0) == 0:
        return +1
    if orf_ori == "+" or orf_ori == 1:
        if int(float(orf_c[0] - (ovorf_c[0] - 1)) % 3.0) == 0:
            return +2
        if int(float(orf_c[0] - (ovorf_c[0] - 2)) % 3.0) == 0:
            return +3
    else:
        if int(float(orf_c[0] - (ovorf_c[0] + 1)) % 3.0) == 0:
            return +2
        if int(float(orf_c[0] - (ovorf_c[0] + 2)) % 3.0) == 0:
            return +3


def detect_frame_antiparallel(orf_c, ovorf_c, orf_ori):
    if int(float(orf_c[0] - ovorf_c[0]) % 3.0) == 0:
        return -1
    if orf_ori == "+" or orf_ori == 1:
        if int(float(orf_c[0] - (ovorf_c[0] + 1)) % 3.0) == 0:
            return -2
        if int(float(orf_c[0] - (ovorf_c[0] + 2)) % 3.0) == 0:
            return -3
    else:
        if int(float(orf_c[0] - (ovorf_c[0] - 1)) % 3.0) == 0:
            return -2
        if int(float(orf_c[0] - (ovorf_c[0] - 2)) % 3.0) == 0:
            return -3


def get_overlap_coords(orf_c, ovorf_c):
    return min(orf_c[1], ovorf_c[1]), max(orf_c[0], ovorf_c[0])


def get_overlap_length(orf_c, ovorf_c):
    overlap_coords = get_overlap_coords(orf_c=orf_c, ovorf_c=ovorf_c)
    return overlap_coords[1] - overlap_coords[0] + 1
