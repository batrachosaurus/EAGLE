import multiprocessing as mp
import urllib2

from EAGLE.lib import worker
from EAGLE.constants import EAGLE_logger


def get_links_from_html(html_link, num_threads=1, n_tries=3, debug=False):
    n_t = 0
    links = mp.Manager().dict()
    params_list = []
    while n_t < n_tries:
        html_file = urllib2.urlopen(html_link)
        params_list.append({'function': _read_html_file_links,
                            'html_file': html_file.read().split("\n"),
                            'links': links,
                            'debug': debug})
        n_t += 1
    pool = mp.Pool(num_threads)
    pool.map(worker, params_list)
    pool.close()
    pool.join()
    links_list = links.keys()
    links_list.sort()
    return links_list


def _read_html_file_links(html_file, links, **kwargs):
    for lines in html_file:
        line = None
        line = lines.strip()
        if not line: continue
        if "<a href" not in line: continue
        if "parent directory" in line.lower() or ".." in line: continue
        line_list = line.split("<a href")
        links[(line_list[1].split('">')[0].strip(' ="/'))] = True
