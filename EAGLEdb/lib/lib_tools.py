def get_links_from_html(html_file):
    links = []
    for lines in html_file:
        line = None
        line = lines.strip()
        if not line: continue
        if "<a href" not in line: continue
        line_list = line.split("<a href")
        links.append(line_list[1].split('">')[0].strip(' ="/'))
    return links
