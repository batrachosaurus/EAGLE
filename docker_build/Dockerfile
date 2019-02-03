FROM python:2

ADD embassy-phylip_3.69.660-3_amd64.deb /src/

RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    muscle \
    emboss \
    hmmer \
    /src/embassy-phylip_3.69.660-3_amd64.deb

ADD fastme /bin/
ADD EAGLE-*.tar.gz /src/

RUN pip install /src/EAGLE-*

RUN mkdir /EAGLE
CMD cd /EAGLE

