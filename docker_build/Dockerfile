FROM python:3

ADD fastme /bin/
ADD msaprobs /bin/
ADD KaKs_Calculator /bin/
ADD embassy-phylip_3.69.660-3_amd64.deb /src/

RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    muscle \
    mafft \
    emboss \
    hmmer \
    /src/embassy-phylip_3.69.660-3_amd64.deb
RUN pip install eagle && mkdir /EAGLE

WORKDIR /EAGLE

