FROM python:2

RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    muscle \
    emboss \
    embassy-phylip \
    hmmer 

ADD fastme /bin/
ADD EAGLE-0.20.6.tar.gz /

RUN pip install /EAGLE-0.20.6.tar.gz
