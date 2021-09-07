FROM python:3.9-slim-buster as base

FROM base as builder

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -y update
RUN apt-get install -y git libgdal-dev build-essential
RUN apt-get install -y --reinstall g++

ENV PATH /code/opendrift/opendrift/scripts:$PATH

RUN mkdir /code
WORKDIR /code

RUN pip install poetry
RUN poetry config virtualenvs.create false

# Cache cartopy maps
# RUN /bin/bash -c "echo -e \"import cartopy\nfor s in ('c', 'l', 'i', 'h', 'f'): cartopy.io.shapereader.gshhs(s)\" | python"

# Install opendrift
ADD . /code
RUN poetry install

# FROM base
# COPY --from=builder /usr/local /usr/local
# COPY --from=builder /root /root
# ADD . /code
# WORKDIR /code

# RUN apt-get -y update
# RUN apt-get install -y libgdal20 libproj12

# Test installation
# RUN /bin/bash -c "echo -e \"import opendrift\" | python"

