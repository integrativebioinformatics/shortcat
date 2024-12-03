FROM mambaorg/micromamba:1.4.9

LABEL image.author.name "Adolfo Rojas Hidalgo"
LABEL image.author.email "adolfo.rojas@ug.uchile.cl"

COPY --chown=$MAMBA_USER:$MAMBA_USER ShortCat.yml /tmp/ShortCat.yml

RUN micromamba --version

RUN micromamba create -n ShortCat 

RUN micromamba install -y -n ShortCat -f /tmp/ShortCat.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/ShortCat/bin:$PATH
RUN rm /tmp/ShortCat.yml
USER root
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*