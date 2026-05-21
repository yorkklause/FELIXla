# Dockerfile
# designed by Kai, implemented by codex

FROM wzhou88/saigetractor:1.4.9

LABEL org.opencontainers.image.title="saigetractor-tractor-hybrid"
LABEL org.opencontainers.image.description="SAIGE-TRACTOR 1.4.9 with tractor_hybrid step2 input and hybrid format utilities"

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      g++ \
      libhts-dev \
      make \
      pkg-config && \
    rm -rf /var/lib/apt/lists/*

COPY saige_step2_adapter/TractorHybrid.hpp /app/src/TractorHybrid.hpp
COPY saige_step2_adapter/TractorHybrid.cpp /app/src/TractorHybrid.cpp
COPY docker/saigetractor-hybrid/patch_saigetractor_hybrid.py /tmp/patch_saigetractor_hybrid.py
COPY Makefile /scripts/Makefile
COPY src /scripts/src

RUN /app/.pixi/envs/default/bin/python /tmp/patch_saigetractor_hybrid.py && \
    export PATH=/app/.pixi/envs/default/bin:$PATH && \
    cd /app && \
    Rscript -e 'Rcpp::compileAttributes("/app")' && \
    R CMD INSTALL . && \
    cd /scripts && \
    make all CXX=g++ CXXFLAGS="-O3 -std=c++17 -Wall -Wextra -Wpedantic" BIN_DIR=/scripts/bin && \
    step2_SPAtests.R --help >/tmp/step2_help.txt && \
    grep -q tractorHybridPrefix /tmp/step2_help.txt && \
    rm -f /tmp/patch_saigetractor_hybrid.py /tmp/step2_help.txt

ENV PATH="/scripts/bin:${PATH}"
