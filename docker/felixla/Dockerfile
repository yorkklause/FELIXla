# syntax=docker/dockerfile:1
# FELIXla image with the complete FELIXla binary and VCF chunk planner.

FROM debian:bookworm-slim AS build
ARG DEBIAN_FRONTEND=noninteractive
ARG FELIXLA_CXXFLAGS="-O3 -std=c++17 -Wall -Wextra -Wpedantic"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bash \
        build-essential \
        ca-certificates \
        libhts-dev \
        make \
        pkg-config \
        python3 && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /src

COPY Makefile README.md ./
COPY src ./src
COPY tests ./tests
COPY testdata ./testdata
COPY resources ./resources

RUN make all CXXFLAGS="${FELIXLA_CXXFLAGS}" && \
    make test CXXFLAGS="${FELIXLA_CXXFLAGS}" && \
    make test-intense CXXFLAGS="${FELIXLA_CXXFLAGS}"

FROM debian:bookworm-slim
ARG DEBIAN_FRONTEND=noninteractive

LABEL org.opencontainers.image.title="FELIXla"
LABEL org.opencontainers.image.description="Single-binary FELIX local-ancestry storage command for scalable ancestry-aware mixed-model analyses"
LABEL org.opencontainers.image.source="https://github.com/yorkklause/FELIXla"
LABEL org.opencontainers.image.licenses="GPL-3.0"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        libhts3 && \
    rm -rf /var/lib/apt/lists/*

COPY --from=build /src/bin/felixla /usr/local/bin/felixla
COPY --from=build /src/bin/vcf_tbi_chunks /usr/local/bin/vcf_tbi_chunks
COPY --from=build /src/README.md /opt/felixla/README.md
COPY --from=build /src/resources/shapeit5_chunks /opt/felixla/resources/shapeit5_chunks

ENV FELIXLA_HOME=/opt/felixla
ENV FELIXLA_RESOURCES=/opt/felixla/resources

WORKDIR /data

ENTRYPOINT ["felixla"]
CMD ["--help"]
