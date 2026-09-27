# syntax=docker/dockerfile:1
# FELIXla image with the complete FELIXla binary and VCF chunk planner.
#
# htslib is built here rather than taken from the distribution, because the
# packaged libhts is not guaranteed to be built against libdeflate, which is
# most of the cost of a conversion. The executables are then linked statically,
# so the runtime stage carries no library dependency at all.

FROM debian:bookworm-slim AS build
ARG DEBIAN_FRONTEND=noninteractive
ARG FELIXLA_CXXFLAGS="-O3 -std=c++17 -Wall -Wextra -Wpedantic"
# Kept in step with HTSLIB_VERSION in .github/workflows/static-release.yml,
# so the image and the released binaries carry the same htslib.
ARG HTSLIB_VERSION=1.23.1

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bash \
        build-essential \
        bzip2 \
        ca-certificates \
        curl \
        libbz2-dev \
        libdeflate-dev \
        liblzma-dev \
        make \
        pkg-config \
        python3 \
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

RUN set -eux; \
    curl -fL \
      "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" \
      -o /tmp/htslib.tar.bz2; \
    mkdir -p /tmp/htslib-src; \
    tar -xjf /tmp/htslib.tar.bz2 -C /tmp/htslib-src --strip-components=1; \
    cd /tmp/htslib-src; \
    ./configure \
      --prefix=/opt/htslib \
      --with-libdeflate \
      --disable-libcurl \
      --disable-s3 \
      --disable-gcs \
      --disable-plugins; \
    make -j"$(nproc)"; \
    make install; \
    rm -rf /tmp/htslib.tar.bz2 /tmp/htslib-src

WORKDIR /src

COPY Makefile README.md ./
COPY src ./src
COPY tests ./tests
COPY testdata ./testdata
COPY resources ./resources

ENV PKG_CONFIG_PATH=/opt/htslib/lib/pkgconfig

# The scalar-path equivalence suite runs a reduced case count here: this image
# is also built for arm64 under emulation, where each case costs many times
# what it does natively. The release workflow runs the full count.
RUN set -eux; \
    make static STATIC_FULLY=1 CXXFLAGS="${FELIXLA_CXXFLAGS}"; \
    make test-static STATIC_FULLY=1 CXXFLAGS="${FELIXLA_CXXFLAGS}"; \
    python3 tests/run_keep_extract_intense.py --bin-dir "/src/bin-static"; \
    python3 tests/run_concat_intense.py --bin-dir "/src/bin-static"; \
    python3 tests/run_scalar_equivalence.py --bin-dir "/src/bin-static" --cases 12; \
    ldd bin-static/felixla 2>&1 | grep -q "not a dynamic executable"; \
    ldd bin-static/vcf_tbi_chunks 2>&1 | grep -q "not a dynamic executable"; \
    strip bin-static/felixla bin-static/vcf_tbi_chunks

FROM debian:bookworm-slim
ARG DEBIAN_FRONTEND=noninteractive

LABEL org.opencontainers.image.title="FELIXla"
LABEL org.opencontainers.image.description="Single-binary FELIX local-ancestry storage command for scalable ancestry-aware mixed-model analyses"
LABEL org.opencontainers.image.source="https://github.com/yorkklause/FELIXla"
LABEL org.opencontainers.image.licenses="GPL-3.0"

COPY --from=build /src/bin-static/felixla /usr/local/bin/felixla
COPY --from=build /src/bin-static/vcf_tbi_chunks /usr/local/bin/vcf_tbi_chunks
COPY --from=build /src/README.md /opt/felixla/README.md
COPY --from=build /src/resources/shapeit5_chunks /opt/felixla/resources/shapeit5_chunks

ENV FELIXLA_HOME=/opt/felixla
ENV FELIXLA_RESOURCES=/opt/felixla/resources

WORKDIR /data

ENTRYPOINT ["felixla"]
CMD ["--help"]
