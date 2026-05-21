# Dockerfile
# designed by Kai, implemented by codex

FROM alpine:3.20

LABEL org.opencontainers.image.title="tractor-hybrid-tools"
LABEL org.opencontainers.image.description="Standalone tractor_hybrid conversion and validation utilities"
LABEL org.opencontainers.image.source="https://github.com/yorkklause/tractor-hybrid-tools"

RUN apk add --no-cache bash

COPY prebuilt/linux-x86_64-static/compare_vcfs /usr/local/bin/compare_vcfs
COPY prebuilt/linux-x86_64-static/estimate_mac_threshold /usr/local/bin/estimate_mac_threshold
COPY prebuilt/linux-x86_64-static/flare_subset_to_tractor_hybrid /usr/local/bin/flare_subset_to_tractor_hybrid
COPY prebuilt/linux-x86_64-static/rfmix_msp_to_tractor_hybrid /usr/local/bin/rfmix_msp_to_tractor_hybrid
COPY prebuilt/linux-x86_64-static/tractor_dosage_vcf_to_hybrid /usr/local/bin/tractor_dosage_vcf_to_hybrid
COPY prebuilt/linux-x86_64-static/tractor_hybrid_to_vcf /usr/local/bin/tractor_hybrid_to_vcf
COPY prebuilt/linux-x86_64-static/README.md /opt/tractor-hybrid-tools/README.md
COPY prebuilt/linux-x86_64-static/SHA256SUMS /opt/tractor-hybrid-tools/SHA256SUMS

WORKDIR /data

CMD ["rfmix_msp_to_tractor_hybrid"]
