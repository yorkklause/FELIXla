# designed by Kai, implemented by codex

CXX ?= g++
PKG_CONFIG ?= pkg-config

BIN_DIR ?= bin
STATIC_BIN_DIR ?= bin-static
PACK_TARGET := $(BIN_DIR)/flare_subset_to_tractor_hybrid
VCF_TARGET := $(BIN_DIR)/tractor_hybrid_to_vcf
CMP_TARGET := $(BIN_DIR)/compare_vcfs
DOSAGE_TARGET := $(BIN_DIR)/tractor_dosage_vcf_to_hybrid
RFMIX_MSP_TARGET := $(BIN_DIR)/rfmix_msp_to_tractor_hybrid
EXTRACT_TARGET := $(BIN_DIR)/tractor_hybrid_extract_region
ADMIX_TARGET := $(BIN_DIR)/calc_tractor_admixture
FELIX_QUERY_TARGET := $(BIN_DIR)/felixla_query
FELIX_CLI_TARGET := $(BIN_DIR)/felixla
SHAPEIT_ARGS_TARGET := $(BIN_DIR)/shapeit_chunks_to_tractor_args
STATIC_PACK_TARGET := $(STATIC_BIN_DIR)/flare_subset_to_tractor_hybrid
STATIC_VCF_TARGET := $(STATIC_BIN_DIR)/tractor_hybrid_to_vcf
STATIC_CMP_TARGET := $(STATIC_BIN_DIR)/compare_vcfs
STATIC_DOSAGE_TARGET := $(STATIC_BIN_DIR)/tractor_dosage_vcf_to_hybrid
STATIC_RFMIX_MSP_TARGET := $(STATIC_BIN_DIR)/rfmix_msp_to_tractor_hybrid
STATIC_EXTRACT_TARGET := $(STATIC_BIN_DIR)/tractor_hybrid_extract_region
STATIC_ADMIX_TARGET := $(STATIC_BIN_DIR)/calc_tractor_admixture
STATIC_FELIX_QUERY_TARGET := $(STATIC_BIN_DIR)/felixla_query
STATIC_FELIX_CLI_TARGET := $(STATIC_BIN_DIR)/felixla
STATIC_SHAPEIT_ARGS_TARGET := $(STATIC_BIN_DIR)/shapeit_chunks_to_tractor_args
PACK_SRC := src/flare_subset_to_tractor_hybrid.cpp
VCF_SRC := src/tractor_hybrid_to_vcf.cpp
CMP_SRC := src/compare_vcfs.cpp
DOSAGE_SRC := src/tractor_dosage_vcf_to_hybrid.cpp
RFMIX_MSP_SRC := src/rfmix_msp_to_tractor_hybrid.cpp
EXTRACT_SRC := src/tractor_hybrid_extract_region.cpp
ADMIX_SRC := src/calc_tractor_admixture.cpp
FELIX_QUERY_SRC := src/felixla_query.cpp
FELIX_CLI_SRC := scripts/felixla
SHAPEIT_ARGS_SRC := scripts/shapeit_chunks_to_tractor_args.sh

PKG_HTSLIB_CFLAGS := $(shell $(PKG_CONFIG) --cflags htslib 2>/dev/null)
PKG_HTSLIB_LIBS := $(shell $(PKG_CONFIG) --libs htslib 2>/dev/null)
PKG_HTSLIB_STATIC_LIBS_RAW := $(shell $(PKG_CONFIG) --libs --static htslib 2>/dev/null)
PKG_HTSLIB_LIBDIR := $(shell $(PKG_CONFIG) --variable=libdir htslib 2>/dev/null)
PKG_HTSLIB_STATIC_ARCHIVE := $(PKG_HTSLIB_LIBDIR)/libhts.a
PKG_LIBDEFLATE_STATIC_LIBS := $(shell $(PKG_CONFIG) --libs --static libdeflate 2>/dev/null)
PKG_LIBDEFLATE_LIBDIR := $(shell $(PKG_CONFIG) --variable=libdir libdeflate 2>/dev/null)
PKG_LIBDEFLATE_STATIC_ARCHIVE := $(PKG_LIBDEFLATE_LIBDIR)/libdeflate.a
ifneq ($(wildcard $(PKG_LIBDEFLATE_STATIC_ARCHIVE)),)
PKG_LIBDEFLATE_STATIC_LIBS := $(PKG_LIBDEFLATE_STATIC_ARCHIVE)
endif
PKG_LIBLZMA_STATIC_LIBS := $(shell $(PKG_CONFIG) --libs --static liblzma 2>/dev/null)
PKG_LIBLZMA_LIBDIR := $(shell $(PKG_CONFIG) --variable=libdir liblzma 2>/dev/null)
PKG_LIBLZMA_STATIC_ARCHIVE := $(PKG_LIBLZMA_LIBDIR)/liblzma.a
ifneq ($(wildcard $(PKG_LIBLZMA_STATIC_ARCHIVE)),)
PKG_LIBLZMA_STATIC_LIBS := $(PKG_LIBLZMA_STATIC_ARCHIVE) $(filter-out -llzma -pthread -lpthread,$(PKG_LIBLZMA_STATIC_LIBS))
endif
PKG_LIBCURL_STATIC_LIBS := $(filter-out -lz,$(shell $(PKG_CONFIG) --libs --static libcurl 2>/dev/null))
HTSLIB_STATIC_DEPS := $(filter-out -lhts,$(PKG_HTSLIB_STATIC_LIBS_RAW))
ifneq ($(strip $(PKG_LIBDEFLATE_STATIC_LIBS)),)
HTSLIB_STATIC_DEPS := $(filter-out -ldeflate,$(HTSLIB_STATIC_DEPS)) $(PKG_LIBDEFLATE_STATIC_LIBS)
endif
ifneq ($(strip $(PKG_LIBLZMA_STATIC_LIBS)),)
HTSLIB_STATIC_DEPS := $(filter-out -llzma,$(HTSLIB_STATIC_DEPS)) $(PKG_LIBLZMA_STATIC_LIBS)
endif
ifneq ($(strip $(PKG_LIBCURL_STATIC_LIBS)),)
HTSLIB_STATIC_DEPS := $(filter-out -lcurl,$(HTSLIB_STATIC_DEPS)) $(PKG_LIBCURL_STATIC_LIBS)
endif
PKG_HTSLIB_STATIC_LIBS :=
ifneq ($(wildcard $(PKG_HTSLIB_STATIC_ARCHIVE)),)
PKG_HTSLIB_STATIC_LIBS := $(PKG_HTSLIB_STATIC_ARCHIVE) $(HTSLIB_STATIC_DEPS)
else
PKG_HTSLIB_STATIC_LIBS := $(PKG_HTSLIB_STATIC_LIBS_RAW)
endif
BREW_HTSLIB_PREFIX := $(shell brew --prefix htslib 2>/dev/null)
BREW_HTSLIB_CFLAGS :=
BREW_HTSLIB_LIBS :=
BREW_HTSLIB_STATIC_LIBS :=
ifneq ($(wildcard $(BREW_HTSLIB_PREFIX)/include/htslib/hts.h),)
BREW_HTSLIB_CFLAGS := -I$(BREW_HTSLIB_PREFIX)/include
BREW_HTSLIB_LIBS := -L$(BREW_HTSLIB_PREFIX)/lib -lhts
endif
ifneq ($(wildcard $(BREW_HTSLIB_PREFIX)/lib/libhts.a),)
BREW_HTSLIB_STATIC_LIBS := $(BREW_HTSLIB_PREFIX)/lib/libhts.a $(HTSLIB_STATIC_DEPS)
endif

HTSLIB_CFLAGS ?= $(if $(PKG_HTSLIB_CFLAGS),$(PKG_HTSLIB_CFLAGS),$(BREW_HTSLIB_CFLAGS))
HTSLIB_LIBS ?= $(if $(PKG_HTSLIB_LIBS),$(PKG_HTSLIB_LIBS),$(BREW_HTSLIB_LIBS))
ifeq ($(strip $(HTSLIB_LIBS)),)
HTSLIB_LIBS := -lhts
endif
HTSLIB_STATIC_LIBS ?= $(if $(PKG_HTSLIB_STATIC_LIBS),$(PKG_HTSLIB_STATIC_LIBS),$(BREW_HTSLIB_STATIC_LIBS))

CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic
STATIC_FULLY ?= 0
STATIC_LDFLAGS ?=
ifeq ($(STATIC_FULLY),1)
STATIC_LDFLAGS += -static
endif
CPPFLAGS += $(HTSLIB_CFLAGS)
LDLIBS += $(HTSLIB_LIBS)

.PHONY: all static clean check-deps check-static-deps test test-static test-intense

all: check-deps $(PACK_TARGET) $(VCF_TARGET) $(CMP_TARGET) $(DOSAGE_TARGET) $(RFMIX_MSP_TARGET) $(EXTRACT_TARGET) $(ADMIX_TARGET) $(FELIX_QUERY_TARGET) $(FELIX_CLI_TARGET) $(SHAPEIT_ARGS_TARGET)

static: check-static-deps $(STATIC_PACK_TARGET) $(STATIC_VCF_TARGET) $(STATIC_CMP_TARGET) $(STATIC_DOSAGE_TARGET) $(STATIC_RFMIX_MSP_TARGET) $(STATIC_EXTRACT_TARGET) $(STATIC_ADMIX_TARGET) $(STATIC_FELIX_QUERY_TARGET) $(STATIC_FELIX_CLI_TARGET) $(STATIC_SHAPEIT_ARGS_TARGET)

check-deps:
	@printf '#include <htslib/hts.h>\n#include <htslib/vcf.h>\n' | \
		$(CXX) $(CPPFLAGS) -x c++ -E - >/dev/null 2>&1 || \
		( echo "ERROR: htslib headers not found. Install htslib or set HTSLIB_CFLAGS/HTSLIB_LIBS." >&2; exit 1 )

check-static-deps: check-deps
	@if [ -z "$(strip $(HTSLIB_STATIC_LIBS))" ]; then \
		echo "ERROR: htslib static library not found. Install libhts.a or set HTSLIB_STATIC_LIBS." >&2; \
		exit 1; \
	fi

$(PACK_TARGET): $(PACK_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(VCF_TARGET): $(VCF_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(CMP_TARGET): $(CMP_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(DOSAGE_TARGET): $(DOSAGE_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(RFMIX_MSP_TARGET): $(RFMIX_MSP_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(EXTRACT_TARGET): $(EXTRACT_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(ADMIX_TARGET): $(ADMIX_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(FELIX_QUERY_TARGET): $(FELIX_QUERY_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

$(FELIX_CLI_TARGET): $(FELIX_CLI_SRC)
	mkdir -p $(BIN_DIR)
	install -m 755 $< $@

$(SHAPEIT_ARGS_TARGET): $(SHAPEIT_ARGS_SRC)
	mkdir -p $(BIN_DIR)
	install -m 755 $< $@

$(STATIC_PACK_TARGET): $(PACK_SRC)
	mkdir -p $(STATIC_BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(STATIC_LDFLAGS) $(HTSLIB_STATIC_LIBS)

$(STATIC_VCF_TARGET): $(VCF_SRC)
	mkdir -p $(STATIC_BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(STATIC_LDFLAGS) $(HTSLIB_STATIC_LIBS)

$(STATIC_CMP_TARGET): $(CMP_SRC)
	mkdir -p $(STATIC_BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(STATIC_LDFLAGS) $(HTSLIB_STATIC_LIBS)

$(STATIC_DOSAGE_TARGET): $(DOSAGE_SRC)
	mkdir -p $(STATIC_BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(STATIC_LDFLAGS) $(HTSLIB_STATIC_LIBS)

$(STATIC_RFMIX_MSP_TARGET): $(RFMIX_MSP_SRC)
	mkdir -p $(STATIC_BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(STATIC_LDFLAGS) $(HTSLIB_STATIC_LIBS)

$(STATIC_EXTRACT_TARGET): $(EXTRACT_SRC)
	mkdir -p $(STATIC_BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(STATIC_LDFLAGS)

$(STATIC_ADMIX_TARGET): $(ADMIX_SRC)
	mkdir -p $(STATIC_BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(STATIC_LDFLAGS) $(HTSLIB_STATIC_LIBS)

$(STATIC_FELIX_QUERY_TARGET): $(FELIX_QUERY_SRC)
	mkdir -p $(STATIC_BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(STATIC_LDFLAGS)

$(STATIC_FELIX_CLI_TARGET): $(FELIX_CLI_SRC)
	mkdir -p $(STATIC_BIN_DIR)
	install -m 755 $< $@

$(STATIC_SHAPEIT_ARGS_TARGET): $(SHAPEIT_ARGS_SRC)
	mkdir -p $(STATIC_BIN_DIR)
	install -m 755 $< $@

test: all
	BIN_DIR="$(abspath $(BIN_DIR))" bash tests/run_tiny.sh

test-intense: all
	python3 tests/run_keep_extract_intense.py --bin-dir "$(abspath $(BIN_DIR))"

test-static: static
	BIN_DIR="$(abspath $(STATIC_BIN_DIR))" bash tests/run_tiny.sh

clean:
	rm -rf $(BIN_DIR) $(STATIC_BIN_DIR)
