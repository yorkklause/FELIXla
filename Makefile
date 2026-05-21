# designed by Kai, implemented by codex

CXX ?= g++
PKG_CONFIG ?= pkg-config

BIN_DIR ?= bin
PACK_TARGET := $(BIN_DIR)/flare_subset_to_tractor_hybrid
VCF_TARGET := $(BIN_DIR)/tractor_hybrid_to_vcf
EST_TARGET := $(BIN_DIR)/estimate_mac_threshold
CMP_TARGET := $(BIN_DIR)/compare_vcfs
DOSAGE_TARGET := $(BIN_DIR)/tractor_dosage_vcf_to_hybrid
PACK_SRC := src/flare_subset_to_tractor_hybrid.cpp
VCF_SRC := src/tractor_hybrid_to_vcf.cpp
EST_SRC := src/estimate_mac_threshold.cpp
CMP_SRC := src/compare_vcfs.cpp
DOSAGE_SRC := src/tractor_dosage_vcf_to_hybrid.cpp

PKG_HTSLIB_CFLAGS := $(shell $(PKG_CONFIG) --cflags htslib 2>/dev/null)
PKG_HTSLIB_LIBS := $(shell $(PKG_CONFIG) --libs htslib 2>/dev/null)
BREW_HTSLIB_PREFIX := $(shell brew --prefix htslib 2>/dev/null)
BREW_HTSLIB_CFLAGS :=
BREW_HTSLIB_LIBS :=
ifneq ($(wildcard $(BREW_HTSLIB_PREFIX)/include/htslib/hts.h),)
BREW_HTSLIB_CFLAGS := -I$(BREW_HTSLIB_PREFIX)/include
BREW_HTSLIB_LIBS := -L$(BREW_HTSLIB_PREFIX)/lib -lhts
endif

HTSLIB_CFLAGS ?= $(if $(PKG_HTSLIB_CFLAGS),$(PKG_HTSLIB_CFLAGS),$(BREW_HTSLIB_CFLAGS))
HTSLIB_LIBS ?= $(if $(PKG_HTSLIB_LIBS),$(PKG_HTSLIB_LIBS),$(BREW_HTSLIB_LIBS))
ifeq ($(strip $(HTSLIB_LIBS)),)
HTSLIB_LIBS := -lhts
endif

CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS += $(HTSLIB_CFLAGS)
LDLIBS += $(HTSLIB_LIBS)

.PHONY: all clean check-deps test

all: check-deps $(PACK_TARGET) $(VCF_TARGET) $(EST_TARGET) $(CMP_TARGET) $(DOSAGE_TARGET)

check-deps:
	@printf '#include <htslib/hts.h>\n#include <htslib/vcf.h>\n' | \
		$(CXX) $(CPPFLAGS) -x c++ -E - >/dev/null 2>&1 || \
		( echo "ERROR: htslib headers not found. Install htslib or set HTSLIB_CFLAGS/HTSLIB_LIBS." >&2; exit 1 )

$(PACK_TARGET): $(PACK_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(VCF_TARGET): $(VCF_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(EST_TARGET): $(EST_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(CMP_TARGET): $(CMP_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(DOSAGE_TARGET): $(DOSAGE_SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

test: all
	bash tests/run_tiny.sh

clean:
	rm -rf $(BIN_DIR)
