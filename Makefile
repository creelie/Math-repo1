# Build the sixteen-runner (K = 15 speeds) pipeline and the independent checks.
#   make            build everything into bin/
#   make selftest   run the cascade self-test and the shift-lemma checks
CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17
# the generator is built as in Allikvere's run scripts
ENGFLAGS := -O3 -march=native -std=c++20 -DBGK_INCREMENTAL_GAIN \
            -DBGK_DEFERRED_PROBES -DBGK_THRESHOLD_WITNESS
K        := 15
BIN      := bin
ENGINES  := $(BIN)/bgk15_nw2 $(BIN)/bgk15_nw3 $(BIN)/bgk15_nw4 $(BIN)/bgk15_nw5

all: $(ENGINES) $(BIN)/cascade_k15 $(BIN)/lift_check $(BIN)/lift_check_shift $(BIN)/lift_check_sym

$(BIN):
	mkdir -p $(BIN)

# Level-one cover generator (Allikvere, unchanged); NW 64-bit words hold the
# NC = (p-1)/2 classes, so NW = 2 serves p < 258, NW = 3 p < 386, and so on.
$(BIN)/bgk15_nw%: src/basegen_k_campaign.cpp src/bgk_incremental_gain.h | $(BIN)
	$(CXX) $(ENGFLAGS) -DK=$(K) -DNW=$* -DMAXN=$$((64*$*)) -o $@ $<

$(BIN)/cascade_k15: src/cascade_filter_k.cpp | $(BIN)
	$(CXX) $(CXXFLAGS) -DK=$(K) -o $@ $<

$(BIN)/lift_check: check/lift_check.cpp | $(BIN)
	$(CXX) $(CXXFLAGS) -DK=$(K) -o $@ $<

$(BIN)/lift_check_shift: check/lift_check_shift.cpp | $(BIN)
	$(CXX) $(CXXFLAGS) -DK=$(K) -o $@ $<

$(BIN)/lift_check_sym: check/lift_check_sym.cpp | $(BIN)
	$(CXX) $(CXXFLAGS) -DK=$(K) -o $@ $<

selftest: all
	$(BIN)/cascade_k15 --selftest
	python3 check/shift_lemma_check.py
	python3 check/bound_n15.py
	python3 check/shape_n15.py

clean:
	rm -rf $(BIN)

.PHONY: all selftest clean
