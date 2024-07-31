all: analyzer univmake chi_square_cc0pi

univmake: univmake.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

analyzer: analyzer.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

chi_square_cc0pi: chi_square_cc0pi.cpp
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

.PHONY: clean

clean:
	$(RM) univmake analyzer chi_square_cc0pi
