all: analyzer analyzer_ru univmake

univmake: univmake.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

analyzer: analyzer.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

analyzer_ru: analyzer_ru.C
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

.PHONY: clean

clean:
	$(RM) univmake analyzer
