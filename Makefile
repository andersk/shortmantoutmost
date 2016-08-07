LEMON_CFLAGS := $(shell pkg-config --cflags lemon)
LEMON_LDFLAGS := $(shell pkg-config --libs lemon)
CXXFLAGS := -O2 -Wall -std=c++11 $(LEMON_CFLAGS)
LDFLAGS := $(LEMON_LDFLAGS)

generate: generate.cc

clean:
	rm -f generate
