CC=g++
CFLAGS=-std=c++14 -c -Wall -O2 -I /Users/danis/Science/local/xtl/include -I /Users/danis/Science/local/xtensor/include 
LDFLAGS=
SOURCES=test.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE= test

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(LDFLAGS) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o test
