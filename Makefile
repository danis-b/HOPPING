CC=g++
CFLAGS=-std=c++11 -c -Wall -O2 -I/opt/local/include  
LDFLAGS=
SOURCES=HOPPING.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE= HOPPING.x

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(LDFLAGS) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o HOPPING.x