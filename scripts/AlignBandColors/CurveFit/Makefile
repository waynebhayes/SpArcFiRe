TARGET = runFit
LIBS = -lm
CC = g++ -std=c++11
CFLAGS = -g -Wall

.PHONY: clean all default

default: $(TARGET)
	chmod a+rwx $(TARGET)
all: $(TARGET)
	chmod a+x $(TARGET)

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.hpp)

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)
