CC = $(LLVM)/clang++
CFLAGS = -O2 -I./ -DNDEBUG=1
LFLAGS = -lstdc++

OBJECTS = \
b2bench.o \
b2Math.o \
b2Settings.o

all: b2bench

%.o: %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

b2bench: $(OBJECTS)
	$(CC) $(LFLAGS) -o $@ $(OBJECTS)

clean:
	rm $(OBJECTS)
