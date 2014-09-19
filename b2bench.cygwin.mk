CC = clang++
CFLAGS = -O2 -IBox2D_v2.2.1 -DNDEBUG=1 -I/usr/lib/gcc/i686-pc-cygwin/4.8.3/include/c++/ -I/usr/lib/gcc/i686-pc-cygwin/4.8.3/include/c++/x86_64-pc-cygwin/ -I/usr/lib/gcc/i686-pc-cygwin/4.8.3/include/c++/backward/ -I/usr/lib/gcc/i686-pc-cygwin/4.8.3/include/c++/i686-pc-cygwin/
LFLAGS = -lstdc++ -lm

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