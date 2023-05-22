CC = gcc -g -Wall
CFLAGS = `pkg-config --cflags glib-2.0`
LIBS = `pkg-config --libs glib-2.0` -lm

ITEST = hefsi.o moller-tri.o itest.o

all: itest

itest: $(ITEST) hefsi.h
	$(CC) $(CFLAGS) $(ITEST) $(LIBS) -o itest

.c.o:
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f *.c~ *.h~ *.o
