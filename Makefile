CC = gcc
CFLAGS = -I/usr/local/include -m64 -g -Wall
LDFLAGS = -L/usr/local/lib 
LDLIBS = -lm

djrandom: rdrand.o djrandom.o djrandommodel.o aes128k128d.o aes128k128d.h djrandommodel.h
	$(CC) $(CFLAGS) $(LDFLAGS) rdrand.o djrandom.o djrandommodel.o aes128k128d.o  -o djrandom $(LDLIBS)

rdrand.o: rdrand.c rdrand.h
	$(CC) -c $(CFLAGS) -o rdrand.o rdrand.c

djrandom.o: djrandom.c aes128k128d.h djrandommodel.h
	$(CC) -c $(CFLAGS) -o djrandom.o djrandom.c

smoothmodel.o: djrandommodel.c aes128k128d.h
	$(CC) -c $(CFLAGS) -o djrandommodel.o djrandommodel.c

aes128k128d.o: aes128k128d.c aes128k128d.h
	$(CC) -c $(CFLAGS) -o aes128k128d.o aes128k128d.c

install:
	cp djrandom /usr/local/bin

clean:
	rm rdrand.o
	rm aes128k128d.o
	rm djrandom.o
	rm djrandommodel.o
	rm djrandom

