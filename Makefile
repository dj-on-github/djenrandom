CC = gcc
CFLAGS = -I/usr/local/include -m64 -g -Wall
LDFLAGS = -L/usr/local/lib 
LDLIBS = -lm

genrandom: rdrand.o genrandom.o genrandommodel.o aes128k128d.o aes128k128d.h genrandommodel.h
	$(CC) $(CFLAGS) $(LDFLAGS) rdrand.o genrandom.o genrandommodel.o aes128k128d.o  -o genrandom $(LDLIBS)

rdrand.o: rdrand.c rdrand.h
	$(CC) -c $(CFLAGS) -o rdrand.o rdrand.c

genrandom.o: genrandom.c aes128k128d.h genrandommodel.h
	$(CC) -c $(CFLAGS) -o genrandom.o genrandom.c

smoothmodel.o: genrandommodel.c aes128k128d.h
	$(CC) -c $(CFLAGS) -o genrandommodel.o genrandommodel.c

aes128k128d.o: aes128k128d.c aes128k128d.h
	$(CC) -c $(CFLAGS) -o aes128k128d.o aes128k128d.c

install:
	cp genrandom /usr/local/bin

clean:
	rm rdrand.o
	rm aes128k128d.o
	rm genrandom.o
	rm genrandommodel.o
	rm genrandom

