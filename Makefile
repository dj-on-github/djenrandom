CC = gcc
CFLAGS = -I/usr/local/include -m64 -g -Wall
LDFLAGS = -L/usr/local/lib 
LDLIBS = -lm

djenrandom: rdrand.o djenrandom.o djenrandommodel.o aes128k128d.o aes128k128d.h djenrandommodel.h
	$(CC) $(CFLAGS) $(LDFLAGS) rdrand.o djenrandom.o djenrandommodel.o aes128k128d.o  -o djenrandom $(LDLIBS)

rdrand.o: rdrand.c rdrand.h
	$(CC) -c $(CFLAGS) -o rdrand.o rdrand.c

djenrandom.o: djenrandom.c aes128k128d.h djenrandommodel.h
	$(CC) -c $(CFLAGS) -o djenrandom.o djenrandom.c

smoothmodel.o: djenrandommodel.c aes128k128d.h
	$(CC) -c $(CFLAGS) -o djenrandommodel.o djenrandommodel.c

aes128k128d.o: aes128k128d.c aes128k128d.h
	$(CC) -c $(CFLAGS) -o aes128k128d.o aes128k128d.c

aesni128k128d.o: aes128k128d.c aes128k128d.h
	$(CC) -c $(CFLAGS) -o aes128k128d.o aes128k128d.c
    
install:
	cp djenrandom /usr/local/bin

clean:
	rm rdrand.o
	rm aes128k128d.o
	rm djenrandom.o
	rm djenrandommodel.o
	rm djenrandom

