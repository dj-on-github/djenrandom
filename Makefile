CC = gcc
CFLAGS = -I/usr/local/include -m64 -Wall -maes
LDFLAGS = -L/usr/local/lib 
LDLIBS = -lm

djenrandom: rdrand.o markov2p.o djenrandom.o djenrandommodel.o aes128k128d.o aes128k128d.h djenrandommodel.h
	$(CC) $(CFLAGS) $(LDFLAGS) rdrand.o markov2p.o djenrandom.o djenrandommodel.o aes128k128d.o  -o djenrandom $(LDLIBS)

rdrand.o: rdrand.c rdrand.h
	$(CC) -c $(CFLAGS) -o rdrand.o rdrand.c

markov2p.o: markov2p.c markov2p.h
	$(CC) -c $(CFLAGS) -o markov2p.o markov2p.c

djenrandom.o: djenrandom.c aes128k128d.h djenrandommodel.h rdrand.h markov2p.h
	$(CC) -c $(CFLAGS) -o djenrandom.o djenrandom.c

smoothmodel.o: djenrandommodel.c aes128k128d.h
	$(CC) -c $(CFLAGS) -o djenrandommodel.o djenrandommodel.c

aesni128k128d.o: aes128k128d.c aes128k128d.h
	$(CC) -c $(CFLAGS) -o aes128k128d.o aes128k128d.c
    
install:
	cp djenrandom /usr/local/bin

clean:
	rm rdrand.o
	rm markov2p.o
	rm aes128k128d.o
	rm djenrandom.o
	rm djenrandommodel.o
	rm djenrandom

