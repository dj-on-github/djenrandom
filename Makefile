CC = gcc
CFLAGS = -I/usr/local/include -Irdrand_stdint -m64 -Wall -maes
LDFLAGS = -L/usr/local/lib 
LDLIBS = -lm

djenrandom: rdrand_stdint.o markov2p.o djenrandom.o djenrandommodel.o aes128k128d.o cmac.o aes128k128d.h djenrandommodel.h cmac.h
	$(CC) $(CFLAGS) $(LDFLAGS) rdrand_stdint.o markov2p.o djenrandom.o djenrandommodel.o aes128k128d.o cmac.o -o djenrandom $(LDLIBS)

rdrand_stdint.o: rdrand_stdint/rdrand_stdint.c rdrand_stdint/rdrand_stdint.h
	$(CC) -c $(CFLAGS) -o rdrand_stdint.o rdrand_stdint/rdrand_stdint.c

markov2p.o: markov2p.c markov2p.h
	$(CC) -c $(CFLAGS) -o markov2p.o markov2p.c

cmac.o: cmac.c cmac.h aes128k128d.h
	$(CC) -c $(CFLAGS) -o cmac.o cmac.c

djenrandom.o: djenrandom.c cmac.h aes128k128d.h djenrandommodel.h rdrand.h markov2p.h
	$(CC) -c $(CFLAGS) -o djenrandom.o djenrandom.c

smoothmodel.o: djenrandommodel.c aes128k128d.h
	$(CC) -c $(CFLAGS) -o djenrandommodel.o djenrandommodel.c

aesni128k128d.o: aes128k128d.c aes128k128d.h
	$(CC) -c $(CFLAGS) -o aes128k128d.o aes128k128d.c
    
install:
	cp djenrandom /usr/local/bin

clean:
	rm rdrand_stdint.o
	rm markov2p.o
	rm aes128k128d.o
	rm djenrandom.o
	rm djenrandommodel.o
	rm djenrandom
	rm cmac.o


