CC=clang -O3 -mlzcnt

all: inverse25519_EEA_vartime.o

inverse25519_EEA_vartime.o: inverse25519_EEA_vartime.c
	$(CC) -c inverse25519_EEA_vartime.c

clean:
	-rm -f inverse25519_EEA_vartime.o
