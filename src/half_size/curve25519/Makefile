CC = clang

CFLAGS_half_size = -O3 -mlzcnt

OBJSIG = curve25519_hEEA_vartime.o curve25519_hEEA_div_vartime.o curve25519_reduce_basis_vartime.o curve25519_hgcd_vartime.o

all: $(OBJSIG)

curve25519_hEEA_vartime.o: curve25519_hEEA_vartime.c
	$(CC) $(CFLAGS_half_size) -c -o curve25519_hEEA_vartime.o curve25519_hEEA_vartime.c
curve25519_hEEA_div_vartime.o: curve25519_hEEA_div_vartime.c
	$(CC) $(CFLAGS_half_size) -c -o curve25519_hEEA_div_vartime.o curve25519_hEEA_div_vartime.c -fPIC

curve25519_reduce_basis_vartime.o: curve25519_reduce_basis_vartime.c
	$(CC) $(CFLAGS_half_size) -c -o curve25519_reduce_basis_vartime.o curve25519_reduce_basis_vartime.c

curve25519_hgcd_vartime.o: curve25519_hgcd_vartime.c gmp-impl_custom.h
	$(CC) $(CFLAGS_half_size) -c -o curve25519_hgcd_vartime.o curve25519_hgcd_vartime.c

clean:
	-rm -f $(OBJSIG)