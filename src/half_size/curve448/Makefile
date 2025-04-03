CC = clang

CFLAGS_half_size = -O3 -mlzcnt

OBJSIG = curve448_hEEA_vartime.o curve448_hEEA_div_vartime.o curve448_reduce_basis_vartime.o curve448_hgcd_vartime.o

all: $(OBJSIG)

curve448_hEEA_vartime.o: curve448_hEEA_vartime.c
	$(CC) $(CFLAGS_half_size) -c -o curve448_hEEA_vartime.o curve448_hEEA_vartime.c
curve448_hEEA_div_vartime.o: curve448_hEEA_div_vartime.c
	$(CC) $(CFLAGS_half_size) -c -o curve448_hEEA_div_vartime.o curve448_hEEA_div_vartime.c -fPIC

curve448_reduce_basis_vartime.o: curve448_reduce_basis_vartime.c
	$(CC) $(CFLAGS_half_size) -c -o curve448_reduce_basis_vartime.o curve448_reduce_basis_vartime.c

curve448_hgcd_vartime.o: curve448_hgcd_vartime.c gmp-impl_custom.h
	$(CC) $(CFLAGS_half_size) -c -o curve448_hgcd_vartime.o curve448_hgcd_vartime.c

clean:
	-rm -f $(OBJSIG)