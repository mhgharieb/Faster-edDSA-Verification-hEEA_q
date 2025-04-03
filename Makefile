CC = clang
CC_ED25519 = gcc

CFLAGS_ed25519 = -O3 -m64 -DED25519_TEST
CFLAGS_half_size = -O3 -mlzcnt
CFLAGS_test = -O3

LD = clang
LDLIBS = -lgmp -lssl -lcrypto

HALFSIZE = src/half_size
CURVE448 = $(HALFSIZE)/curve448
CURVE25519 = $(HALFSIZE)/curve25519
ED25519 = src/ed25519-donna

INVERSE25519 = src/inverse25519
BINGCD=$(INVERSE25519)/bingcd/src
EEA_q=$(INVERSE25519)/EEA_q
SAFEGCD=$(INVERSE25519)/inverse25519skylake-20210110

OBJCURVE448 = $(CURVE448)/curve448_hEEA_vartime.o \
	$(CURVE448)/curve448_hEEA_div_vartime.o \
	$(CURVE448)/curve448_reduce_basis_vartime.o \
	$(CURVE448)/curve448_hgcd_vartime.o

OBJCURVE25519 = $(CURVE25519)/curve25519_hEEA_vartime.o \
	$(CURVE25519)/curve25519_hEEA_div_vartime.o \
	$(CURVE25519)/curve25519_reduce_basis_vartime.o \
	$(CURVE25519)/curve25519_hgcd_vartime.o

OBJVERIFICATION = $(CURVE25519)/curve25519_hEEA_vartime.o \
	$(CURVE25519)/curve25519_hgcd_vartime.o \
	$(ED25519)/ed25519.o

OBJINVERSE25519 = $(EEA_q)/inverse25519_EEA_vartime.o \
	$(BINGCD)/gf25519.o \
	$(SAFEGCD)/asm.o \
	$(SAFEGCD)/table.o


all: halfSize ed25519 inverse25519 testHalfSizeEd25519 testHalfSizeEd448 testSingle testBatch testInverse25519

halfSize:
	$(MAKE) -C $(HALFSIZE)

ed25519: $(ED25519)/ed25519.c
	$(CC_ED25519) $(CFLAGS_ed25519) -c -o $(ED25519)/ed25519.o $(ED25519)/ed25519.c

inverse25519:
	$(MAKE) -C $(INVERSE25519)

testHalfSizeEd448: $(OBJCURVE448) test/test_halfSize_ed448.c
	$(CC) $(CFLAGS_test) -o test_halfSize_ed448 test/test_halfSize_ed448.c $(OBJCURVE448) $(LDLIBS)

testHalfSizeEd25519: $(OBJCURVE25519) test/test_halfSize_ed25519.c
	$(CC) $(CFLAGS_test) -o test_halfSize_ed25519 test/test_halfSize_ed25519.c $(OBJCURVE25519) $(LDLIBS)

testSingle: $(OBJVERIFICATION) test/test_singleVerification.c
	$(CC) $(CFLAGS_test) -o test_singleVerification  test/test_singleVerification.c $(OBJVERIFICATION) $(LDLIBS)

testBatch: $(OBJVERIFICATION) test/test_batchVerification.c
	$(CC) $(CFLAGS_test) -o test_batchVerification  test/test_batchVerification.c $(OBJVERIFICATION) $(LDLIBS)

testInverse25519: $(OBJINVERSE25519) test/test_inverse25519.c
	$(CC) $(CFLAGS_test) -o test_inverse25519 test/test_inverse25519.c $(OBJINVERSE25519) $(LDLIBS)



clean:
	$(MAKE) -C $(HALFSIZE) clean && \
	$(MAKE) -C $(INVERSE25519) clean && \
	rm -f test_halfSize_ed448 test_halfSize_ed25519 test_singleVerification test_batchVerification test_inverse25519 $(OBJVERIFICATION)