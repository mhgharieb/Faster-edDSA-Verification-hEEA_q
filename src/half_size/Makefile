CURVE448=curve448
CURVE25519=curve25519

.PHONY: all curve448 curve25519

all: curve448 curve25519

curve448:
	$(MAKE) -C $(CURVE448)

curve25519:
	$(MAKE) -C $(CURVE25519)

clean:
	$(MAKE) -C $(CURVE448) clean && \
	$(MAKE) -C $(CURVE25519)  clean