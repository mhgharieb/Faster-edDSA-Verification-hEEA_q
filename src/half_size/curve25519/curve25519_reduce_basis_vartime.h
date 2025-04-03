#include <stdint.h>
#include <string.h>

void curve25519_reduce_basis_vartime(
	uint64_t *c0, uint64_t *c1, const uint64_t *b);