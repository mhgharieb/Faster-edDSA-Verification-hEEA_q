#include <stdint.h>

void curve25519_hgcd_vartime(
	uint64_t *, uint64_t *, const uint64_t *, int *, int *);

void curve25519_hgcd_vartime_enhance1(
	uint64_t *, uint64_t *, const uint64_t *, int *, int *);

void curve25519_hgcd_vartime_enhance2(
	uint64_t *, uint64_t *, const uint64_t *, int *, int *);

	