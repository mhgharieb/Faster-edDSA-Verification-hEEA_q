#include <stdint.h>
#include <immintrin.h>
/*
 * Read the cycle counter. The 'lfence' call should guarantee enough
 * serialization without adding too much overhead (contrary to what,
 * say, 'cpuid' would do).
 */
static inline uint64_t
get_ticks(void)
{
	/*
	 * GCC seems not to have the __rdtsc() intrinsic.
	 */
#if defined __GNUC__ && !defined __clang__
	uint32_t hi, lo;

	_mm_lfence();
	__asm__ __volatile__ ("rdtsc" : "=d" (hi), "=a" (lo) : : );
	return ((uint64_t)hi << 32) | (uint64_t)lo;
#else
	_mm_lfence();
	return __rdtsc();
#endif
}