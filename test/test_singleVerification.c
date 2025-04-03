#include <stdio.h>
#include <string.h>
#include "test-ticks.h"
#include "../src/ed25519-donna/ed25519.h"

#define number_of_samples 10000
#define number_of_rounds 20


struct benchmark_result{
	uint64_t best;
	uint64_t median;
	double average;
};

static int cmp_int64(const void *v1, const void *v2)
{
	int64_t x1, x2;

	x1 = *(const int64_t *)v1;
	x2 = *(const int64_t *)v2;
	if (x1 < x2) {
		return -1;
	} else if (x1 == x2) {
		return 0;
	} else {
		return 1;
	}
}


int test_instance(size_t test_count){
    ed25519_secret_key sks[test_count];
	ed25519_public_key pks[test_count];
	ed25519_signature sigs[test_count];
	unsigned char messages[test_count][128];
	size_t message_lengths[test_count];
	const unsigned char *message_pointers[test_count];
	const unsigned char *pk_pointers[test_count];
	const unsigned char *sig_pointers[test_count];
	int ret;
	size_t i, j;
	uint64_t t_begin;

	uint64_t t[test_count];
	uint64_t total_t = 0;
	struct benchmark_result benchmark_old; 
	struct benchmark_result benchmark_old_2B; 
	struct benchmark_result benchmark_hEEA; 
	struct benchmark_result benchmark_hEEA_samePre; 
	struct benchmark_result benchmark_gmp_hgcd; 

	/* generate keys */
	for (i = 0; i < test_count; i++) {
		ed25519_randombytes_unsafe(sks[i], sizeof(sks[i]));
		ed25519_publickey(sks[i], pks[i]);
		pk_pointers[i] = pks[i];
	}

	/* generate messages */
	ed25519_randombytes_unsafe(messages, sizeof(messages));
	for (i = 0; i < test_count; i++) {
		message_pointers[i] = messages[i];
		message_lengths[i] = (i & 127) + 1;
	}

	/* sign messages */
	for (i = 0; i < test_count; i++) {
		ed25519_sign(message_pointers[i], message_lengths[i], sks[i], pks[i], sigs[i]);
		sig_pointers[i] = sigs[i];
	}


	/* verify messages */	
			
	/* New single verification using QSM_B_B' with hgcd_enhance2 */
	total_t = 0;
	for (i = 0; i < test_count; i++) {
		t_begin = get_ticks();
		for(j=0;j<number_of_rounds;j++){
			ret = ed25519_sign_open_hgcd(message_pointers[i], message_lengths[i], pks[i], sigs[i]);
		}
		t[i] = get_ticks() - t_begin;
		total_t += t[i];
		if (ret){
			fprintf(stderr, "ERR: failed to open message %lu using new open using hgcd\n", i);
			exit(EXIT_FAILURE);
		}
	}
	qsort(t, test_count, sizeof(uint64_t), cmp_int64);
	benchmark_gmp_hgcd.best =  t[0]/number_of_rounds;
	benchmark_gmp_hgcd.median = t[test_count/2]/number_of_rounds;
	benchmark_gmp_hgcd.average =  total_t/(number_of_rounds*(double)test_count);
	
	/* Old single verification */
	total_t = 0;
	for (i = 0; i < test_count; i++) {
		t_begin = get_ticks();
		for(j=0;j<number_of_rounds;j++){
			ret = ed25519_sign_open(message_pointers[i], message_lengths[i], pks[i], sigs[i]);
		}
		t[i] = get_ticks() - t_begin;
		total_t += t[i];
		if (ret){
			fprintf(stderr, "ERR: failed to open message %lu\n", i);
			exit(EXIT_FAILURE);

		}
	}
	qsort(t, test_count, sizeof(uint64_t), cmp_int64);
	benchmark_old.best =  t[0]/number_of_rounds;
	benchmark_old.median = t[test_count/2]/number_of_rounds;
	benchmark_old.average =  total_t/(number_of_rounds*(double)test_count);

	/* Old single verification using DSM_B_doublePre in which the precomputed table for B is doubled */
	total_t = 0;
	for (i = 0; i < test_count; i++) {
		t_begin = get_ticks();
		for(j=0;j<number_of_rounds;j++){
			ret = ed25519_sign_open_2B(message_pointers[i], message_lengths[i], pks[i], sigs[i]);
		}
		t[i] = get_ticks() - t_begin;
		total_t += t[i];
		if (ret){
			fprintf(stderr, "ERR: failed to open message %lu\n", i);
			exit(EXIT_FAILURE);

		}
	}
	qsort(t, test_count, sizeof(uint64_t), cmp_int64);
	benchmark_old_2B.best =  t[0]/number_of_rounds;
	benchmark_old_2B.median = t[test_count/2]/number_of_rounds;
	benchmark_old_2B.average =  total_t/(number_of_rounds*(double)test_count);

	/* New single verification using QSM_B with hEEA */
	total_t = 0;
	for (i = 0; i < test_count; i++) {	
		t_begin = get_ticks();
		for(j=0;j<number_of_rounds;j++){
			ret = ed25519_sign_open_hEEA_samePre(message_pointers[i], message_lengths[i], pks[i], sigs[i]);
		}
		t[i] = get_ticks() - t_begin;
		total_t += t[i];
		if (ret){
			fprintf(stderr, "ERR: failed to open message %lu using new open using hEEA with a same Pre-table\n", i);
			exit(EXIT_FAILURE);
		}
	}
	qsort(t, test_count, sizeof(uint64_t), cmp_int64);
	benchmark_hEEA_samePre.best =  t[0]/number_of_rounds;
	benchmark_hEEA_samePre.median = t[test_count/2]/number_of_rounds;
	benchmark_hEEA_samePre.average =  total_t/(number_of_rounds*(double)test_count);

	/* New single verification using QSM_B_B' with hEEA */
	total_t = 0;
	for (i = 0; i < test_count; i++) {	
		t_begin = get_ticks();
		for(j=0;j<number_of_rounds;j++){
			ret = ed25519_sign_open_hEEA(message_pointers[i], message_lengths[i], pks[i], sigs[i]);
		}
		t[i] = get_ticks() - t_begin;
		total_t += t[i];
		if (ret){
			fprintf(stderr, "ERR: failed to open message %lu using new open using hEEA_q\n", i);
			exit(EXIT_FAILURE);
		}	
	}
	qsort(t, test_count, sizeof(uint64_t), cmp_int64);
	benchmark_hEEA.best =  t[0]/number_of_rounds;
	benchmark_hEEA.median = t[test_count/2]/number_of_rounds;
	benchmark_hEEA.average =  total_t/(number_of_rounds*(double)test_count);


	printf("───────────────────────────────────────────────────────────────────────────────\n");
	printf("Time (ticks)|       DSM_B       |  QSM_B (hEEA_q)  | Speed up     | Improvement\n");
	printf("───────────────────────────────────────────────────────────────────────────────\n");
	printf("Best        | %-17lu | %-16lu | %-12.4f | %.2f %%\n", benchmark_old.best, benchmark_hEEA_samePre.best, (double)benchmark_old.best/(double)benchmark_hEEA_samePre.best,(double)(benchmark_old.best - benchmark_hEEA_samePre.best)/((double)benchmark_old.best) * 100);
	printf("Median      | %-17lu | %-16lu | %-12.4f | %.2f %%\n", benchmark_old.median, benchmark_hEEA_samePre.median, (double)benchmark_old.median/(double)benchmark_hEEA_samePre.median,(double)(benchmark_old.median - benchmark_hEEA_samePre.median)/((double)benchmark_old.median) * 100);
	printf("Average     | %-17.2f | %-16.2f | %-12.4f | %.2f %%\n", benchmark_old.average, benchmark_hEEA_samePre.average, (double)benchmark_old.average/(double)benchmark_hEEA_samePre.average,(double)(benchmark_old.average - benchmark_hEEA_samePre.average)/((double)benchmark_old.average) * 100);
	printf("───────────────────────────────────────────────────────────────────────────────\n");

	printf("───────────────────────────────────────────────────────────────────────────────\n");
	printf("Time (ticks)|  DSM_B_doublePre  | QSM_B_B' (hEEA_q)| Speed up     | Improvement\n");
	printf("───────────────────────────────────────────────────────────────────────────────\n");
	printf("Best        | %-17lu | %-16lu | %-12.4f | %.2f %%\n", benchmark_old_2B.best, benchmark_hEEA.best, (double)benchmark_old_2B.best/(double)benchmark_hEEA.best,(double)(benchmark_old_2B.best - benchmark_hEEA.best)/((double)benchmark_old_2B.best) * 100);
	printf("Median      | %-17lu | %-16lu | %-12.4f | %.2f %%\n", benchmark_old_2B.median, benchmark_hEEA.median, (double)benchmark_old_2B.median/(double)benchmark_hEEA.median,(double)(benchmark_old_2B.median - benchmark_hEEA.median)/((double)benchmark_old_2B.median) * 100);
	printf("Average     | %-17.2f | %-16.2f | %-12.4f | %.2f %%\n", benchmark_old_2B.average, benchmark_hEEA.average, (double)benchmark_old_2B.average/(double)benchmark_hEEA.average,(double)(benchmark_old_2B.average - benchmark_hEEA.average)/((double)benchmark_old_2B.average) * 100);
	printf("───────────────────────────────────────────────────────────────────────────────\n");

	printf("───────────────────────────────────────────────────────────────────────────────\n");
	printf("Time (ticks)| QSM_B_B' (hEEA_q) |QSM_B_B' (hgcd_enhance2)|Speed up| Improvement\n");
	printf("───────────────────────────────────────────────────────────────────────────────\n");
	printf("Best        | %-17lu | %-22lu | %-6.4f | %.2f %%\n", benchmark_hEEA.best, benchmark_gmp_hgcd.best, (double)benchmark_hEEA.best/(double)benchmark_gmp_hgcd.best,((double)benchmark_hEEA.best - (double)benchmark_gmp_hgcd.best)/((double)benchmark_hEEA.best) * 100);
	printf("Median      | %-17lu | %-22lu | %-6.4f | %.2f %%\n", benchmark_hEEA.median, benchmark_gmp_hgcd.median, (double)benchmark_hEEA.median/(double)benchmark_gmp_hgcd.median,((double)benchmark_hEEA.median - (double)benchmark_gmp_hgcd.median)/((double)benchmark_hEEA.median) * 100);
	printf("Average     | %-17.2f | %-22.2f | %-6.4f | %.2f %%\n", benchmark_hEEA.average, benchmark_gmp_hgcd.average, (double)benchmark_hEEA.average/(double)benchmark_gmp_hgcd.average,(double)(benchmark_hEEA.average - benchmark_gmp_hgcd.average)/((double)benchmark_hEEA.average) * 100);
	printf("───────────────────────────────────────────────────────────────────────────────\n");


	return 0;

}

int main(){
	printf("Benchmark of individual verification:\n");
	printf("Number of samples = %i \n", number_of_samples);
	printf("Number of rounds = %i \n", number_of_rounds);
	test_instance(number_of_samples);
	printf("Done!\n");
		
}