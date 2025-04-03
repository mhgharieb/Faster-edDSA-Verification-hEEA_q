#include <stdio.h>
#include <string.h>
#include "../src/ed25519-donna/ed25519.h"

#include "test-ticks.h"

#define max_batch_size 128
#define number_of_samples 100
#define number_of_rounds 10


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


int test_batch_instance(int test_count, int batch_size){
    ed25519_secret_key sks[batch_size];
	ed25519_public_key pks[batch_size];
	ed25519_signature sigs[batch_size];
	unsigned char messages[batch_size][128];
	size_t message_lengths[batch_size];
	const unsigned char *message_pointers[batch_size];
	const unsigned char *pk_pointers[batch_size];
	const unsigned char *sig_pointers[batch_size];
	int valid[batch_size], ret, validret;
	size_t i, j;
	uint64_t t_begin;
	double total_time_new=0,total_time_new_hgcd=0, total_time_old=0;

	for (j=0; j<test_count; j++){
		/* generate keys */
		for (i = 0; i < batch_size; i++) {
			ed25519_randombytes_unsafe(sks[i], sizeof(sks[i]));
			ed25519_publickey(sks[i], pks[i]);
			pk_pointers[i] = pks[i];
		}

		/* generate messages */
		ed25519_randombytes_unsafe(messages, sizeof(messages));
		for (i = 0; i < batch_size; i++) {
			message_pointers[i] = messages[i];
			message_lengths[i] = (i & 127) + 1;
		}

		/* sign messages */
		for (i = 0; i < batch_size; i++) {
			ed25519_sign(message_pointers[i], message_lengths[i], sks[i], pks[i], sigs[i]);
			sig_pointers[i] = sigs[i];
		}

    	/* batch verify */

		/* Old approach */
		t_begin = get_ticks();
		for(i=0;i<number_of_rounds;i++)
			ret = ed25519_sign_open_batch(message_pointers, message_lengths, pk_pointers, sig_pointers, batch_size, valid);
		total_time_old += get_ticks() - t_begin;
		if (ret){
			fprintf(stderr, "ERR: Old Batch verification failed\n");
			exit(EXIT_FAILURE);
		}

		/* New batch verification using the new randamization with hEEA */
		t_begin = get_ticks();
		for(i=0;i<number_of_rounds;i++)
			ret = ed25519_sign_open_batch_hEEA(message_pointers, message_lengths, pk_pointers, sig_pointers, batch_size, valid);
		total_time_new += get_ticks() - t_begin;
		if (ret){
			fprintf(stderr, "ERR: New Batch verification failed\n");
			exit(EXIT_FAILURE);
		}

		/* New batch verification using the new randamization with hgcd */
		t_begin = get_ticks();
		for(i=0;i<number_of_rounds;i++)
			ret = ed25519_sign_open_batch_hgcd(message_pointers, message_lengths, pk_pointers, sig_pointers, batch_size, valid);
		total_time_new_hgcd += get_ticks() - t_begin;
		if (ret){
			fprintf(stderr, "ERR: New Batch verification failed\n");
			exit(EXIT_FAILURE);
		}

	}

	total_time_old = total_time_old/((double)(number_of_rounds*test_count*batch_size));
	total_time_new = total_time_new/((double)(number_of_rounds*test_count*batch_size));
	total_time_new_hgcd = total_time_new_hgcd/((double)(number_of_rounds*test_count*batch_size));



	
	// printf("──────────────────────────────────────────────────────────────────────\n");
	// printf("Batch size   | Old approach | New approach | Speed up     | Impovement\n");
	// printf("──────────────────────────────────────────────────────────────────────\n");
	// printf("Best         | %-12lu | %-12lu | %-12.2f | %.2f %%\n", t_old[0], t_new[0], (double)t_old[0]/(double)t_new[0],(double)(t_old[0] - t_new[0])/((double)t_old[0]) * 100);
	// printf("Median       | %-12lu | %-12lu | %-12.2f | %.2f %%\n", t_old[test_count/2], t_new[test_count/2], (double)t_old[test_count/2]/(double)t_new[test_count/2], (double)(t_old[test_count/2] - t_new[test_count/2])/((double)t_old[test_count/2]) * 100);
	printf("%-10i | %-12.2f | %-14.2f | %-8.4f | %-5.2f %%     || %-18.2f | %-12.4f | %.2f %%\n", batch_size, total_time_old, total_time_new ,total_time_old/total_time_new ,(total_time_old - total_time_new)/total_time_old * 100, total_time_new_hgcd, total_time_new/total_time_new_hgcd, (total_time_new-total_time_new_hgcd)/total_time_new_hgcd * 100);
	// printf("──────────────────────────────────────────────────────────────────────\n");

	return 0;

}


int main(){
	printf("Benchmark of batch verification:\n");
	printf("Number of samples = %i \n", number_of_samples);
	printf("Number of rounds = %i \n", number_of_rounds);
	
	printf("───────────────────────────────────────Average Time (ticks/verification)──────────────────────────────────────────────\n");
	printf("──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n");
	printf("Batch size | Old approach | New using hEEA | Speed up | Improvement || New using GMP_hgcd | Speed up     | Improvement\n");
	printf("──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n");
	for(int batch_size=4;batch_size<=max_batch_size;batch_size++)
		test_batch_instance(number_of_samples, batch_size);
	
	printf("──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────\n");

	printf("Done!\n");	
}
