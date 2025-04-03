#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <openssl/rand.h>
#include "test-ticks.h"
#include "../src/inverse25519/EEA_q/inverse25519_EEA_vartime.h"
#include "../src/inverse25519/inverse25519skylake-20210110/inverse25519.h"
#include "../src/inverse25519/bingcd/src/gf25519.h"

#define number_of_samples 10000

struct benchmark_result{
	uint64_t best;
	uint64_t median;
	double average;
};

void gmp_import(mpz_t z,const unsigned char *s,unsigned long long slen)
{
  mpz_import(z,slen,-1,1,0,0,s);
}

int gmp_export(unsigned char *s,unsigned long long slen,mpz_t z)
{
  unsigned long long i;
  if (mpz_sizeinbase(z,256) > slen) return -1;
  for (i = 0;i < slen;++i) s[i] = 0;
  mpz_export(s,0,-1,1,0,0,z);
  return 0;
}



static void
gf_to_mpz(mpz_t r, const gf *a)
{
	mpz_set_ui(r, (uint32_t)(a->v3 >> 32));
	mpz_mul_2exp(r, r, 32);
	mpz_add_ui(r, r, (uint32_t)a->v3);

	mpz_mul_2exp(r, r, 32);
	mpz_add_ui(r, r, (uint32_t)(a->v2 >> 32));
	mpz_mul_2exp(r, r, 32);
	mpz_add_ui(r, r, (uint32_t)a->v2);

	mpz_mul_2exp(r, r, 32);
	mpz_add_ui(r, r, (uint32_t)(a->v1 >> 32));
	mpz_mul_2exp(r, r, 32);
	mpz_add_ui(r, r, (uint32_t)a->v1);

	mpz_mul_2exp(r, r, 32);
	mpz_add_ui(r, r, (uint32_t)(a->v0 >> 32));
	mpz_mul_2exp(r, r, 32);
	mpz_add_ui(r, r, (uint32_t)a->v0);
}

static void
mpz_2_gf(gf *out, const mpz_t in_mpz)
{
	out->v0 = in_mpz->_mp_d[0];
	out->v1 = in_mpz->_mp_d[1];
	out->v2 = in_mpz->_mp_d[2];
	out->v3 = in_mpz->_mp_d[3];
	
}




/* u64 to mpz */
void u64_2_mpz(mpz_t out, const uint64_t* in, size_t l) 
{
	mpz_import(out, l, -1, sizeof(uint64_t), 0, 0, in);
}

/* mpz to u64 */
void mpz_2_u64(uint64_t *out, const mpz_t in, size_t l) {
    size_t count,i;
    mpz_export(out, &count, -1, sizeof(uint64_t), 0, 0, in);
    if (count < l){
        for(i=count; i<l; i++){
            out[i] = 0;
        };
    };
};


void rand_mpz(mpz_t b_mpz, gmp_randstate_t state, const mpz_t L){
	for(;;){
		mpz_urandomb(b_mpz, state, 255);
		if(mpz_cmp_ui(b_mpz, 0))
			break;
	}
	
	mpz_mod(b_mpz, b_mpz, L);
}

/* check if 1 = x_inv * x mod p */
int check_correctness(const mpz_t x_inv_mpz, const mpz_t x_mpz, const mpz_t p_mpz){

	mpz_t tmp;
	mpz_init(tmp);
	mpz_mul(tmp, x_inv_mpz, x_mpz);
	mpz_mod(tmp, tmp, p_mpz);
	return mpz_cmp_ui(tmp, 1) ? 0:1;	

}

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

int test_instance(int test_count){

	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, time(NULL));

	mpz_t p_mpz;
	mpz_init(p_mpz);
	mpz_set_str(p_mpz, "57896044618658097711785492504343953926634992332820282019728792003956564819949", 10);


	mpz_t x_mpz[test_count], x_inv_mpz;
	unsigned char x_c[32], x_inv_c[32];
	gf x_gf, x_inv_gf;
	uint64_t x[4], x_inv[4];
	mpz_init(x_inv_mpz);

	uint64_t t_begin;

	uint64_t t[test_count];
	uint64_t total_t = 0;
	struct benchmark_result benchmark_safegcd; 
	struct benchmark_result benchmark_EEA; 
	struct benchmark_result benchmark_bingcd; 
	
	for(size_t j=0; j<test_count; j++)
	{
		/* pick a random b, s.t. 0 < x < p */
		mpz_init(x_mpz[j]);
		rand_mpz(x_mpz[j], state, p_mpz);
	}


	/* Test bingcd */
	total_t = 0;
	for(size_t j=0; j<test_count; j++)
	{
		
		mpz_2_gf(&x_gf, x_mpz[j]);
    	
		t_begin = get_ticks();

		for(int i=0;i<20;i++){
    	gf_inv(&x_inv_gf, &x_gf);
    	gf_inv(&x_gf, &x_inv_gf);
    	gf_inv(&x_inv_gf, &x_gf);
		gf_inv(&x_gf, &x_inv_gf);
    	gf_inv(&x_inv_gf, &x_gf);
		}

		t[j] = get_ticks() - t_begin;
		total_t += t[j]; 

		gf_to_mpz(x_inv_mpz, &x_inv_gf);
		if (!check_correctness(x_inv_mpz, x_mpz[j], p_mpz)){
			fprintf(stderr, "ERR: wrong reduction result using `gf_inv` of bingcd\n");
			exit(EXIT_FAILURE);
		}
	}
	qsort(t, test_count, sizeof(uint64_t), cmp_int64);
	benchmark_bingcd.best =  t[0]/100;
	benchmark_bingcd.median = t[test_count/2]/100;
	benchmark_bingcd.average =  total_t/(100 * (double)test_count);



	/* Test EEA_q */
	total_t = 0;
	for(size_t j=0; j<test_count; j++)
	{
		
		x_inv[0] = x_inv[1] = x_inv[2] = x_inv[3] = 0ul;
		mpz_2_u64(x, x_mpz[j], 4);
    	
		t_begin = get_ticks();
		for(int i=0;i<20;i++){
			inverse25519_EEA_vartime(x_inv, x);
			inverse25519_EEA_vartime(x, x_inv);
			inverse25519_EEA_vartime(x_inv, x);
			inverse25519_EEA_vartime(x, x_inv);
			inverse25519_EEA_vartime(x_inv, x);
		}
		t[j] = get_ticks() - t_begin;
		total_t += t[j]; 

		u64_2_mpz(x_inv_mpz, x_inv, 4);
		if (!check_correctness(x_inv_mpz, x_mpz[j], p_mpz)){
			fprintf(stderr, "ERR: wrong reduction result using `inverse25519_EEA_vartime`\n");
			exit(EXIT_FAILURE);
		}
	}
	qsort(t, test_count, sizeof(uint64_t), cmp_int64);
	benchmark_EEA.best =  t[0]/100;
	benchmark_EEA.median = t[test_count/2]/100;
	benchmark_EEA.average =  total_t/(100*(double)test_count);



	/* Test safegcd */
	total_t = 0;
	for(size_t j=0; j<test_count; j++)
	{
		gmp_export(x_c, 32, x_mpz[j]);
    	
		t_begin = get_ticks();
		for(int i=0;i<20;i++){
			inverse25519(x_inv_c, x_c);
			inverse25519(x_c, x_inv_c);
			inverse25519(x_inv_c, x_c);
			inverse25519(x_c, x_inv_c);
			inverse25519(x_inv_c, x_c);
		}
		t[j] = get_ticks() - t_begin;
		total_t += t[j]; 
		gmp_import(x_inv_mpz, x_inv_c, 32);
		if (!check_correctness(x_inv_mpz, x_mpz[j], p_mpz)){
			fprintf(stderr, "ERR: wrong reduction result using `inverse25519` in  safegcd\n");
			exit(EXIT_FAILURE);
		}
	}
	qsort(t, test_count, sizeof(uint64_t), cmp_int64);
	benchmark_safegcd.best =  t[0]/100;
	benchmark_safegcd.median = t[test_count/2]/100;
	benchmark_safegcd.average =  total_t/(100*(double)test_count);

	    
	printf("───────────────────────────────────────────────────────────────────────\n");
	printf("Time (ticks)|   binGCD     |     EEA_q    | Speed up     | Improvement\n");
	printf("───────────────────────────────────────────────────────────────────────\n");
	printf("Best        | %-12lu | %-12lu | %-12.4f | %.2f %%\n", benchmark_bingcd.best, benchmark_EEA.best, (double)benchmark_bingcd.best/(double)benchmark_EEA.best,((double)benchmark_bingcd.best - (double)benchmark_EEA.best)/((double)benchmark_bingcd.best) * 100);
	printf("Median      | %-12lu | %-12lu | %-12.4f | %.2f %%\n", benchmark_bingcd.median, benchmark_EEA.median, (double)benchmark_bingcd.median/(double)benchmark_EEA.median,((double)benchmark_bingcd.median - (double)benchmark_EEA.median)/((double)benchmark_bingcd.median) * 100);
	printf("Average     | %-12.2f | %-12.2f | %-12.4f | %.2f %%\n", benchmark_bingcd.average, benchmark_EEA.average, (double)benchmark_bingcd.average/(double)benchmark_EEA.average,((double)benchmark_bingcd.average - (double)benchmark_EEA.average)/((double)benchmark_bingcd.average) * 100);
	printf("───────────────────────────────────────────────────────────────────────\n");

	printf("───────────────────────────────────────────────────────────────────────\n");
	printf("Time (ticks)|   safeGCD    |     EEA_q    | Speed up     | Improvement\n");
	printf("───────────────────────────────────────────────────────────────────────\n");
	printf("Best        | %-12lu | %-12lu | %-12.4f | %.2f %%\n", benchmark_safegcd.best, benchmark_EEA.best, (double)benchmark_safegcd.best/(double)benchmark_EEA.best,((double)benchmark_safegcd.best - (double)benchmark_EEA.best)/((double)benchmark_safegcd.best) * 100);
	printf("Median      | %-12lu | %-12lu | %-12.4f | %.2f %%\n", benchmark_safegcd.median, benchmark_EEA.median, (double)benchmark_safegcd.median/(double)benchmark_EEA.median,((double)benchmark_safegcd.median - (double)benchmark_EEA.median)/((double)benchmark_safegcd.median) * 100);
	printf("Average     | %-12.2f | %-12.2f | %-12.4f | %.2f %%\n", benchmark_safegcd.average, benchmark_EEA.average, (double)benchmark_safegcd.average/(double)benchmark_EEA.average,((double)benchmark_safegcd.average - (double)benchmark_EEA.average)/((double)benchmark_safegcd.average) * 100);
	printf("───────────────────────────────────────────────────────────────────────\n");


	return 0;
}


int main(){
	
	printf("Benchmark of inverse25519:\n");
	printf("Number of samples = %i \n", number_of_samples);
	test_instance(number_of_samples);
	printf("Done!\n");
		
}
