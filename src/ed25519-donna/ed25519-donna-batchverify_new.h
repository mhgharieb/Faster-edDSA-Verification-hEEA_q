/*
 * Copyright (c) 2024 Muhammad ElSheikh (mhgaelsh@uwaterloo.ca). All rights reserved.
 * This file is part of the source code for the paper "Accelerating EdDSA Signature Verification with Faster Scalar Size Halving".
 * Licensed under the MIT License (for original modifications and new files).
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 * FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
 * THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * Ed25519 new batch verification
 * Updated version from ed25519-donna-batchverify.h
*/




// #define max_batch_size_new 64
// #define heap_batch_size_new max_batch_size_new
#define bignum256modm_128bitsScalar_most_significant_limb_index (128 / bignum256modm_bits_per_limb)

static void
heap_build_new(batch_heap *heap, size_t count) {
	heap->heap[0] = 0;
	heap->size = 0;
	while (heap->size < count){
		// heap_insert_next_new(heap);

		size_t node = heap->size, parent;
		heap_index_t *pheap = heap->heap;
		bignum256modm *scalars = heap->scalars;

		/* insert at the bottom */
		pheap[node] = (heap_index_t)node;

		/* sift node up to its sorted spot */
		parent = (node - 1) / 2;

		while (node && lt256_modm_batch(scalars[pheap[parent]], scalars[pheap[node]], bignum256modm_128bitsScalar_most_significant_limb_index)) {
			heap_swap(pheap, parent, node);
			node = parent;
			parent = (node - 1) / 2;
		}
		heap->size++;
		
	}

}

static void
ge25519_multi_scalarmult_vartime_final_new(ge25519 *r, ge25519 *point, bignum256modm scalar) {
	const bignum256modm_element_t topbit = ((bignum256modm_element_t)1 << (bignum256modm_bits_per_limb - 1));
	size_t limb = limb128bits;
	bignum256modm_element_t flag;

	if (isone256_modm_batch(scalar)) {
		/* this will happen most of the time after bos-carter */
		*r = *point;
		return;
	} else if (iszero256_modm_batch(scalar)) {
		/* this will only happen if all scalars == 0 */
		memset(r, 0, sizeof(*r));
		r->y[0] = 1;
		r->z[0] = 1;
		return;
	}

	*r = *point;

	/* find the limb where first bit is set */
	while (!scalar[limb])
		limb--;

	/* find the first bit */
	flag = topbit;
	while ((scalar[limb] & flag) == 0)
		flag >>= 1;

	/* exponentiate */
	flag >>= 1; // since *r = *p, pass the 1st iteration
	for (;;) {
		ge25519_double(r, r);
		if (scalar[limb] & flag)
			ge25519_add(r, r, point);

		flag >>= 1;
		if (!flag) {
			if (!limb--)
				break;
			flag = topbit;
		}
	}
}

/* count must be >= 3 */
static void
ge25519_multi_scalarmult_vartime_new(heap_index_t *max_index, batch_heap *heap, size_t count) {
	heap_index_t max1, max2;

	size_t limbsize =  bignum256modm_128bitsScalar_most_significant_limb_index;
	
	// heap_build(heap, ((count + 1) / 2) | 1);
	heap_build_new(heap, count);

	for (;;) {
	heap_get_top2(heap, &max1, &max2, limbsize);

	/* only one scalar remaining, we're done */
	if (iszero256_modm_batch(heap->scalars[max2]))
		break;

	/* exhausted another limb? */
	if (!heap->scalars[max1][limbsize])
		limbsize -= 1;

	/* can we extend to the 128 bit scalars? */
	// if (!extended && isatmost128bits256_modm_batch(heap->scalars[max1])) {
		// heap_extend(heap, count);
	// 	heap_get_top2(heap, &max1, &max2, limbsize);
		// extended = 1;
	// }

	sub256_modm_batch(heap->scalars[max1], heap->scalars[max1], heap->scalars[max2], limbsize);
	ge25519_add(&heap->points[max2], &heap->points[max2], &heap->points[max1]);
	heap_updated_root(heap, limbsize);
	}

	// ge25519_multi_scalarmult_vartime_new_final(r, &heap->points[max1], heap->scalars[max1]);
	*max_index = max1;
}


/* New batch verification using the new randamization with hEEA */
int
ED25519_FN(ed25519_sign_open_batch_hEEA) (const unsigned char **m, size_t *mlen, const unsigned char **pk, const unsigned char **RS, size_t num, int *valid) {
	heap_index_t max_index;
	batch_heap ALIGN(16) batch;
	ge25519 ALIGN(16) sumR;
	ge25519 ALIGN(16) sumBA;
	ge25519 ALIGN(16) sumBAR;
	ge25519 *tmp_point;
	bignum25519 tmp_cordinate;
	size_t i, batchsize, heap_size, limbsize = bignum256modm_limb_size - 1;
	unsigned char hram[64], U_chr[32];
	bignum256modm U, Uinv, h, v, TauS, sumTauS = {0}, Zero = {0}, rho[max_batch_size], tau_last, rho_last, S2={0};
	int tau_isNegative[max_batch_size], rho_isNegative[max_batch_size], tau_last_isNegative, rho_last_isNegative;
	int ret = 0;

	for (i = 0; i < num; i++)
		valid[i] = 1;

	while (num > 3) {
		batchsize = (num > max_batch_size) ? max_batch_size : num;

		/* heap_size must be odd */
		heap_size = batchsize + !(batchsize % 2);

		/* pick a random U, s.t. 0 < U < el */
		ed25519_randombytes_unsafe(U_chr, 32);
		expand256_modm(U, U_chr, 32);

		/* Uinv * U = 1 mod el */
		inv256_modm(Uinv, U);

		for (i = 0; i < batchsize; i++) {

			/* compute h <-- H(R_i,A_i,m_i) */
			ed25519_hram(hram, RS[i], pk[i], m[i], mlen[i]);
			expand256_modm(h, hram, 64);
			
			/* compute v <-- Uinv * h mod el */
			mul256_modm(v, Uinv, h);
			
			/* compute rho and tau s.t. tau * v = rho mod el ==> tau * h = U * rho mod el */
			curve25519_half_size_scalar_vartime_hEEA(rho[i], batch.scalars[i], v, &rho_isNegative[i], &tau_isNegative[i]);
			
			/* extract S */
			expand256_modm(TauS, RS[i] + 32, 32);
			
			/* compute tau * S */
			mul256_modm(TauS, TauS, batch.scalars[i]);

			/* compute sum of tau * S */
			if (tau_isNegative[i])
				sub256_modm_batch(TauS, modm_m, TauS, limbsize);
			add256_modm(sumTauS, sumTauS, TauS);
		};

		/* unpacking (-R_i) and adjust the sign based on the sign of tau_i */
		for (i = 0; i < batchsize; i++){
			if (tau_isNegative[i]){
				if (!ge25519_unpack_positive_vartime(&batch.points[i], RS[i]))
					goto fallback;
			}else{
				if (!ge25519_unpack_negative_vartime(&batch.points[i], RS[i]))
					goto fallback;
			}
		}
		/* if batchsize is even, add (point, scalar) = (neutral,{0}) to the heap */
		if (heap_size != batchsize){
			memcpy(batch.scalars[batchsize], Zero, sizeof(Zero));
			tmp_point = &batch.points[batchsize];
			memset(tmp_point, 0, sizeof(*tmp_point));
			tmp_point->y[0] = 1;
			tmp_point->z[0] = 1;
		}

		/* comute sumR <-- sum([tau_i](-R_i)) */
		ge25519_multi_scalarmult_vartime_new(&max_index, &batch, heap_size);
		ge25519_multi_scalarmult_vartime_final_new(&sumR, &batch.points[max_index], batch.scalars[max_index]);

		/* unpacking (-A_i) and adjust the sign based on the sign of rho_i */
		for (i = 0; i < batchsize; i++){
			if (rho_isNegative[i]){
				if (!ge25519_unpack_positive_vartime(&batch.points[i], pk[i]))
				goto fallback;
			}else{
				if (!ge25519_unpack_negative_vartime(&batch.points[i], pk[i]))
				goto fallback;
			}
		}
		
		/* Assign rho to batch.scalars */
		for (i = 0; i < batchsize; i++) {
			for (size_t j = 0; j<bignum256modm_limb_size; j++)
				batch.scalars[i][j] = rho[i][j];
		}

		/* if batchsize is even, add (point, scalar) = (neutral,{0}) to the heap */
		if (heap_size != batchsize){
			memcpy(batch.scalars[batchsize], Zero, sizeof(Zero));
			tmp_point = &batch.points[batchsize];
			memset(tmp_point, 0, sizeof(*tmp_point));
			tmp_point->y[0] = 1;
			tmp_point->z[0] = 1;
		}

		/* comput sumA = [k]P = sum([rho_i](-A_i)) */
		ge25519_multi_scalarmult_vartime_new(&max_index, &batch, heap_size);
		
		/* compute U <-- Uk = U * k mod el */
		mul256_modm(U, U, batch.scalars[max_index]);

		/* compute rho_last and tau_last s.t. tau_last * Uk = rho_last mod el */
		curve25519_half_size_scalar_vartime_hEEA(rho_last, tau_last, U, &rho_last_isNegative, &tau_last_isNegative);
		
		/*
		[sumTauS]B + sumR + [U]sumA =? 0  ==> [tau_last * sumTauS]B + [tau_last]sumR + [rho_last]P =? 0
		----------------------------------------------------------------------
		tau_last_isNegative | rho_last_isNegative | [tau_last * sumTauS]B + [tau_last]sumR + [rho_last]P
		----------------------------------------------------------------------
		         0          |         0           | [|tau_last| * sumTauS]B + [|tau_last|]sumR + [|rho_last|]P
		----------------------------------------------------------------------
		         0          |         1           | [|tau_last| * sumTauS]B + [|tau_last|]sumR - [|rho_last|]P
		----------------------------------------------------------------------
		         1          |         0           | [-|tau_last| * sumTauS]B + [-|tau_last|]sumR + [|rho_last|]P ==> 
			  		        |                     | [|tau_last| * sumTauS] B + [|tau_last|]sumR - [|rho_last|]P
		----------------------------------------------------------------------
		         1          |         1           | [-|tau_last| * sumTauS]B + [-|tau_last|]sumR + [-|rho_last|]P ==> 
			  	        	|                     | [|tau_last| * sumTauS] B + [|tau_last|]sumR + [|rho_last|]P
		----------------------------------------------------------------------
		*/
		/* Adjust the sign of P based on the previous table*/
		if (tau_last_isNegative != rho_last_isNegative){
			/* P <-- -P */
			tmp_point = &batch.points[max_index];
			curve25519_copy(tmp_cordinate, tmp_point->x);
			curve25519_neg(tmp_point->x, tmp_cordinate);
			curve25519_copy(tmp_cordinate, tmp_point->t);
			curve25519_neg(tmp_point->t, tmp_cordinate);
		}

		/* sumTauS <-- tau_last * sumTauS mod el */
		mul256_modm(sumTauS, sumTauS, tau_last);

		/* split S to S1 and S2, s.t. S = (S2<<126 + S1), bl(S1) = 126, bl(S2) = 127 */
		S2[0] = (sumTauS[2] >> 14) | ((sumTauS[3] & 0x3FFF) << 42);
		S2[1] = (sumTauS[3] >> 14) | ((sumTauS[4] & 0x3FFF) << 42);
		S2[2] = sumTauS[4] >> 14;
		sumTauS[2] &= 0x3FFF;
		sumTauS[3] = 0;
		sumTauS[4] = 0;

		/*
		* [sumTauS]B + sumR + [U]sumA =? 0 ==> 
		* [sumTauS]B + sumR + [Uk]P =? 0 ==>
		* [tau_last * sumTauS]B + [tau_last]sumR + [rho_last]P =? 0 ==>
		* [S1]B + [S2]([2^126]B) + [tau_last]sumR + [rho_last]P =? 0 
		*/
		ge25519_quadruple_scalarmult_vartime(&sumBAR, &sumR, &batch.points[max_index], tau_last, rho_last, sumTauS, S2);

		/* Check if the sum is 0 */
		if (!ge25519_is_neutral_vartime(&sumBAR)) {
			ret |= 2;

			fallback:
			for (i = 0; i < batchsize; i++) {
				valid[i] = ed25519_sign_open_hEEA (m[i], mlen[i], pk[i], RS[i]) ? 0 : 1;
				ret |= (valid[i] ^ 1);
			}
		}

		m += batchsize;
		mlen += batchsize;
		pk += batchsize;
		RS += batchsize;
		num -= batchsize;
		valid += batchsize;
	}

	for (i = 0; i < num; i++) {
		valid[i] = ed25519_sign_open_hEEA (m[i], mlen[i], pk[i], RS[i]) ? 0 : 1;
		ret |= (valid[i] ^ 1);
	}

	return ret;
}

/* New batch verification using the new randamization with hgcd_enhance2 */
int
ED25519_FN(ed25519_sign_open_batch_hgcd) (const unsigned char **m, size_t *mlen, const unsigned char **pk, const unsigned char **RS, size_t num, int *valid) {
	heap_index_t max_index;
	batch_heap ALIGN(16) batch;
	ge25519 ALIGN(16) sumR;
	ge25519 ALIGN(16) sumBA;
	ge25519 ALIGN(16) sumBAR;
	ge25519 *tmp_point;
	bignum25519 tmp_cordinate;
	size_t i, batchsize, heap_size, limbsize = bignum256modm_limb_size - 1;
	unsigned char hram[64], U_chr[32];
	bignum256modm U, Uinv, h, v, TauS, sumTauS = {0}, Zero = {0}, rho[max_batch_size], tau_last, rho_last, S2={0};
	int tau_isNegative[max_batch_size], rho_isNegative[max_batch_size], tau_last_isNegative, rho_last_isNegative;
	int ret = 0;

	for (i = 0; i < num; i++)
		valid[i] = 1;

	while (num > 3) {
		batchsize = (num > max_batch_size) ? max_batch_size : num;

		/* heap_size must be odd */
		heap_size = batchsize + !(batchsize % 2);

		/* pick a random U, s.t. 0 < U < el */
		ed25519_randombytes_unsafe(U_chr, 32);
		expand256_modm(U, U_chr, 32);

		/* Uinv * U = 1 mod el */
		inv256_modm(Uinv, U);

		for (i = 0; i < batchsize; i++) {

			/* compute h <-- H(R_i,A_i,m_i) */
			ed25519_hram(hram, RS[i], pk[i], m[i], mlen[i]);
			expand256_modm(h, hram, 64);
			
			/* compute v <-- Uinv * h mod el */
			mul256_modm(v, Uinv, h);
			
			/* compute rho and tau s.t. tau * v = rho mod el ==> tau * h = U * rho mod el */
			// curve25519_half_size_scalar_vartime_hEEA(rho[i], batch.scalars[i], v, &rho_isNegative[i], &tau_isNegative[i]);
			curve25519_half_size_scalar_vartime_hgcd(rho[i], batch.scalars[i], v, &rho_isNegative[i], &tau_isNegative[i]);
			
			/* extract S */
			expand256_modm(TauS, RS[i] + 32, 32);
			
			/* compute tau * S */
			mul256_modm(TauS, TauS, batch.scalars[i]);

			/* compute sum of tau * S */
			if (tau_isNegative[i])
				sub256_modm_batch(TauS, modm_m, TauS, limbsize);
			add256_modm(sumTauS, sumTauS, TauS);
		};

		/* unpacking (-R_i) and adjust the sign based on the sign of tau_i */
		for (i = 0; i < batchsize; i++){
			if (tau_isNegative[i]){
				if (!ge25519_unpack_positive_vartime(&batch.points[i], RS[i]))
					goto fallback;
			}else{
				if (!ge25519_unpack_negative_vartime(&batch.points[i], RS[i]))
					goto fallback;
			}
		}
		/* if batchsize is even, add (point, scalar) = (neutral,{0}) to the heap */
		if (heap_size != batchsize){
			memcpy(batch.scalars[batchsize], Zero, sizeof(Zero));
			tmp_point = &batch.points[batchsize];
			memset(tmp_point, 0, sizeof(*tmp_point));
			tmp_point->y[0] = 1;
			tmp_point->z[0] = 1;
		}

		/* comute sumR <-- sum([tau_i](-R_i)) */
		ge25519_multi_scalarmult_vartime_new(&max_index, &batch, heap_size);
		ge25519_multi_scalarmult_vartime_final_new(&sumR, &batch.points[max_index], batch.scalars[max_index]);

		/* unpacking (-A_i) and adjust the sign based on the sign of rho_i */
		for (i = 0; i < batchsize; i++){
			if (rho_isNegative[i]){
				if (!ge25519_unpack_positive_vartime(&batch.points[i], pk[i]))
				goto fallback;
			}else{
				if (!ge25519_unpack_negative_vartime(&batch.points[i], pk[i]))
				goto fallback;
			}
		}
		
		/* Assign rho to batch.scalars */
		for (i = 0; i < batchsize; i++) {
			for (size_t j = 0; j<bignum256modm_limb_size; j++)
				batch.scalars[i][j] = rho[i][j];
		}

		/* if batchsize is even, add (point, scalar) = (neutral,{0}) to the heap */
		if (heap_size != batchsize){
			memcpy(batch.scalars[batchsize], Zero, sizeof(Zero));
			tmp_point = &batch.points[batchsize];
			memset(tmp_point, 0, sizeof(*tmp_point));
			tmp_point->y[0] = 1;
			tmp_point->z[0] = 1;
		}

		/* comput sumA = [k]P = sum([rho_i](-A_i)) */
		ge25519_multi_scalarmult_vartime_new(&max_index, &batch, heap_size);
		
		/* compute U <-- Uk = U * k mod el */
		mul256_modm(U, U, batch.scalars[max_index]);

		/* compute rho_last and tau_last s.t. tau_last * Uk = rho_last mod el */
		// curve25519_half_size_scalar_vartime_hEEA(rho_last, tau_last, U, &rho_last_isNegative, &tau_last_isNegative);
		curve25519_half_size_scalar_vartime_hgcd(rho_last, tau_last, U, &rho_last_isNegative, &tau_last_isNegative);
		
		/*
		[sumTauS]B + sumR + [U]sumA =? 0  ==> [tau_last * sumTauS]B + [tau_last]sumR + [rho_last]P =? 0
		----------------------------------------------------------------------
		tau_last_isNegative | rho_last_isNegative | [tau_last * sumTauS]B + [tau_last]sumR + [rho_last]P
		----------------------------------------------------------------------
		         0          |         0           | [|tau_last| * sumTauS]B + [|tau_last|]sumR + [|rho_last|]P
		----------------------------------------------------------------------
		         0          |         1           | [|tau_last| * sumTauS]B + [|tau_last|]sumR - [|rho_last|]P
		----------------------------------------------------------------------
		         1          |         0           | [-|tau_last| * sumTauS]B + [-|tau_last|]sumR + [|rho_last|]P ==> 
			  		        |                     | [|tau_last| * sumTauS] B + [|tau_last|]sumR - [|rho_last|]P
		----------------------------------------------------------------------
		         1          |         1           | [-|tau_last| * sumTauS]B + [-|tau_last|]sumR + [-|rho_last|]P ==> 
			  	        	|                     | [|tau_last| * sumTauS] B + [|tau_last|]sumR + [|rho_last|]P
		----------------------------------------------------------------------
		*/
		/* Adjust the sign of P based on the previous table*/
		if (tau_last_isNegative != rho_last_isNegative){
			/* P <-- -P */
			tmp_point = &batch.points[max_index];
			curve25519_copy(tmp_cordinate, tmp_point->x);
			curve25519_neg(tmp_point->x, tmp_cordinate);
			curve25519_copy(tmp_cordinate, tmp_point->t);
			curve25519_neg(tmp_point->t, tmp_cordinate);
		}

		/* sumTauS <-- tau_last * sumTauS mod el */
		mul256_modm(sumTauS, sumTauS, tau_last);

		/* split S to S1 and S2, s.t. S = (S2<<126 + S1), bl(S1) = 126, bl(S2) = 127 */
		S2[0] = (sumTauS[2] >> 14) | ((sumTauS[3] & 0x3FFF) << 42);
		S2[1] = (sumTauS[3] >> 14) | ((sumTauS[4] & 0x3FFF) << 42);
		S2[2] = sumTauS[4] >> 14;
		sumTauS[2] &= 0x3FFF;
		sumTauS[3] = 0;
		sumTauS[4] = 0;

		/*
		* [sumTauS]B + sumR + [U]sumA =? 0 ==> 
		* [sumTauS]B + sumR + [Uk]P =? 0 ==>
		* [tau_last * sumTauS]B + [tau_last]sumR + [rho_last]P =? 0 ==>
		* [S1]B + [S2]([2^126]B) + [tau_last]sumR + [rho_last]P =? 0 
		*/
		ge25519_quadruple_scalarmult_vartime(&sumBAR, &sumR, &batch.points[max_index], tau_last, rho_last, sumTauS, S2);

		/* Check if the sum is 0 */
		if (!ge25519_is_neutral_vartime(&sumBAR)) {
			ret |= 2;

			fallback:
			for (i = 0; i < batchsize; i++) {
				// valid[i] = ed25519_sign_open_hEEA (m[i], mlen[i], pk[i], RS[i]) ? 0 : 1;
				valid[i] = ed25519_sign_open_hgcd (m[i], mlen[i], pk[i], RS[i]) ? 0 : 1;
				ret |= (valid[i] ^ 1);
			}
		}

		m += batchsize;
		mlen += batchsize;
		pk += batchsize;
		RS += batchsize;
		num -= batchsize;
		valid += batchsize;
	}

	for (i = 0; i < num; i++) {
		// valid[i] = ed25519_sign_open_hEEA (m[i], mlen[i], pk[i], RS[i]) ? 0 : 1;
		valid[i] = ed25519_sign_open_hgcd (m[i], mlen[i], pk[i], RS[i]) ? 0 : 1;
		ret |= (valid[i] ^ 1);
	}

	return ret;
}