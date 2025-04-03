/*
 * Copyright (c) 2019 Thomas Pornin
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
 * Modified from https://github.com/pornin/curve9767/blob/master/src/scalar_amd64.c
*/

#include "curve448_reduce_basis_vartime.h"

/*
 * We use the intrinsic functions:
 *   _lzcnt_u32(), _lzcnt_u64(), _addcarry_u64(), _subborrow_u64().
 */
#include <immintrin.h>

/*
 * For lattice basis reduction, we use an other representation with
 * 64-bit limbs, in 64-bit words (no padding bit). "Small integers"
 * fit on 7 limbs, "large integers" use 14 limbs. Values are signed;
 * two's complement is used for negative values.
 */

/*
 * Multiply two small values, result in a large value. The two source
 * values MUST be nonnegative and in the proper range.
 */
static void
mul(uint64_t *d, const uint64_t *a, const uint64_t *b)
{
	size_t i, j;
	uint64_t t[14];

	memset(t, 0, sizeof t);
	for (i = 0; i < 7; i ++) {
		uint64_t cc;

		cc = 0;
		for (j = 0; j < 7; j ++) {
			unsigned __int128 z;

			z = (unsigned __int128)a[i] * (unsigned __int128)b[j]
				+ (unsigned __int128)t[i + j] + cc;
			t[i + j] = (uint64_t)z;
			cc = (uint64_t)(z >> 64);
		}
		t[i + 7] = cc;
	}
	memcpy(d, t, sizeof t);
}

/*
 * This macro defines a function with prototype:
 *  static inline void name(uint64_t *a, const uint64_t *b, unsigned s)
 * It adds lshift(b,s) to a (if intrinsic_op is _addcarry_u64), or
 * subtracts lshift(b,s) from a (if intrinsic_op is _subborrow_u64).
 * Both values consist of 'size' limbs. Truncation happens if the result
 * does not fit in the output.
 */
#define DEF_OP_LSHIFT(name, intrinsic_op, size) \
static void \
name(uint64_t *a, const uint64_t *b, unsigned s) \
{ \
	uint64_t b2[(size)]; \
	unsigned char cc; \
	unsigned long long w, w2, e; \
	size_t i; \
 \
 	if (s >= 64) { \
		unsigned k; \
 \
		k = s >> 6; \
		s &= 63; \
		if (k >= (size)) { \
			return; \
		} \
		memset(b2, 0, k * sizeof b2[0]); \
		memcpy(b2 + k, b, ((size) - k) * sizeof b2[0]); \
		b = b2; \
	} \
	if (s == 0) { \
		cc = 0; \
		for (i = 0; i < (size); i ++) { \
			cc = intrinsic_op(cc, a[i], b[i], &w); \
			a[i] = w; \
		} \
	} else { \
		cc = 0; \
		e = 0; \
		for (i = 0; i < (size); i ++) { \
			w = b[i]; \
			cc = intrinsic_op(cc, a[i], (w << s) | e, &w2); \
			e = w >> (64 - s); \
			a[i] = w2; \
		} \
	} \
}

DEF_OP_LSHIFT(add_lshift_small, _addcarry_u64, 4)
DEF_OP_LSHIFT(sub_lshift_small, _subborrow_u64, 4)
DEF_OP_LSHIFT(add_lshift_large, _addcarry_u64, 14)
DEF_OP_LSHIFT(sub_lshift_large, _subborrow_u64, 14)



/*
 * Given a scalar b, compute signed integers c0 and c1 with the following
 * characteristics:
 *   - c1 != 0
 *   - c0 = c1 * b mod n
 *   - |c0| < 2^223
 *   - |c1| < 2^223
 *
 * THIS FUNCTION IS NOT CONSTANT-TIME. It is meant to be used as part
 * of Schnorr signature verification, when the signature value, the
 * signed message, and the public key, are all assumed to be public.
 */

// void
// curve9767_inner_reduce_basis_vartime(
// 	uint8_t *c0, uint8_t *c1, const curve9767_scalar *b)
void
curve448_reduce_basis_vartime(
	uint64_t *c0, uint64_t *c1, const uint64_t *b)
{
	unsigned long long u0_0, u0_1, u0_2, u0_3;
	unsigned long long u1_0, u1_1, u1_2, u1_3;
	unsigned long long v0_0, v0_1, v0_2, v0_3;
	unsigned long long v1_0, v1_1, v1_2, v1_3;
	unsigned long long nu_0, nu_1, nu_2, nu_3, nu_4, nu_5, nu_6, nu_7, nu_8, nu_9, nu_10, nu_11, nu_12, nu_13;
	unsigned long long nv_0, nv_1, nv_2, nv_3, nv_4, nv_5, nv_6, nv_7, nv_8, nv_9, nv_10, nv_11, nv_12, nv_13;
	unsigned long long sp_0, sp_1, sp_2, sp_3, sp_4, sp_5, sp_6, sp_7, sp_8, sp_9, sp_10, sp_11, sp_12, sp_13;
	unsigned long long w;
	uint64_t tmp_b[7], tmp_l[14];
	unsigned char cc;

	
	static const uint64_t order_u64[] = {
		0x2378c292ab5844f3,
		0x216cc2728dc58f55,
		0xc44edb49aed63690,
		0xffffffff7cca23e9,
		0xffffffffffffffff,
		0xffffffffffffffff,
		0x3fffffffffffffff
	};

	/*
	 * Normalize b to 0..n-1.
	 */
	// scalar_normalize(bw, b->v.w16);
	// to_small(tmp_b, bw);
	// memcpy(tmp_b, b, sizeof(b));
	for(int i =0; i<7;i++) {
			tmp_b[i] = b[i];
	}

	/*
	 * Init:
	 *   u = [n, 0]
	 *   v = [b, 1]
	 */
	u0_0 = order_u64[0];
	u0_1 = order_u64[1];
	u0_2 = order_u64[2];
	u0_3 = order_u64[3];
	u1_0 = 0;
	u1_1 = 0;
	u1_2 = 0;
	u1_3 = 0;
	v0_0 = tmp_b[0];
	v0_1 = tmp_b[1];
	v0_2 = tmp_b[2];
	v0_3 = tmp_b[3];
	v1_0 = 1;
	v1_1 = 0;
	v1_2 = 0;
	v1_3 = 0;

	/*
	 * nu = <u, u> = n^2
	 * nv = <v, v> = b^2 + 1
	 * sp = <u, v> = n*b
	 */
	

	nu_0 = 0xc1adfbb81ba1fea9;
	nu_1 = 0xb91bf53749e0a8b2;
	nu_2 = 0x4525492be764d815;
	nu_3 = 0x4ae17cf6a2b8716d;
	nu_4 = 0xf1a9cc14ba3c47c4;
	nu_5 = 0x92052bcb7e4d070a;
	nu_6 = 0xc3402a939f823b72;
	nu_7 = 0x91bc614955ac2279;
	nu_8 = 0x10b6613946e2c7aa;
	nu_9 = 0xe2276da4d76b1b48;
	nu_10 = 0xffffffffbe6511f4;
	nu_11 = 0xffffffffffffffff;
	nu_12 = 0xffffffffffffffff;
	nu_13 = 0xfffffffffffffff;

	mul(tmp_l, tmp_b, tmp_b);
	
	cc = _addcarry_u64(0, tmp_l[0], 1, &nv_0);
	cc = _addcarry_u64(cc, tmp_l[1], 0, &nv_1);
	cc = _addcarry_u64(cc, tmp_l[2], 0, &nv_2);
	cc = _addcarry_u64(cc, tmp_l[3], 0, &nv_3);
	cc = _addcarry_u64(cc, tmp_l[4], 0, &nv_4);
	cc = _addcarry_u64(cc, tmp_l[5], 0, &nv_5);
	cc = _addcarry_u64(cc, tmp_l[6], 0, &nv_6);
	cc = _addcarry_u64(cc, tmp_l[7], 0, &nv_7);
	cc = _addcarry_u64(cc, tmp_l[8], 0, &nv_8);
	cc = _addcarry_u64(cc, tmp_l[9], 0, &nv_9);
	cc = _addcarry_u64(cc, tmp_l[10], 0, &nv_10);
	cc = _addcarry_u64(cc, tmp_l[11], 0, &nv_11);
	cc = _addcarry_u64(cc, tmp_l[12], 0, &nv_12);
	(void)_addcarry_u64(cc, tmp_l[13], 0, &nv_13);
		

	mul(tmp_l, order_u64, tmp_b);

	sp_0 = tmp_l[0];
	sp_1 = tmp_l[1];
	sp_2 = tmp_l[2];
	sp_3 = tmp_l[3];
	sp_4 = tmp_l[4];
	sp_5 = tmp_l[5];
	sp_6 = tmp_l[6];
	sp_7 = tmp_l[7];
	sp_8 = tmp_l[8];
	sp_9 = tmp_l[9];
	sp_10 = tmp_l[10];
	sp_11 = tmp_l[11];
	sp_12 = tmp_l[12];
	sp_13 = tmp_l[13];

	/*
	 * Algorithm:
	 *   - If nu < nv, Then:
	 *         swap(u, v)
	 *         swap(nu, nv)
	 *   - If bitlength(nv) <= target, Then:
	 *         c0 <- v0
	 *         c1 <- v1
	 *         return success
	 *   - Set s <- max(0, bitlength(sp) - bitlength(nv))
	 *   - If sp > 0, Then:
	 *         u <- u - lshift(v, s)
	 *         nu <- nu + lshift(nv, 2*s) - lshift(sp, s+1)
	 *         sp <- sp - lshift(nv, s)
	 *     Else:
	 *         u <- u + lshift(v, s)
	 *         nu <- nu + lshift(nv, 2*s) + lshift(sp, s+1)
	 *         sp <- sp + lshift(nv, s)
	 */
	for (;;) {
		unsigned s;
		unsigned long long m;
		unsigned bl_nv, bl_sp;

		/*
		 * If nu < nv, then swap(u,v) and swap(nu,nv).
		 */
		cc = _subborrow_u64(0, nu_0, nv_0, &w);
		cc = _subborrow_u64(cc, nu_1, nv_1, &w);
		cc = _subborrow_u64(cc, nu_2, nv_2, &w);
		cc = _subborrow_u64(cc, nu_3, nv_3, &w);
		cc = _subborrow_u64(cc, nu_4, nv_4, &w);
		cc = _subborrow_u64(cc, nu_5, nv_5, &w);
		cc = _subborrow_u64(cc, nu_6, nv_6, &w);
		cc = _subborrow_u64(cc, nu_7, nv_7, &w);
		cc = _subborrow_u64(cc, nu_8, nv_8, &w);
		cc = _subborrow_u64(cc, nu_9, nv_9, &w);
		cc = _subborrow_u64(cc, nu_10, nv_10, &w);
		cc = _subborrow_u64(cc, nu_11, nv_11, &w);
		cc = _subborrow_u64(cc, nu_12, nv_12, &w);
		cc = _subborrow_u64(cc, nu_13, nv_13, &w);
		m = -(unsigned long long)cc;


#define SWAPCOND(x, y, ctlmask)   do { \
		unsigned long long swap_mask = (ctlmask) & ((x) ^ (y)); \
		(x) ^= swap_mask; \
		(y) ^= swap_mask; \
	} while (0)

		SWAPCOND(u0_0, v0_0, m);
		SWAPCOND(u0_1, v0_1, m);
		SWAPCOND(u0_2, v0_2, m);
		SWAPCOND(u0_3, v0_3, m);
		
		SWAPCOND(u1_0, v1_0, m);
		SWAPCOND(u1_1, v1_1, m);
		SWAPCOND(u1_2, v1_2, m);
		SWAPCOND(u1_3, v1_3, m);

		SWAPCOND(nu_0, nv_0, m);
		SWAPCOND(nu_1, nv_1, m);
		SWAPCOND(nu_2, nv_2, m);
		SWAPCOND(nu_3, nv_3, m);
		SWAPCOND(nu_4, nv_4, m);
		SWAPCOND(nu_5, nv_5, m);
		SWAPCOND(nu_6, nv_6, m);
		SWAPCOND(nu_7, nv_7, m);
		SWAPCOND(nu_8, nv_8, m);
		SWAPCOND(nu_9, nv_9, m);
		SWAPCOND(nu_10, nv_10, m);
		SWAPCOND(nu_11, nv_11, m);
		SWAPCOND(nu_12, nv_12, m);
		SWAPCOND(nu_13, nv_13, m);

		/*
		 * u is now the largest vector; its square norm is nu.
		 * We know that:
		 *   N(u-v) = N(u) + N(v) - 2*<u,v>
		 *   N(u+v) = N(u) + N(v) + 2*<u,v>
		 * Since all squared norms are positive, and since
		 * N(u) >= N(v), it is guaranteed that |sp| <= nu.
		 * Thus, if nu fits on 12 words (including the sign bit),
		 * then we can jump to phase 2, where 'large' values use
		 * 12 words.
		 */
		if ((nu_13 | nu_12) == 0 && nu_11 <= 0x7FFFFFFFFFFFFFFFull) {
			break;
		}

#define BITLENGTH(size, bb)   do { \
		unsigned bitlength_acc = 14; \
		int bitlength_flag = 1; \
		unsigned long long bitlength_mask = -(bb ## 13 >> 63); \
		unsigned long long bitlength_top = 0; \
		unsigned long long bitlength_word; \
		bitlength_word = bb ## 13 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 12 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 11 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 10 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 9 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 8 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 7 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 6 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 5 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 4 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 3 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 2 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 1 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 0 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		(size) = (bitlength_acc << 6) - _lzcnt_u64(bitlength_top); \
	} while (0)

		BITLENGTH(bl_nv, nv_);
		BITLENGTH(bl_sp, sp_);

		/*
		 * If v is small enough, return.
		 */
		if (bl_nv <= 447) {
			c0[0] = v0_0;
			c0[1] = v0_1;
			c0[2] = v0_2;
			c0[3] = v0_3;
			
			c1[0] = v1_0;
			c1[1] = v1_1;
			c1[2] = v1_2;
			c1[3] = v1_3;

			return;
		}

		/*
		 * Compute shift amount.
		 */
		s = bl_sp - bl_nv;
		s &= ~-(s >> 31);

		/*
		 * It is very rare that s > 31. We handle it with some
		 * generic code; branch prediction will soon learn that
		 * this path is normally not taken.
		 */
		if (s > 31) {
			uint64_t tu0[4], tu1[4], tv0[4], tv1[4];
			uint64_t tnu[14], tnv[14], tsp[14];

			tu0[0] = u0_0;
			tu0[1] = u0_1;
			tu0[2] = u0_2;
			tu0[3] = u0_3;

			tu1[0] = u1_0;
			tu1[1] = u1_1;
			tu1[2] = u1_2;
			tu1[3] = u1_3;

			tv0[0] = v0_0;
			tv0[1] = v0_1;
			tv0[2] = v0_2;
			tv0[3] = v0_3;

			tv1[0] = v1_0;
			tv1[1] = v1_1;
			tv1[2] = v1_2;
			tv1[3] = v1_3;

			tnu[0] = nu_0;
			tnu[1] = nu_1;
			tnu[2] = nu_2;
			tnu[3] = nu_3;
			tnu[4] = nu_4;
			tnu[5] = nu_5;
			tnu[6] = nu_6;
			tnu[7] = nu_7;
			tnu[8] = nu_8;
			tnu[9] = nu_9;
			tnu[10] = nu_10;
			tnu[11] = nu_11;
			tnu[12] = nu_12;
			tnu[13] = nu_13;

			tnv[0] = nv_0;
			tnv[1] = nv_1;
			tnv[2] = nv_2;
			tnv[3] = nv_3;
			tnv[4] = nv_4;
			tnv[5] = nv_5;
			tnv[6] = nv_6;
			tnv[7] = nv_7;
			tnv[8] = nv_8;
			tnv[9] = nv_9;
			tnv[10] = nv_10;
			tnv[11] = nv_11;
			tnv[12] = nv_12;
			tnv[13] = nv_13;

			tsp[0] = sp_0;
			tsp[1] = sp_1;
			tsp[2] = sp_2;
			tsp[3] = sp_3;
			tsp[4] = sp_4;
			tsp[5] = sp_5;
			tsp[6] = sp_6;
			tsp[7] = sp_7;
			tsp[8] = sp_8;
			tsp[9] = sp_9;
			tsp[10] = sp_10;
			tsp[11] = sp_11;
			tsp[12] = sp_12;
			tsp[13] = sp_13;

			if ((sp_13 >> 63) == 0) {
				sub_lshift_small(tu0, tv0, s);
				sub_lshift_small(tu1, tv1, s);
				add_lshift_large(tnu, tnv, 2 * s);
				sub_lshift_large(tnu, tsp, s + 1);
				sub_lshift_large(tsp, tnv, s);
			} else {
				add_lshift_small(tu0, tv0, s);
				add_lshift_small(tu1, tv1, s);
				add_lshift_large(tnu, tnv, 2 * s);
				add_lshift_large(tnu, tsp, s + 1);
				add_lshift_large(tsp, tnv, s);
			}


			u0_0 = tu0[0];
			u0_1 = tu0[1];
			u0_2 = tu0[2];
			u0_3 = tu0[3];

			u1_0 = tu1[0];
			u1_1 = tu1[1];
			u1_2 = tu1[2];
			u1_3 = tu1[3];

			v0_0 = tv0[0];
			v0_1 = tv0[1];
			v0_2 = tv0[2];
			v0_3 = tv0[3];

			v1_0 = tv1[0];
			v1_1 = tv1[1];
			v1_2 = tv1[2];
			v1_3 = tv1[3];

			nu_0 = tnu[0];
			nu_1 = tnu[1];
			nu_2 = tnu[2];
			nu_3 = tnu[3];
			nu_4 = tnu[4];
			nu_5 = tnu[5];
			nu_6 = tnu[6];
			nu_7 = tnu[7];
			nu_8 = tnu[8];
			nu_9 = tnu[9];
			nu_10 = tnu[10];
			nu_11 = tnu[11];
			nu_12 = tnu[12];
			nu_13 = tnu[13];

			

			nv_0 = tnv[0];
			nv_1 = tnv[1];
			nv_2 = tnv[2];
			nv_3 = tnv[3];
			nv_4 = tnv[4];
			nv_5 = tnv[5];
			nv_6 = tnv[6];
			nv_7 = tnv[7];
			nv_8 = tnv[8];
			nv_9 = tnv[9];
			nv_10 = tnv[10];
			nv_11 = tnv[11];
			nv_12 = tnv[12];
			nv_13 = tnv[13];

			sp_0 = tsp[0];
			sp_1 = tsp[1];
			sp_2 = tsp[2];
			sp_3 = tsp[3];
			sp_4 = tsp[4];
			sp_5 = tsp[5];
			sp_6 = tsp[6];
			sp_7 = tsp[7];
			sp_8 = tsp[8];
			sp_9 = tsp[9];
			sp_10 = tsp[10];
			sp_11 = tsp[11];
			sp_12 = tsp[12];
			sp_13 = tsp[13];
			continue;
		}

		if ((sp_13 >> 63) == 0) {
			if (s == 0) {
				cc = _subborrow_u64(0, u0_0, v0_0, &u0_0);
				cc = _subborrow_u64(cc, u0_1, v0_1, &u0_1);
				cc = _subborrow_u64(cc, u0_2, v0_2, &u0_2);
				(void)_subborrow_u64(cc, u0_3, v0_3, &u0_3);

				cc = _subborrow_u64(0, u1_0, v1_0, &u1_0);
				cc = _subborrow_u64(cc, u1_1, v1_1, &u1_1);
				cc = _subborrow_u64(cc, u1_2, v1_2, &u1_2);
				(void)_subborrow_u64(cc, u1_3, v1_3, &u1_3);

				cc = _addcarry_u64(0, nu_0, nv_0, &nu_0);
				cc = _addcarry_u64(cc, nu_1, nv_1, &nu_1);
				cc = _addcarry_u64(cc, nu_2, nv_2, &nu_2);
				cc = _addcarry_u64(cc, nu_3, nv_3, &nu_3);
				cc = _addcarry_u64(cc, nu_4, nv_4, &nu_4);
				cc = _addcarry_u64(cc, nu_5, nv_5, &nu_5);
				cc = _addcarry_u64(cc, nu_6, nv_6, &nu_6);
				cc = _addcarry_u64(cc, nu_7, nv_7, &nu_7);
				cc = _addcarry_u64(cc, nu_8, nv_8, &nu_8);
				cc = _addcarry_u64(cc, nu_9, nv_9, &nu_9);
				cc = _addcarry_u64(cc, nu_10, nv_10, &nu_10);
				cc = _addcarry_u64(cc, nu_11, nv_11, &nu_11);
				cc = _addcarry_u64(cc, nu_12, nv_12, &nu_12);
				(void)_addcarry_u64(cc, nu_13, nv_13, &nu_13);

				cc = _subborrow_u64(0, nu_0,
					(sp_0 << 1), &nu_0);
				cc = _subborrow_u64(cc, nu_1,
					(sp_1 << 1) | (sp_0 >> 63), &nu_1);
				cc = _subborrow_u64(cc, nu_2,
					(sp_2 << 1) | (sp_1 >> 63), &nu_2);
				cc = _subborrow_u64(cc, nu_3,
					(sp_3 << 1) | (sp_2 >> 63), &nu_3);
				cc = _subborrow_u64(cc, nu_4,
					(sp_4 << 1) | (sp_3 >> 63), &nu_4);
				cc = _subborrow_u64(cc, nu_5,
					(sp_5 << 1) | (sp_4 >> 63), &nu_5);
				cc = _subborrow_u64(cc, nu_6,
					(sp_6 << 1) | (sp_5 >> 63), &nu_6);
				cc = _subborrow_u64(cc, nu_7,
					(sp_7 << 1) | (sp_6 >> 63), &nu_7);
				cc = _subborrow_u64(cc, nu_8,
					(sp_8 << 1) | (sp_7 >> 63), &nu_8);
				cc = _subborrow_u64(cc, nu_9,
					(sp_9 << 1) | (sp_8 >> 63), &nu_9);
				cc = _subborrow_u64(cc, nu_10,
					(sp_10 << 1) | (sp_9 >> 63), &nu_10);
				cc = _subborrow_u64(cc, nu_11,
					(sp_11 << 1) | (sp_10 >> 63), &nu_11);
				cc = _subborrow_u64(cc, nu_12,
					(sp_12 << 1) | (sp_11 >> 63), &nu_12);
				(void)_subborrow_u64(cc, nu_13,
					(sp_13 << 1) | (sp_12 >> 63), &nu_13);

				cc = _subborrow_u64(0, sp_0, nv_0, &sp_0);
				cc = _subborrow_u64(cc, sp_1, nv_1, &sp_1);
				cc = _subborrow_u64(cc, sp_2, nv_2, &sp_2);
				cc = _subborrow_u64(cc, sp_3, nv_3, &sp_3);
				cc = _subborrow_u64(cc, sp_4, nv_4, &sp_4);
				cc = _subborrow_u64(cc, sp_5, nv_5, &sp_5);
				cc = _subborrow_u64(cc, sp_6, nv_6, &sp_6);
				cc = _subborrow_u64(cc, sp_7, nv_7, &sp_7);
				cc = _subborrow_u64(cc, sp_8, nv_8, &sp_8);
				cc = _subborrow_u64(cc, sp_9, nv_9, &sp_9);
				cc = _subborrow_u64(cc, sp_10, nv_10, &sp_10);
				cc = _subborrow_u64(cc, sp_11, nv_11, &sp_11);
				cc = _subborrow_u64(cc, sp_12, nv_12, &sp_12);
				(void)_subborrow_u64(cc, sp_13, nv_13, &sp_13);

			} else {
				unsigned rs, s1, rs1, s2, rs2;

				rs = 64 - s;
				s1 = s + 1;
				rs1 = 64 - s1;
				s2 = s << 1;
				rs2 = 64 - s2;

				cc = _subborrow_u64(0, u0_0,
					(v0_0 << s), &u0_0);
				cc = _subborrow_u64(cc, u0_1,
					(v0_1 << s) | (v0_0 >> rs), &u0_1);
				cc = _subborrow_u64(cc, u0_2,
					(v0_2 << s) | (v0_1 >> rs), &u0_2);
				(void)_subborrow_u64(cc, u0_3,
					(v0_3 << s) | (v0_2 >> rs), &u0_3);


				cc = _subborrow_u64(0, u1_0,
					(v1_0 << s), &u1_0);
				cc = _subborrow_u64(cc, u1_1,
					(v1_1 << s) | (v1_0 >> rs), &u1_1);
				cc = _subborrow_u64(cc, u1_2,
					(v1_2 << s) | (v1_1 >> rs), &u1_2);
				(void)_subborrow_u64(cc, u1_3,
					(v1_3 << s) | (v1_2 >> rs), &u1_3);

				cc = _addcarry_u64(0, nu_0,
					(nv_0 << s2), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(nv_1 << s2) | (nv_0 >> rs2), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(nv_2 << s2) | (nv_1 >> rs2), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(nv_3 << s2) | (nv_2 >> rs2), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(nv_4 << s2) | (nv_3 >> rs2), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(nv_5 << s2) | (nv_4 >> rs2), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(nv_6 << s2) | (nv_5 >> rs2), &nu_6);
				cc = _addcarry_u64(cc, nu_7,
					(nv_7 << s2) | (nv_6 >> rs2), &nu_7);
				cc = _addcarry_u64(cc, nu_8,
					(nv_8 << s2) | (nv_7 >> rs2), &nu_8);
				cc = _addcarry_u64(cc, nu_9,
					(nv_9 << s2) | (nv_8 >> rs2), &nu_9);
				cc = _addcarry_u64(cc, nu_10,
					(nv_10 << s2) | (nv_9 >> rs2), &nu_10);
				cc = _addcarry_u64(cc, nu_11,
					(nv_11 << s2) | (nv_10 >> rs2), &nu_11);
				cc = _addcarry_u64(cc, nu_12,
					(nv_12 << s2) | (nv_11 >> rs2), &nu_12);
				(void)_addcarry_u64(cc, nu_13,
					(nv_13 << s2) | (nv_12 >> rs2), &nu_13);

				cc = _subborrow_u64(0, nu_0,
					(sp_0 << s1), &nu_0);
				cc = _subborrow_u64(cc, nu_1,
					(sp_1 << s1) | (sp_0 >> rs1), &nu_1);
				cc = _subborrow_u64(cc, nu_2,
					(sp_2 << s1) | (sp_1 >> rs1), &nu_2);
				cc = _subborrow_u64(cc, nu_3,
					(sp_3 << s1) | (sp_2 >> rs1), &nu_3);
				cc = _subborrow_u64(cc, nu_4,
					(sp_4 << s1) | (sp_3 >> rs1), &nu_4);
				cc = _subborrow_u64(cc, nu_5,
					(sp_5 << s1) | (sp_4 >> rs1), &nu_5);
				cc = _subborrow_u64(cc, nu_6,
					(sp_6 << s1) | (sp_5 >> rs1), &nu_6);
				cc = _subborrow_u64(cc, nu_7,
					(sp_7 << s1) | (sp_6 >> rs1), &nu_7);
				cc = _subborrow_u64(cc, nu_8,
					(sp_8 << s1) | (sp_7 >> rs1), &nu_8);
				cc = _subborrow_u64(cc, nu_9,
					(sp_9 << s1) | (sp_8 >> rs1), &nu_9);
				cc = _subborrow_u64(cc, nu_10,
					(sp_10 << s1) | (sp_9 >> rs1), &nu_10);
				cc = _subborrow_u64(cc, nu_11,
					(sp_11 << s1) | (sp_10 >> rs1), &nu_11);
				cc = _subborrow_u64(cc, nu_12,
					(sp_12 << s1) | (sp_11 >> rs1), &nu_12);
				(void)_subborrow_u64(cc, nu_13,
					(sp_13 << s1) | (sp_12 >> rs1), &nu_13);

				cc = _subborrow_u64(0, sp_0,
					(nv_0 << s), &sp_0);
				cc = _subborrow_u64(cc, sp_1,
					(nv_1 << s) | (nv_0 >> rs), &sp_1);
				cc = _subborrow_u64(cc, sp_2,
					(nv_2 << s) | (nv_1 >> rs), &sp_2);
				cc = _subborrow_u64(cc, sp_3,
					(nv_3 << s) | (nv_2 >> rs), &sp_3);
				cc = _subborrow_u64(cc, sp_4,
					(nv_4 << s) | (nv_3 >> rs), &sp_4);
				cc = _subborrow_u64(cc, sp_5,
					(nv_5 << s) | (nv_4 >> rs), &sp_5);
				cc = _subborrow_u64(cc, sp_6,
					(nv_6 << s) | (nv_5 >> rs), &sp_6);
				cc = _subborrow_u64(cc, sp_7,
					(nv_7 << s) | (nv_6 >> rs), &sp_7);
				cc = _subborrow_u64(cc, sp_8,
					(nv_8 << s) | (nv_7 >> rs), &sp_8);
				cc = _subborrow_u64(cc, sp_9,
					(nv_9 << s) | (nv_8 >> rs), &sp_9);
				cc = _subborrow_u64(cc, sp_10,
					(nv_10 << s) | (nv_9 >> rs), &sp_10);
				cc = _subborrow_u64(cc, sp_11,
					(nv_11 << s) | (nv_10 >> rs), &sp_11);
				cc = _subborrow_u64(cc, sp_12,
					(nv_12 << s) | (nv_11 >> rs), &sp_12);
				(void)_subborrow_u64(cc, sp_13,
					(nv_13 << s) | (nv_12 >> rs), &sp_13);
			}
		} else {
			if (s == 0) {
				cc = _addcarry_u64(0, u0_0, v0_0, &u0_0);
				cc = _addcarry_u64(cc, u0_1, v0_1, &u0_1);
				cc = _addcarry_u64(cc, u0_2, v0_2, &u0_2);
				(void)_addcarry_u64(cc, u0_3, v0_3, &u0_3);

				cc = _addcarry_u64(0, u1_0, v1_0, &u1_0);
				cc = _addcarry_u64(cc, u1_1, v1_1, &u1_1);
				cc = _addcarry_u64(cc, u1_2, v1_2, &u1_2);
				(void)_addcarry_u64(cc, u1_3, v1_3, &u1_3);

				cc = _addcarry_u64(0, nu_0, nv_0, &nu_0);
				cc = _addcarry_u64(cc, nu_1, nv_1, &nu_1);
				cc = _addcarry_u64(cc, nu_2, nv_2, &nu_2);
				cc = _addcarry_u64(cc, nu_3, nv_3, &nu_3);
				cc = _addcarry_u64(cc, nu_4, nv_4, &nu_4);
				cc = _addcarry_u64(cc, nu_5, nv_5, &nu_5);
				cc = _addcarry_u64(cc, nu_6, nv_6, &nu_6);
				cc = _addcarry_u64(cc, nu_7, nv_7, &nu_7);
				cc = _addcarry_u64(cc, nu_8, nv_8, &nu_8);
				cc = _addcarry_u64(cc, nu_9, nv_9, &nu_9);
				cc = _addcarry_u64(cc, nu_10, nv_10, &nu_10);
				cc = _addcarry_u64(cc, nu_11, nv_11, &nu_11);
				cc = _addcarry_u64(cc, nu_12, nv_12, &nu_12);
				(void)_addcarry_u64(cc, nu_13, nv_13, &nu_13);

				cc = _addcarry_u64(0, nu_0,
					(sp_0 << 1), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(sp_1 << 1) | (sp_0 >> 63), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(sp_2 << 1) | (sp_1 >> 63), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(sp_3 << 1) | (sp_2 >> 63), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(sp_4 << 1) | (sp_3 >> 63), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(sp_5 << 1) | (sp_4 >> 63), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(sp_6 << 1) | (sp_5 >> 63), &nu_6);
				cc = _addcarry_u64(cc, nu_7,
					(sp_7 << 1) | (sp_6 >> 63), &nu_7);
				cc = _addcarry_u64(cc, nu_8,
					(sp_8 << 1) | (sp_7 >> 63), &nu_8);
				cc = _addcarry_u64(cc, nu_9,
					(sp_9 << 1) | (sp_8 >> 63), &nu_9);
				cc = _addcarry_u64(cc, nu_10,
					(sp_10 << 1) | (sp_9 >> 63), &nu_10);
				cc = _addcarry_u64(cc, nu_11,
					(sp_11 << 1) | (sp_10 >> 63), &nu_11);
				cc = _addcarry_u64(cc, nu_12,
					(sp_12 << 1) | (sp_11 >> 63), &nu_12);
				(void)_addcarry_u64(cc, nu_13,
					(sp_13 << 1) | (sp_12 >> 63), &nu_13);



				cc = _addcarry_u64(0, sp_0, nv_0, &sp_0);
				cc = _addcarry_u64(cc, sp_1, nv_1, &sp_1);
				cc = _addcarry_u64(cc, sp_2, nv_2, &sp_2);
				cc = _addcarry_u64(cc, sp_3, nv_3, &sp_3);
				cc = _addcarry_u64(cc, sp_4, nv_4, &sp_4);
				cc = _addcarry_u64(cc, sp_5, nv_5, &sp_5);
				cc = _addcarry_u64(cc, sp_6, nv_6, &sp_6);
				cc = _addcarry_u64(cc, sp_7, nv_7, &sp_7);
				cc = _addcarry_u64(cc, sp_8, nv_8, &sp_8);
				cc = _addcarry_u64(cc, sp_9, nv_9, &sp_9);
				cc = _addcarry_u64(cc, sp_10, nv_10, &sp_10);
				cc = _addcarry_u64(cc, sp_11, nv_11, &sp_11);
				cc = _addcarry_u64(cc, sp_12, nv_12, &sp_12);
				(void)_addcarry_u64(cc, sp_13, nv_13, &sp_13);

			} else {
				unsigned rs, s1, rs1, s2, rs2;

				rs = 64 - s;
				s1 = s + 1;
				rs1 = 64 - s1;
				s2 = s << 1;
				rs2 = 64 - s2;

				cc = _addcarry_u64(0, u0_0,
					(v0_0 << s), &u0_0);
				cc = _addcarry_u64(cc, u0_1,
					(v0_1 << s) | (v0_0 >> rs), &u0_1);
				cc = _addcarry_u64(cc, u0_2,
					(v0_2 << s) | (v0_1 >> rs), &u0_2);
				(void)_addcarry_u64(cc, u0_3,
					(v0_3 << s) | (v0_2 >> rs), &u0_3);


				cc = _addcarry_u64(0, u1_0,
					(v1_0 << s), &u1_0);
				cc = _addcarry_u64(cc, u1_1,
					(v1_1 << s) | (v1_0 >> rs), &u1_1);
				cc = _addcarry_u64(cc, u1_2,
					(v1_2 << s) | (v1_1 >> rs), &u1_2);
				(void)_addcarry_u64(cc, u1_3,
					(v1_3 << s) | (v1_2 >> rs), &u1_3);

				cc = _addcarry_u64(0, nu_0,
					(nv_0 << s2), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(nv_1 << s2) | (nv_0 >> rs2), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(nv_2 << s2) | (nv_1 >> rs2), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(nv_3 << s2) | (nv_2 >> rs2), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(nv_4 << s2) | (nv_3 >> rs2), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(nv_5 << s2) | (nv_4 >> rs2), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(nv_6 << s2) | (nv_5 >> rs2), &nu_6);
				cc = _addcarry_u64(cc, nu_7,
					(nv_7 << s2) | (nv_6 >> rs2), &nu_7);
				cc = _addcarry_u64(cc, nu_8,
					(nv_8 << s2) | (nv_7 >> rs2), &nu_8);
				cc = _addcarry_u64(cc, nu_9,
					(nv_9 << s2) | (nv_8 >> rs2), &nu_9);
				cc = _addcarry_u64(cc, nu_10,
					(nv_10 << s2) | (nv_9 >> rs2), &nu_10);
				cc = _addcarry_u64(cc, nu_11,
					(nv_11 << s2) | (nv_10 >> rs2), &nu_11);
				cc = _addcarry_u64(cc, nu_12,
					(nv_12 << s2) | (nv_11 >> rs2), &nu_12);
				(void)_addcarry_u64(cc, nu_13,
					(nv_13 << s2) | (nv_12 >> rs2), &nu_13);

				cc = _addcarry_u64(0, nu_0,
					(sp_0 << s1), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(sp_1 << s1) | (sp_0 >> rs1), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(sp_2 << s1) | (sp_1 >> rs1), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(sp_3 << s1) | (sp_2 >> rs1), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(sp_4 << s1) | (sp_3 >> rs1), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(sp_5 << s1) | (sp_4 >> rs1), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(sp_6 << s1) | (sp_5 >> rs1), &nu_6);
				cc = _addcarry_u64(cc, nu_7,
					(sp_7 << s1) | (sp_6 >> rs1), &nu_7);
				cc = _addcarry_u64(cc, nu_8,
					(sp_8 << s1) | (sp_7 >> rs1), &nu_8);
				cc = _addcarry_u64(cc, nu_9,
					(sp_9 << s1) | (sp_8 >> rs1), &nu_9);
				cc = _addcarry_u64(cc, nu_10,
					(sp_10 << s1) | (sp_9 >> rs1), &nu_10);
				cc = _addcarry_u64(cc, nu_11,
					(sp_11 << s1) | (sp_10 >> rs1), &nu_11);
				cc = _addcarry_u64(cc, nu_12,
					(sp_12 << s1) | (sp_11 >> rs1), &nu_12);
				(void)_addcarry_u64(cc, nu_13,
					(sp_13 << s1) | (sp_12 >> rs1), &nu_13);

				cc = _addcarry_u64(0, sp_0,
					(nv_0 << s), &sp_0);
				cc = _addcarry_u64(cc, sp_1,
					(nv_1 << s) | (nv_0 >> rs), &sp_1);
				cc = _addcarry_u64(cc, sp_2,
					(nv_2 << s) | (nv_1 >> rs), &sp_2);
				cc = _addcarry_u64(cc, sp_3,
					(nv_3 << s) | (nv_2 >> rs), &sp_3);
				cc = _addcarry_u64(cc, sp_4,
					(nv_4 << s) | (nv_3 >> rs), &sp_4);
				cc = _addcarry_u64(cc, sp_5,
					(nv_5 << s) | (nv_4 >> rs), &sp_5);
				cc = _addcarry_u64(cc, sp_6,
					(nv_6 << s) | (nv_5 >> rs), &sp_6);
				cc = _addcarry_u64(cc, sp_7,
					(nv_7 << s) | (nv_6 >> rs), &sp_7);
				cc = _addcarry_u64(cc, sp_8,
					(nv_8 << s) | (nv_7 >> rs), &sp_8);
				cc = _addcarry_u64(cc, sp_9,
					(nv_9 << s) | (nv_8 >> rs), &sp_9);
				cc = _addcarry_u64(cc, sp_10,
					(nv_10 << s) | (nv_9 >> rs), &sp_10);
				cc = _addcarry_u64(cc, sp_11,
					(nv_11 << s) | (nv_10 >> rs), &sp_11);
				cc = _addcarry_u64(cc, sp_12,
					(nv_12 << s) | (nv_11 >> rs), &sp_12);
				(void)_addcarry_u64(cc, sp_13,
					(nv_13 << s) | (nv_12 >> rs), &sp_13);
			}
		}
	}

	/*
	 * The part below is reached when large values have shrunk to 12
	 * words.
	 */

	for (;;) {
		unsigned s;
		unsigned long long m;
		unsigned bl_nv, bl_sp;

		/*
		 * If nu < nv, then swap(u,v) and swap(nu,nv).
		 */
		cc = _subborrow_u64(0, nu_0, nv_0, &w);
		cc = _subborrow_u64(cc, nu_1, nv_1, &w);
		cc = _subborrow_u64(cc, nu_2, nv_2, &w);
		cc = _subborrow_u64(cc, nu_3, nv_3, &w);
		cc = _subborrow_u64(cc, nu_4, nv_4, &w);
		cc = _subborrow_u64(cc, nu_5, nv_5, &w);
		cc = _subborrow_u64(cc, nu_6, nv_6, &w);
		cc = _subborrow_u64(cc, nu_7, nv_7, &w);
		cc = _subborrow_u64(cc, nu_8, nv_8, &w);
		cc = _subborrow_u64(cc, nu_9, nv_9, &w);
		cc = _subborrow_u64(cc, nu_10, nv_10, &w);
		cc = _subborrow_u64(cc, nu_11, nv_11, &w);
		// cc = _subborrow_u64(cc, nu_12, nv_12, &w);
		// cc = _subborrow_u64(cc, nu_13, nv_13, &w);
		m = -(unsigned long long)cc;


		SWAPCOND(u0_0, v0_0, m);
		SWAPCOND(u0_1, v0_1, m);
		SWAPCOND(u0_2, v0_2, m);
		SWAPCOND(u0_3, v0_3, m);
		
		SWAPCOND(u1_0, v1_0, m);
		SWAPCOND(u1_1, v1_1, m);
		SWAPCOND(u1_2, v1_2, m);
		SWAPCOND(u1_3, v1_3, m);

		SWAPCOND(nu_0, nv_0, m);
		SWAPCOND(nu_1, nv_1, m);
		SWAPCOND(nu_2, nv_2, m);
		SWAPCOND(nu_3, nv_3, m);
		SWAPCOND(nu_4, nv_4, m);
		SWAPCOND(nu_5, nv_5, m);
		SWAPCOND(nu_6, nv_6, m);
		SWAPCOND(nu_7, nv_7, m);
		SWAPCOND(nu_8, nv_8, m);
		SWAPCOND(nu_9, nv_9, m);
		SWAPCOND(nu_10, nv_10, m);
		SWAPCOND(nu_11, nv_11, m);
		// SWAPCOND(nu_12, nv_12, m);
		// SWAPCOND(nu_13, nv_13, m);

		/*
		 * u is now the largest vector; its square norm is nu.
		 * We know that:
		 *   N(u-v) = N(u) + N(v) - 2*<u,v>
		 *   N(u+v) = N(u) + N(v) + 2*<u,v>
		 * Since all squared norms are positive, and since
		 * N(u) >= N(v), it is guaranteed that |sp| <= nu.
		 * Thus, if nu fits on 10 words (including the sign bit),
		 * then we can jump to phase 3, where 'large' values use
		 * 10 words.
		 */
		if ((nu_11 | nu_10) == 0 && nu_9 <= 0x7FFFFFFFFFFFFFFFull) {
			break;
		}

#define BITLENGTH12(size, bb)   do { \
		unsigned bitlength_acc = 12; \
		int bitlength_flag = 1; \
		unsigned long long bitlength_mask = -(bb ## 11 >> 63); \
		unsigned long long bitlength_top = 0; \
		unsigned long long bitlength_word; \
		bitlength_word = bb ## 11 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 10 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 9 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 8 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 7 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 6 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 5 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 4 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 2 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 1 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 0 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		(size) = (bitlength_acc << 6) - _lzcnt_u64(bitlength_top); \
	} while (0)

		BITLENGTH12(bl_nv, nv_);
		BITLENGTH12(bl_sp, sp_);

		/*
		 * If v is small enough, return.
		 */
		if (bl_nv <= 447) {
			c0[0] = v0_0;
			c0[1] = v0_1;
			c0[2] = v0_2;
			c0[3] = v0_3;
			
			c1[0] = v1_0;
			c1[1] = v1_1;
			c1[2] = v1_2;
			c1[3] = v1_3;

			return;
		}

		/*
		 * Compute shift amount.
		 */
		s = bl_sp - bl_nv;
		s &= ~-(s >> 31);

		/*
		 * It is very rare that s > 31. We handle it with some
		 * generic code; branch prediction will soon learn that
		 * this path is normally not taken.
		 */
		if (s > 31) {
			uint64_t tu0[4], tu1[4], tv0[4], tv1[4];
			uint64_t tnu[14], tnv[14], tsp[14];

			tu0[0] = u0_0;
			tu0[1] = u0_1;
			tu0[2] = u0_2;
			tu0[3] = u0_3;

			tu1[0] = u1_0;
			tu1[1] = u1_1;
			tu1[2] = u1_2;
			tu1[3] = u1_3;

			tv0[0] = v0_0;
			tv0[1] = v0_1;
			tv0[2] = v0_2;
			tv0[3] = v0_3;

			tv1[0] = v1_0;
			tv1[1] = v1_1;
			tv1[2] = v1_2;
			tv1[3] = v1_3;

			tnu[0] = nu_0;
			tnu[1] = nu_1;
			tnu[2] = nu_2;
			tnu[3] = nu_3;
			tnu[4] = nu_4;
			tnu[5] = nu_5;
			tnu[6] = nu_6;
			tnu[7] = nu_7;
			tnu[8] = nu_8;
			tnu[9] = nu_9;
			tnu[10] = nu_10;
			tnu[11] = nu_11;
			tnu[12] = 0;
			tnu[13] = 0;

			tnv[0] = nv_0;
			tnv[1] = nv_1;
			tnv[2] = nv_2;
			tnv[3] = nv_3;
			tnv[4] = nv_4;
			tnv[5] = nv_5;
			tnv[6] = nv_6;
			tnv[7] = nv_7;
			tnv[8] = nv_8;
			tnv[9] = nv_9;
			tnv[10] = nv_10;
			tnv[11] = nv_11;
			tnv[12] = 0;
			tnv[13] = 0;

			tsp[0] = sp_0;
			tsp[1] = sp_1;
			tsp[2] = sp_2;
			tsp[3] = sp_3;
			tsp[4] = sp_4;
			tsp[5] = sp_5;
			tsp[6] = sp_6;
			tsp[7] = sp_7;
			tsp[8] = sp_8;
			tsp[9] = sp_9;
			tsp[10] = sp_10;
			tsp[11] = sp_11;
			tsp[12] = -(sp_11 >> 63);
			tsp[13] = tsp[12];

			if ((sp_13 >> 63) == 0) {
				sub_lshift_small(tu0, tv0, s);
				sub_lshift_small(tu1, tv1, s);
				add_lshift_large(tnu, tnv, 2 * s);
				sub_lshift_large(tnu, tsp, s + 1);
				sub_lshift_large(tsp, tnv, s);
			} else {
				add_lshift_small(tu0, tv0, s);
				add_lshift_small(tu1, tv1, s);
				add_lshift_large(tnu, tnv, 2 * s);
				add_lshift_large(tnu, tsp, s + 1);
				add_lshift_large(tsp, tnv, s);
			}


			u0_0 = tu0[0];
			u0_1 = tu0[1];
			u0_2 = tu0[2];
			u0_3 = tu0[3];

			u1_0 = tu1[0];
			u1_1 = tu1[1];
			u1_2 = tu1[2];
			u1_3 = tu1[3];

			v0_0 = tv0[0];
			v0_1 = tv0[1];
			v0_2 = tv0[2];
			v0_3 = tv0[3];

			v1_0 = tv1[0];
			v1_1 = tv1[1];
			v1_2 = tv1[2];
			v1_3 = tv1[3];

			nu_0 = tnu[0];
			nu_1 = tnu[1];
			nu_2 = tnu[2];
			nu_3 = tnu[3];
			nu_4 = tnu[4];
			nu_5 = tnu[5];
			nu_6 = tnu[6];
			nu_7 = tnu[7];
			nu_8 = tnu[8];
			nu_9 = tnu[9];
			nu_10 = tnu[10];
			nu_11 = tnu[11];
			// nu_12 = tnu[12];
			// nu_13 = tnu[13];

			

			nv_0 = tnv[0];
			nv_1 = tnv[1];
			nv_2 = tnv[2];
			nv_3 = tnv[3];
			nv_4 = tnv[4];
			nv_5 = tnv[5];
			nv_6 = tnv[6];
			nv_7 = tnv[7];
			nv_8 = tnv[8];
			nv_9 = tnv[9];
			nv_10 = tnv[10];
			nv_11 = tnv[11];
			// nv_12 = tnv[12];
			// nv_13 = tnv[13];

			sp_0 = tsp[0];
			sp_1 = tsp[1];
			sp_2 = tsp[2];
			sp_3 = tsp[3];
			sp_4 = tsp[4];
			sp_5 = tsp[5];
			sp_6 = tsp[6];
			sp_7 = tsp[7];
			sp_8 = tsp[8];
			sp_9 = tsp[9];
			sp_10 = tsp[10];
			sp_11 = tsp[11];
			// sp_12 = tsp[12];
			// sp_13 = tsp[13];
			continue;
		}

		if ((sp_11 >> 63) == 0) {
			if (s == 0) {
				cc = _subborrow_u64(0, u0_0, v0_0, &u0_0);
				cc = _subborrow_u64(cc, u0_1, v0_1, &u0_1);
				cc = _subborrow_u64(cc, u0_2, v0_2, &u0_2);
				(void)_subborrow_u64(cc, u0_3, v0_3, &u0_3);

				cc = _subborrow_u64(0, u1_0, v1_0, &u1_0);
				cc = _subborrow_u64(cc, u1_1, v1_1, &u1_1);
				cc = _subborrow_u64(cc, u1_2, v1_2, &u1_2);
				(void)_subborrow_u64(cc, u1_3, v1_3, &u1_3);

				cc = _addcarry_u64(0, nu_0, nv_0, &nu_0);
				cc = _addcarry_u64(cc, nu_1, nv_1, &nu_1);
				cc = _addcarry_u64(cc, nu_2, nv_2, &nu_2);
				cc = _addcarry_u64(cc, nu_3, nv_3, &nu_3);
				cc = _addcarry_u64(cc, nu_4, nv_4, &nu_4);
				cc = _addcarry_u64(cc, nu_5, nv_5, &nu_5);
				cc = _addcarry_u64(cc, nu_6, nv_6, &nu_6);
				cc = _addcarry_u64(cc, nu_7, nv_7, &nu_7);
				cc = _addcarry_u64(cc, nu_8, nv_8, &nu_8);
				cc = _addcarry_u64(cc, nu_9, nv_9, &nu_9);
				cc = _addcarry_u64(cc, nu_10, nv_10, &nu_10);
				(void)_addcarry_u64(cc, nu_11, nv_11, &nu_11);
				// cc = _addcarry_u64(cc, nu_12, nv_12, &nu_12);
				// (void)_addcarry_u64(cc, nu_13, nv_13, &nu_13);

				cc = _subborrow_u64(0, nu_0,
					(sp_0 << 1), &nu_0);
				cc = _subborrow_u64(cc, nu_1,
					(sp_1 << 1) | (sp_0 >> 63), &nu_1);
				cc = _subborrow_u64(cc, nu_2,
					(sp_2 << 1) | (sp_1 >> 63), &nu_2);
				cc = _subborrow_u64(cc, nu_3,
					(sp_3 << 1) | (sp_2 >> 63), &nu_3);
				cc = _subborrow_u64(cc, nu_4,
					(sp_4 << 1) | (sp_3 >> 63), &nu_4);
				cc = _subborrow_u64(cc, nu_5,
					(sp_5 << 1) | (sp_4 >> 63), &nu_5);
				cc = _subborrow_u64(cc, nu_6,
					(sp_6 << 1) | (sp_5 >> 63), &nu_6);
				cc = _subborrow_u64(cc, nu_7,
					(sp_7 << 1) | (sp_6 >> 63), &nu_7);
				cc = _subborrow_u64(cc, nu_8,
					(sp_8 << 1) | (sp_7 >> 63), &nu_8);
				cc = _subborrow_u64(cc, nu_9,
					(sp_9 << 1) | (sp_8 >> 63), &nu_9);
				cc = _subborrow_u64(cc, nu_10,
					(sp_10 << 1) | (sp_9 >> 63), &nu_10);
				(void)_subborrow_u64(cc, nu_11,
					(sp_11 << 1) | (sp_10 >> 63), &nu_11);
				// cc = _subborrow_u64(cc, nu_12,
				// 	(sp_12 << 1) | (sp_11 >> 63), &nu_12);
				// (void)_subborrow_u64(cc, nu_13,
				// 	(sp_13 << 1) | (sp_12 >> 63), &nu_13);

				cc = _subborrow_u64(0, sp_0, nv_0, &sp_0);
				cc = _subborrow_u64(cc, sp_1, nv_1, &sp_1);
				cc = _subborrow_u64(cc, sp_2, nv_2, &sp_2);
				cc = _subborrow_u64(cc, sp_3, nv_3, &sp_3);
				cc = _subborrow_u64(cc, sp_4, nv_4, &sp_4);
				cc = _subborrow_u64(cc, sp_5, nv_5, &sp_5);
				cc = _subborrow_u64(cc, sp_6, nv_6, &sp_6);
				cc = _subborrow_u64(cc, sp_7, nv_7, &sp_7);
				cc = _subborrow_u64(cc, sp_8, nv_8, &sp_8);
				cc = _subborrow_u64(cc, sp_9, nv_9, &sp_9);
				cc = _subborrow_u64(cc, sp_10, nv_10, &sp_10);
				(void)_subborrow_u64(cc, sp_11, nv_11, &sp_11);
				// cc = _subborrow_u64(cc, sp_12, nv_12, &sp_12);
				// (void)_subborrow_u64(cc, sp_13, nv_13, &sp_13);

			} else {
				unsigned rs, s1, rs1, s2, rs2;

				rs = 64 - s;
				s1 = s + 1;
				rs1 = 64 - s1;
				s2 = s << 1;
				rs2 = 64 - s2;

				cc = _subborrow_u64(0, u0_0,
					(v0_0 << s), &u0_0);
				cc = _subborrow_u64(cc, u0_1,
					(v0_1 << s) | (v0_0 >> rs), &u0_1);
				cc = _subborrow_u64(cc, u0_2,
					(v0_2 << s) | (v0_1 >> rs), &u0_2);
				(void)_subborrow_u64(cc, u0_3,
					(v0_3 << s) | (v0_2 >> rs), &u0_3);


				cc = _subborrow_u64(0, u1_0,
					(v1_0 << s), &u1_0);
				cc = _subborrow_u64(cc, u1_1,
					(v1_1 << s) | (v1_0 >> rs), &u1_1);
				cc = _subborrow_u64(cc, u1_2,
					(v1_2 << s) | (v1_1 >> rs), &u1_2);
				(void)_subborrow_u64(cc, u1_3,
					(v1_3 << s) | (v1_2 >> rs), &u1_3);

				cc = _addcarry_u64(0, nu_0,
					(nv_0 << s2), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(nv_1 << s2) | (nv_0 >> rs2), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(nv_2 << s2) | (nv_1 >> rs2), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(nv_3 << s2) | (nv_2 >> rs2), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(nv_4 << s2) | (nv_3 >> rs2), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(nv_5 << s2) | (nv_4 >> rs2), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(nv_6 << s2) | (nv_5 >> rs2), &nu_6);
				cc = _addcarry_u64(cc, nu_7,
					(nv_7 << s2) | (nv_6 >> rs2), &nu_7);
				cc = _addcarry_u64(cc, nu_8,
					(nv_8 << s2) | (nv_7 >> rs2), &nu_8);
				cc = _addcarry_u64(cc, nu_9,
					(nv_9 << s2) | (nv_8 >> rs2), &nu_9);
				cc = _addcarry_u64(cc, nu_10,
					(nv_10 << s2) | (nv_9 >> rs2), &nu_10);
				(void)_addcarry_u64(cc, nu_11,
					(nv_11 << s2) | (nv_10 >> rs2), &nu_11);
				// cc = _addcarry_u64(cc, nu_12,
				// 	(nv_12 << s2) | (nv_11 >> rs2), &nu_12);
				// (void)_addcarry_u64(cc, nu_13,
				// 	(nv_13 << s2) | (nv_12 >> rs2), &nu_13);

				cc = _subborrow_u64(0, nu_0,
					(sp_0 << s1), &nu_0);
				cc = _subborrow_u64(cc, nu_1,
					(sp_1 << s1) | (sp_0 >> rs1), &nu_1);
				cc = _subborrow_u64(cc, nu_2,
					(sp_2 << s1) | (sp_1 >> rs1), &nu_2);
				cc = _subborrow_u64(cc, nu_3,
					(sp_3 << s1) | (sp_2 >> rs1), &nu_3);
				cc = _subborrow_u64(cc, nu_4,
					(sp_4 << s1) | (sp_3 >> rs1), &nu_4);
				cc = _subborrow_u64(cc, nu_5,
					(sp_5 << s1) | (sp_4 >> rs1), &nu_5);
				cc = _subborrow_u64(cc, nu_6,
					(sp_6 << s1) | (sp_5 >> rs1), &nu_6);
				cc = _subborrow_u64(cc, nu_7,
					(sp_7 << s1) | (sp_6 >> rs1), &nu_7);
				cc = _subborrow_u64(cc, nu_8,
					(sp_8 << s1) | (sp_7 >> rs1), &nu_8);
				cc = _subborrow_u64(cc, nu_9,
					(sp_9 << s1) | (sp_8 >> rs1), &nu_9);
				cc = _subborrow_u64(cc, nu_10,
					(sp_10 << s1) | (sp_9 >> rs1), &nu_10);
				(void)_subborrow_u64(cc, nu_11,
					(sp_11 << s1) | (sp_10 >> rs1), &nu_11);
				// cc = _subborrow_u64(cc, nu_12,
				// 	(sp_12 << s1) | (sp_11 >> rs1), &nu_12);
				// (void)_subborrow_u64(cc, nu_13,
				// 	(sp_13 << s1) | (sp_12 >> rs1), &nu_13);

				cc = _subborrow_u64(0, sp_0,
					(nv_0 << s), &sp_0);
				cc = _subborrow_u64(cc, sp_1,
					(nv_1 << s) | (nv_0 >> rs), &sp_1);
				cc = _subborrow_u64(cc, sp_2,
					(nv_2 << s) | (nv_1 >> rs), &sp_2);
				cc = _subborrow_u64(cc, sp_3,
					(nv_3 << s) | (nv_2 >> rs), &sp_3);
				cc = _subborrow_u64(cc, sp_4,
					(nv_4 << s) | (nv_3 >> rs), &sp_4);
				cc = _subborrow_u64(cc, sp_5,
					(nv_5 << s) | (nv_4 >> rs), &sp_5);
				cc = _subborrow_u64(cc, sp_6,
					(nv_6 << s) | (nv_5 >> rs), &sp_6);
				cc = _subborrow_u64(cc, sp_7,
					(nv_7 << s) | (nv_6 >> rs), &sp_7);
				cc = _subborrow_u64(cc, sp_8,
					(nv_8 << s) | (nv_7 >> rs), &sp_8);
				cc = _subborrow_u64(cc, sp_9,
					(nv_9 << s) | (nv_8 >> rs), &sp_9);
				cc = _subborrow_u64(cc, sp_10,
					(nv_10 << s) | (nv_9 >> rs), &sp_10);
				(void)_subborrow_u64(cc, sp_11,
					(nv_11 << s) | (nv_10 >> rs), &sp_11);
				// cc = _subborrow_u64(cc, sp_12,
				// 	(nv_12 << s) | (nv_11 >> rs), &sp_12);
				// (void)_subborrow_u64(cc, sp_13,
				// 	(nv_13 << s) | (nv_12 >> rs), &sp_13);
			}
		} else {
			if (s == 0) {
				cc = _addcarry_u64(0, u0_0, v0_0, &u0_0);
				cc = _addcarry_u64(cc, u0_1, v0_1, &u0_1);
				cc = _addcarry_u64(cc, u0_2, v0_2, &u0_2);
				(void)_addcarry_u64(cc, u0_3, v0_3, &u0_3);

				cc = _addcarry_u64(0, u1_0, v1_0, &u1_0);
				cc = _addcarry_u64(cc, u1_1, v1_1, &u1_1);
				cc = _addcarry_u64(cc, u1_2, v1_2, &u1_2);
				(void)_addcarry_u64(cc, u1_3, v1_3, &u1_3);

				cc = _addcarry_u64(0, nu_0, nv_0, &nu_0);
				cc = _addcarry_u64(cc, nu_1, nv_1, &nu_1);
				cc = _addcarry_u64(cc, nu_2, nv_2, &nu_2);
				cc = _addcarry_u64(cc, nu_3, nv_3, &nu_3);
				cc = _addcarry_u64(cc, nu_4, nv_4, &nu_4);
				cc = _addcarry_u64(cc, nu_5, nv_5, &nu_5);
				cc = _addcarry_u64(cc, nu_6, nv_6, &nu_6);
				cc = _addcarry_u64(cc, nu_7, nv_7, &nu_7);
				cc = _addcarry_u64(cc, nu_8, nv_8, &nu_8);
				cc = _addcarry_u64(cc, nu_9, nv_9, &nu_9);
				cc = _addcarry_u64(cc, nu_10, nv_10, &nu_10);
				(void)_addcarry_u64(cc, nu_11, nv_11, &nu_11);
				// cc = _addcarry_u64(cc, nu_12, nv_12, &nu_12);
				// (void)_addcarry_u64(cc, nu_13, nv_13, &nu_13);

				cc = _addcarry_u64(0, nu_0,
					(sp_0 << 1), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(sp_1 << 1) | (sp_0 >> 63), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(sp_2 << 1) | (sp_1 >> 63), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(sp_3 << 1) | (sp_2 >> 63), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(sp_4 << 1) | (sp_3 >> 63), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(sp_5 << 1) | (sp_4 >> 63), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(sp_6 << 1) | (sp_5 >> 63), &nu_6);
				cc = _addcarry_u64(cc, nu_7,
					(sp_7 << 1) | (sp_6 >> 63), &nu_7);
				cc = _addcarry_u64(cc, nu_8,
					(sp_8 << 1) | (sp_7 >> 63), &nu_8);
				cc = _addcarry_u64(cc, nu_9,
					(sp_9 << 1) | (sp_8 >> 63), &nu_9);
				cc = _addcarry_u64(cc, nu_10,
					(sp_10 << 1) | (sp_9 >> 63), &nu_10);
				(void)_addcarry_u64(cc, nu_11,
					(sp_11 << 1) | (sp_10 >> 63), &nu_11);
				// cc = _addcarry_u64(cc, nu_12,
				// 	(sp_12 << 1) | (sp_11 >> 63), &nu_12);
				// (void)_addcarry_u64(cc, nu_13,
				// 	(sp_13 << 1) | (sp_12 >> 63), &nu_13);


				cc = _addcarry_u64(0, sp_0, nv_0, &sp_0);
				cc = _addcarry_u64(cc, sp_1, nv_1, &sp_1);
				cc = _addcarry_u64(cc, sp_2, nv_2, &sp_2);
				cc = _addcarry_u64(cc, sp_3, nv_3, &sp_3);
				cc = _addcarry_u64(cc, sp_4, nv_4, &sp_4);
				cc = _addcarry_u64(cc, sp_5, nv_5, &sp_5);
				cc = _addcarry_u64(cc, sp_6, nv_6, &sp_6);
				cc = _addcarry_u64(cc, sp_7, nv_7, &sp_7);
				cc = _addcarry_u64(cc, sp_8, nv_8, &sp_8);
				cc = _addcarry_u64(cc, sp_9, nv_9, &sp_9);
				cc = _addcarry_u64(cc, sp_10, nv_10, &sp_10);
				(void)_addcarry_u64(cc, sp_11, nv_11, &sp_11);
				// cc = _addcarry_u64(cc, sp_12, nv_12, &sp_12);
				// (void)_addcarry_u64(cc, sp_13, nv_13, &sp_13);

			} else {
				unsigned rs, s1, rs1, s2, rs2;

				rs = 64 - s;
				s1 = s + 1;
				rs1 = 64 - s1;
				s2 = s << 1;
				rs2 = 64 - s2;

				cc = _addcarry_u64(0, u0_0,
					(v0_0 << s), &u0_0);
				cc = _addcarry_u64(cc, u0_1,
					(v0_1 << s) | (v0_0 >> rs), &u0_1);
				cc = _addcarry_u64(cc, u0_2,
					(v0_2 << s) | (v0_1 >> rs), &u0_2);
				(void)_addcarry_u64(cc, u0_3,
					(v0_3 << s) | (v0_2 >> rs), &u0_3);


				cc = _addcarry_u64(0, u1_0,
					(v1_0 << s), &u1_0);
				cc = _addcarry_u64(cc, u1_1,
					(v1_1 << s) | (v1_0 >> rs), &u1_1);
				cc = _addcarry_u64(cc, u1_2,
					(v1_2 << s) | (v1_1 >> rs), &u1_2);
				(void)_addcarry_u64(cc, u1_3,
					(v1_3 << s) | (v1_2 >> rs), &u1_3);

				cc = _addcarry_u64(0, nu_0,
					(nv_0 << s2), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(nv_1 << s2) | (nv_0 >> rs2), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(nv_2 << s2) | (nv_1 >> rs2), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(nv_3 << s2) | (nv_2 >> rs2), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(nv_4 << s2) | (nv_3 >> rs2), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(nv_5 << s2) | (nv_4 >> rs2), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(nv_6 << s2) | (nv_5 >> rs2), &nu_6);
				cc = _addcarry_u64(cc, nu_7,
					(nv_7 << s2) | (nv_6 >> rs2), &nu_7);
				cc = _addcarry_u64(cc, nu_8,
					(nv_8 << s2) | (nv_7 >> rs2), &nu_8);
				cc = _addcarry_u64(cc, nu_9,
					(nv_9 << s2) | (nv_8 >> rs2), &nu_9);
				cc = _addcarry_u64(cc, nu_10,
					(nv_10 << s2) | (nv_9 >> rs2), &nu_10);
				(void)_addcarry_u64(cc, nu_11,
					(nv_11 << s2) | (nv_10 >> rs2), &nu_11);
				// cc = _addcarry_u64(cc, nu_12,
				// 	(nv_12 << s2) | (nv_11 >> rs2), &nu_12);
				// (void)_addcarry_u64(cc, nu_13,
				// 	(nv_13 << s2) | (nv_12 >> rs2), &nu_13);

				cc = _addcarry_u64(0, nu_0,
					(sp_0 << s1), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(sp_1 << s1) | (sp_0 >> rs1), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(sp_2 << s1) | (sp_1 >> rs1), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(sp_3 << s1) | (sp_2 >> rs1), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(sp_4 << s1) | (sp_3 >> rs1), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(sp_5 << s1) | (sp_4 >> rs1), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(sp_6 << s1) | (sp_5 >> rs1), &nu_6);
				cc = _addcarry_u64(cc, nu_7,
					(sp_7 << s1) | (sp_6 >> rs1), &nu_7);
				cc = _addcarry_u64(cc, nu_8,
					(sp_8 << s1) | (sp_7 >> rs1), &nu_8);
				cc = _addcarry_u64(cc, nu_9,
					(sp_9 << s1) | (sp_8 >> rs1), &nu_9);
				cc = _addcarry_u64(cc, nu_10,
					(sp_10 << s1) | (sp_9 >> rs1), &nu_10);
				(void)_addcarry_u64(cc, nu_11,
					(sp_11 << s1) | (sp_10 >> rs1), &nu_11);
				// cc = _addcarry_u64(cc, nu_12,
				// 	(sp_12 << s1) | (sp_11 >> rs1), &nu_12);
				// (void)_addcarry_u64(cc, nu_13,
				// 	(sp_13 << s1) | (sp_12 >> rs1), &nu_13);

				cc = _addcarry_u64(0, sp_0,
					(nv_0 << s), &sp_0);
				cc = _addcarry_u64(cc, sp_1,
					(nv_1 << s) | (nv_0 >> rs), &sp_1);
				cc = _addcarry_u64(cc, sp_2,
					(nv_2 << s) | (nv_1 >> rs), &sp_2);
				cc = _addcarry_u64(cc, sp_3,
					(nv_3 << s) | (nv_2 >> rs), &sp_3);
				cc = _addcarry_u64(cc, sp_4,
					(nv_4 << s) | (nv_3 >> rs), &sp_4);
				cc = _addcarry_u64(cc, sp_5,
					(nv_5 << s) | (nv_4 >> rs), &sp_5);
				cc = _addcarry_u64(cc, sp_6,
					(nv_6 << s) | (nv_5 >> rs), &sp_6);
				cc = _addcarry_u64(cc, sp_7,
					(nv_7 << s) | (nv_6 >> rs), &sp_7);
				cc = _addcarry_u64(cc, sp_8,
					(nv_8 << s) | (nv_7 >> rs), &sp_8);
				cc = _addcarry_u64(cc, sp_9,
					(nv_9 << s) | (nv_8 >> rs), &sp_9);
				cc = _addcarry_u64(cc, sp_10,
					(nv_10 << s) | (nv_9 >> rs), &sp_10);
				(void)_addcarry_u64(cc, sp_11,
					(nv_11 << s) | (nv_10 >> rs), &sp_11);
				// cc = _addcarry_u64(cc, sp_12,
				// 	(nv_12 << s) | (nv_11 >> rs), &sp_12);
				// (void)_addcarry_u64(cc, sp_13,
				// 	(nv_13 << s) | (nv_12 >> rs), &sp_13);
			}
		}
	}



	/*
	 * The part below is reached when large values have shrunk to 10
	 * words.
	 */

	for (;;) {
		unsigned s;
		unsigned long long m;
		unsigned bl_nv, bl_sp;

		/*
		 * If nu < nv, then swap(u,v) and swap(nu,nv).
		 */
		cc = _subborrow_u64(0, nu_0, nv_0, &w);
		cc = _subborrow_u64(cc, nu_1, nv_1, &w);
		cc = _subborrow_u64(cc, nu_2, nv_2, &w);
		cc = _subborrow_u64(cc, nu_3, nv_3, &w);
		cc = _subborrow_u64(cc, nu_4, nv_4, &w);
		cc = _subborrow_u64(cc, nu_5, nv_5, &w);
		cc = _subborrow_u64(cc, nu_6, nv_6, &w);
		cc = _subborrow_u64(cc, nu_7, nv_7, &w);
		cc = _subborrow_u64(cc, nu_8, nv_8, &w);
		cc = _subborrow_u64(cc, nu_9, nv_9, &w);
		// cc = _subborrow_u64(cc, nu_10, nv_10, &w);
		// cc = _subborrow_u64(cc, nu_11, nv_11, &w);
		// cc = _subborrow_u64(cc, nu_12, nv_12, &w);
		// cc = _subborrow_u64(cc, nu_13, nv_13, &w);
		m = -(unsigned long long)cc;


		SWAPCOND(u0_0, v0_0, m);
		SWAPCOND(u0_1, v0_1, m);
		SWAPCOND(u0_2, v0_2, m);
		SWAPCOND(u0_3, v0_3, m);
		
		SWAPCOND(u1_0, v1_0, m);
		SWAPCOND(u1_1, v1_1, m);
		SWAPCOND(u1_2, v1_2, m);
		SWAPCOND(u1_3, v1_3, m);

		SWAPCOND(nu_0, nv_0, m);
		SWAPCOND(nu_1, nv_1, m);
		SWAPCOND(nu_2, nv_2, m);
		SWAPCOND(nu_3, nv_3, m);
		SWAPCOND(nu_4, nv_4, m);
		SWAPCOND(nu_5, nv_5, m);
		SWAPCOND(nu_6, nv_6, m);
		SWAPCOND(nu_7, nv_7, m);
		SWAPCOND(nu_8, nv_8, m);
		SWAPCOND(nu_9, nv_9, m);
		// SWAPCOND(nu_10, nv_10, m);
		// SWAPCOND(nu_11, nv_11, m);
		// SWAPCOND(nu_12, nv_12, m);
		// SWAPCOND(nu_13, nv_13, m);

		/*
		 * u is now the largest vector; its square norm is nu.
		 * We know that:
		 *   N(u-v) = N(u) + N(v) - 2*<u,v>
		 *   N(u+v) = N(u) + N(v) + 2*<u,v>
		 * Since all squared norms are positive, and since
		 * N(u) >= N(v), it is guaranteed that |sp| <= nu.
		 * Thus, if nu fits on 8 words (including the sign bit),
		 * then we can jump to phase 4, where 'large' values use
		 * 8 words.
		 */
		if ((nu_9 | nu_8) == 0 && nu_7 <= 0x7FFFFFFFFFFFFFFFull) {
			break;
		}

#define BITLENGTH10(size, bb)   do { \
		unsigned bitlength_acc = 10; \
		int bitlength_flag = 1; \
		unsigned long long bitlength_mask = -(bb ## 9 >> 63); \
		unsigned long long bitlength_top = 0; \
		unsigned long long bitlength_word; \
		bitlength_word = bb ## 9 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 8 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 7 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 6 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 5 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 4 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 3 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 2 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 1 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 0 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		(size) = (bitlength_acc << 6) - _lzcnt_u64(bitlength_top); \
	} while (0)

		BITLENGTH10(bl_nv, nv_);
		BITLENGTH10(bl_sp, sp_);

		/*
		 * If v is small enough, return.
		 */
		if (bl_nv <= 447) {
			c0[0] = v0_0;
			c0[1] = v0_1;
			c0[2] = v0_2;
			c0[3] = v0_3;
			
			c1[0] = v1_0;
			c1[1] = v1_1;
			c1[2] = v1_2;
			c1[3] = v1_3;

			return;
		}

		/*
		 * Compute shift amount.
		 */
		s = bl_sp - bl_nv;
		s &= ~-(s >> 31);

		/*
		 * It is very rare that s > 31. We handle it with some
		 * generic code; branch prediction will soon learn that
		 * this path is normally not taken.
		 */
		if (s > 31) {
			uint64_t tu0[4], tu1[4], tv0[4], tv1[4];
			uint64_t tnu[14], tnv[14], tsp[14];

			tu0[0] = u0_0;
			tu0[1] = u0_1;
			tu0[2] = u0_2;
			tu0[3] = u0_3;

			tu1[0] = u1_0;
			tu1[1] = u1_1;
			tu1[2] = u1_2;
			tu1[3] = u1_3;

			tv0[0] = v0_0;
			tv0[1] = v0_1;
			tv0[2] = v0_2;
			tv0[3] = v0_3;

			tv1[0] = v1_0;
			tv1[1] = v1_1;
			tv1[2] = v1_2;
			tv1[3] = v1_3;

			tnu[0] = nu_0;
			tnu[1] = nu_1;
			tnu[2] = nu_2;
			tnu[3] = nu_3;
			tnu[4] = nu_4;
			tnu[5] = nu_5;
			tnu[6] = nu_6;
			tnu[7] = nu_7;
			tnu[8] = nu_8;
			tnu[9] = nu_9;
			tnu[10] = 0;
			tnu[11] = 0;
			tnu[12] = 0;
			tnu[13] = 0;

			tnv[0] = nv_0;
			tnv[1] = nv_1;
			tnv[2] = nv_2;
			tnv[3] = nv_3;
			tnv[4] = nv_4;
			tnv[5] = nv_5;
			tnv[6] = nv_6;
			tnv[7] = nv_7;
			tnv[8] = nv_8;
			tnv[9] = nv_9;
			tnv[10] = 0;
			tnv[11] = 0;
			tnv[12] = 0;
			tnv[13] = 0;

			tsp[0] = sp_0;
			tsp[1] = sp_1;
			tsp[2] = sp_2;
			tsp[3] = sp_3;
			tsp[4] = sp_4;
			tsp[5] = sp_5;
			tsp[6] = sp_6;
			tsp[7] = sp_7;
			tsp[8] = sp_8;
			tsp[9] = sp_9;
			tsp[10] = -(sp_9 >> 63);
			tsp[11] = tsp[10];
			tsp[12] = tsp[10];
			tsp[13] = tsp[10];

			if ((sp_13 >> 63) == 0) {
				sub_lshift_small(tu0, tv0, s);
				sub_lshift_small(tu1, tv1, s);
				add_lshift_large(tnu, tnv, 2 * s);
				sub_lshift_large(tnu, tsp, s + 1);
				sub_lshift_large(tsp, tnv, s);
			} else {
				add_lshift_small(tu0, tv0, s);
				add_lshift_small(tu1, tv1, s);
				add_lshift_large(tnu, tnv, 2 * s);
				add_lshift_large(tnu, tsp, s + 1);
				add_lshift_large(tsp, tnv, s);
			}


			u0_0 = tu0[0];
			u0_1 = tu0[1];
			u0_2 = tu0[2];
			u0_3 = tu0[3];

			u1_0 = tu1[0];
			u1_1 = tu1[1];
			u1_2 = tu1[2];
			u1_3 = tu1[3];

			v0_0 = tv0[0];
			v0_1 = tv0[1];
			v0_2 = tv0[2];
			v0_3 = tv0[3];

			v1_0 = tv1[0];
			v1_1 = tv1[1];
			v1_2 = tv1[2];
			v1_3 = tv1[3];

			nu_0 = tnu[0];
			nu_1 = tnu[1];
			nu_2 = tnu[2];
			nu_3 = tnu[3];
			nu_4 = tnu[4];
			nu_5 = tnu[5];
			nu_6 = tnu[6];
			nu_7 = tnu[7];
			nu_8 = tnu[8];
			nu_9 = tnu[9];
			// nu_10 = tnu[10];
			// nu_11 = tnu[11];
			// nu_12 = tnu[12];
			// nu_13 = tnu[13];

			

			nv_0 = tnv[0];
			nv_1 = tnv[1];
			nv_2 = tnv[2];
			nv_3 = tnv[3];
			nv_4 = tnv[4];
			nv_5 = tnv[5];
			nv_6 = tnv[6];
			nv_7 = tnv[7];
			nv_8 = tnv[8];
			nv_9 = tnv[9];
			// nv_10 = tnv[10];
			// nv_11 = tnv[11];
			// nv_12 = tnv[12];
			// nv_13 = tnv[13];

			sp_0 = tsp[0];
			sp_1 = tsp[1];
			sp_2 = tsp[2];
			sp_3 = tsp[3];
			sp_4 = tsp[4];
			sp_5 = tsp[5];
			sp_6 = tsp[6];
			sp_7 = tsp[7];
			sp_8 = tsp[8];
			sp_9 = tsp[9];
			// sp_10 = tsp[10];
			// sp_11 = tsp[11];
			// sp_12 = tsp[12];
			// sp_13 = tsp[13];
			continue;
		}

		if ((sp_9 >> 63) == 0) {
			if (s == 0) {
				cc = _subborrow_u64(0, u0_0, v0_0, &u0_0);
				cc = _subborrow_u64(cc, u0_1, v0_1, &u0_1);
				cc = _subborrow_u64(cc, u0_2, v0_2, &u0_2);
				(void)_subborrow_u64(cc, u0_3, v0_3, &u0_3);

				cc = _subborrow_u64(0, u1_0, v1_0, &u1_0);
				cc = _subborrow_u64(cc, u1_1, v1_1, &u1_1);
				cc = _subborrow_u64(cc, u1_2, v1_2, &u1_2);
				(void)_subborrow_u64(cc, u1_3, v1_3, &u1_3);

				cc = _addcarry_u64(0, nu_0, nv_0, &nu_0);
				cc = _addcarry_u64(cc, nu_1, nv_1, &nu_1);
				cc = _addcarry_u64(cc, nu_2, nv_2, &nu_2);
				cc = _addcarry_u64(cc, nu_3, nv_3, &nu_3);
				cc = _addcarry_u64(cc, nu_4, nv_4, &nu_4);
				cc = _addcarry_u64(cc, nu_5, nv_5, &nu_5);
				cc = _addcarry_u64(cc, nu_6, nv_6, &nu_6);
				cc = _addcarry_u64(cc, nu_7, nv_7, &nu_7);
				cc = _addcarry_u64(cc, nu_8, nv_8, &nu_8);
				(void)_addcarry_u64(cc, nu_9, nv_9, &nu_9);
				// cc = _addcarry_u64(cc, nu_10, nv_10, &nu_10);
				// (void)_addcarry_u64(cc, nu_11, nv_11, &nu_11);
				// cc = _addcarry_u64(cc, nu_12, nv_12, &nu_12);
				// (void)_addcarry_u64(cc, nu_13, nv_13, &nu_13);

				cc = _subborrow_u64(0, nu_0,
					(sp_0 << 1), &nu_0);
				cc = _subborrow_u64(cc, nu_1,
					(sp_1 << 1) | (sp_0 >> 63), &nu_1);
				cc = _subborrow_u64(cc, nu_2,
					(sp_2 << 1) | (sp_1 >> 63), &nu_2);
				cc = _subborrow_u64(cc, nu_3,
					(sp_3 << 1) | (sp_2 >> 63), &nu_3);
				cc = _subborrow_u64(cc, nu_4,
					(sp_4 << 1) | (sp_3 >> 63), &nu_4);
				cc = _subborrow_u64(cc, nu_5,
					(sp_5 << 1) | (sp_4 >> 63), &nu_5);
				cc = _subborrow_u64(cc, nu_6,
					(sp_6 << 1) | (sp_5 >> 63), &nu_6);
				cc = _subborrow_u64(cc, nu_7,
					(sp_7 << 1) | (sp_6 >> 63), &nu_7);
				cc = _subborrow_u64(cc, nu_8,
					(sp_8 << 1) | (sp_7 >> 63), &nu_8);
				(void)_subborrow_u64(cc, nu_9,
					(sp_9 << 1) | (sp_8 >> 63), &nu_9);
				// cc = _subborrow_u64(cc, nu_10,
				// 	(sp_10 << 1) | (sp_9 >> 63), &nu_10);
				// (void)_subborrow_u64(cc, nu_11,
				// 	(sp_11 << 1) | (sp_10 >> 63), &nu_11);
				// cc = _subborrow_u64(cc, nu_12,
				// 	(sp_12 << 1) | (sp_11 >> 63), &nu_12);
				// (void)_subborrow_u64(cc, nu_13,
				// 	(sp_13 << 1) | (sp_12 >> 63), &nu_13);

				cc = _subborrow_u64(0, sp_0, nv_0, &sp_0);
				cc = _subborrow_u64(cc, sp_1, nv_1, &sp_1);
				cc = _subborrow_u64(cc, sp_2, nv_2, &sp_2);
				cc = _subborrow_u64(cc, sp_3, nv_3, &sp_3);
				cc = _subborrow_u64(cc, sp_4, nv_4, &sp_4);
				cc = _subborrow_u64(cc, sp_5, nv_5, &sp_5);
				cc = _subborrow_u64(cc, sp_6, nv_6, &sp_6);
				cc = _subborrow_u64(cc, sp_7, nv_7, &sp_7);
				cc = _subborrow_u64(cc, sp_8, nv_8, &sp_8);
				(void)_subborrow_u64(cc, sp_9, nv_9, &sp_9);
				// cc = _subborrow_u64(cc, sp_10, nv_10, &sp_10);
				// (void)_subborrow_u64(cc, sp_11, nv_11, &sp_11);
				// cc = _subborrow_u64(cc, sp_12, nv_12, &sp_12);
				// (void)_subborrow_u64(cc, sp_13, nv_13, &sp_13);

			} else {
				unsigned rs, s1, rs1, s2, rs2;

				rs = 64 - s;
				s1 = s + 1;
				rs1 = 64 - s1;
				s2 = s << 1;
				rs2 = 64 - s2;

				cc = _subborrow_u64(0, u0_0,
					(v0_0 << s), &u0_0);
				cc = _subborrow_u64(cc, u0_1,
					(v0_1 << s) | (v0_0 >> rs), &u0_1);
				cc = _subborrow_u64(cc, u0_2,
					(v0_2 << s) | (v0_1 >> rs), &u0_2);
				(void)_subborrow_u64(cc, u0_3,
					(v0_3 << s) | (v0_2 >> rs), &u0_3);


				cc = _subborrow_u64(0, u1_0,
					(v1_0 << s), &u1_0);
				cc = _subborrow_u64(cc, u1_1,
					(v1_1 << s) | (v1_0 >> rs), &u1_1);
				cc = _subborrow_u64(cc, u1_2,
					(v1_2 << s) | (v1_1 >> rs), &u1_2);
				(void)_subborrow_u64(cc, u1_3,
					(v1_3 << s) | (v1_2 >> rs), &u1_3);

				cc = _addcarry_u64(0, nu_0,
					(nv_0 << s2), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(nv_1 << s2) | (nv_0 >> rs2), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(nv_2 << s2) | (nv_1 >> rs2), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(nv_3 << s2) | (nv_2 >> rs2), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(nv_4 << s2) | (nv_3 >> rs2), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(nv_5 << s2) | (nv_4 >> rs2), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(nv_6 << s2) | (nv_5 >> rs2), &nu_6);
				cc = _addcarry_u64(cc, nu_7,
					(nv_7 << s2) | (nv_6 >> rs2), &nu_7);
				cc = _addcarry_u64(cc, nu_8,
					(nv_8 << s2) | (nv_7 >> rs2), &nu_8);
				(void)_addcarry_u64(cc, nu_9,
					(nv_9 << s2) | (nv_8 >> rs2), &nu_9);
				// cc = _addcarry_u64(cc, nu_10,
				// 	(nv_10 << s2) | (nv_9 >> rs2), &nu_10);
				// (void)_addcarry_u64(cc, nu_11,
				// 	(nv_11 << s2) | (nv_10 >> rs2), &nu_11);
				// cc = _addcarry_u64(cc, nu_12,
				// 	(nv_12 << s2) | (nv_11 >> rs2), &nu_12);
				// (void)_addcarry_u64(cc, nu_13,
				// 	(nv_13 << s2) | (nv_12 >> rs2), &nu_13);

				cc = _subborrow_u64(0, nu_0,
					(sp_0 << s1), &nu_0);
				cc = _subborrow_u64(cc, nu_1,
					(sp_1 << s1) | (sp_0 >> rs1), &nu_1);
				cc = _subborrow_u64(cc, nu_2,
					(sp_2 << s1) | (sp_1 >> rs1), &nu_2);
				cc = _subborrow_u64(cc, nu_3,
					(sp_3 << s1) | (sp_2 >> rs1), &nu_3);
				cc = _subborrow_u64(cc, nu_4,
					(sp_4 << s1) | (sp_3 >> rs1), &nu_4);
				cc = _subborrow_u64(cc, nu_5,
					(sp_5 << s1) | (sp_4 >> rs1), &nu_5);
				cc = _subborrow_u64(cc, nu_6,
					(sp_6 << s1) | (sp_5 >> rs1), &nu_6);
				cc = _subborrow_u64(cc, nu_7,
					(sp_7 << s1) | (sp_6 >> rs1), &nu_7);
				cc = _subborrow_u64(cc, nu_8,
					(sp_8 << s1) | (sp_7 >> rs1), &nu_8);
				(void)_subborrow_u64(cc, nu_9,
					(sp_9 << s1) | (sp_8 >> rs1), &nu_9);
				// cc = _subborrow_u64(cc, nu_10,
				// 	(sp_10 << s1) | (sp_9 >> rs1), &nu_10);
				// (void)_subborrow_u64(cc, nu_11,
				// 	(sp_11 << s1) | (sp_10 >> rs1), &nu_11);
				// cc = _subborrow_u64(cc, nu_12,
				// 	(sp_12 << s1) | (sp_11 >> rs1), &nu_12);
				// (void)_subborrow_u64(cc, nu_13,
				// 	(sp_13 << s1) | (sp_12 >> rs1), &nu_13);

				cc = _subborrow_u64(0, sp_0,
					(nv_0 << s), &sp_0);
				cc = _subborrow_u64(cc, sp_1,
					(nv_1 << s) | (nv_0 >> rs), &sp_1);
				cc = _subborrow_u64(cc, sp_2,
					(nv_2 << s) | (nv_1 >> rs), &sp_2);
				cc = _subborrow_u64(cc, sp_3,
					(nv_3 << s) | (nv_2 >> rs), &sp_3);
				cc = _subborrow_u64(cc, sp_4,
					(nv_4 << s) | (nv_3 >> rs), &sp_4);
				cc = _subborrow_u64(cc, sp_5,
					(nv_5 << s) | (nv_4 >> rs), &sp_5);
				cc = _subborrow_u64(cc, sp_6,
					(nv_6 << s) | (nv_5 >> rs), &sp_6);
				cc = _subborrow_u64(cc, sp_7,
					(nv_7 << s) | (nv_6 >> rs), &sp_7);
				cc = _subborrow_u64(cc, sp_8,
					(nv_8 << s) | (nv_7 >> rs), &sp_8);
				(void)_subborrow_u64(cc, sp_9,
					(nv_9 << s) | (nv_8 >> rs), &sp_9);
				// cc = _subborrow_u64(cc, sp_10,
				// 	(nv_10 << s) | (nv_9 >> rs), &sp_10);
				// (void)_subborrow_u64(cc, sp_11,
				// 	(nv_11 << s) | (nv_10 >> rs), &sp_11);
				// cc = _subborrow_u64(cc, sp_12,
				// 	(nv_12 << s) | (nv_11 >> rs), &sp_12);
				// (void)_subborrow_u64(cc, sp_13,
				// 	(nv_13 << s) | (nv_12 >> rs), &sp_13);
			}
		} else {
			if (s == 0) {
				cc = _addcarry_u64(0, u0_0, v0_0, &u0_0);
				cc = _addcarry_u64(cc, u0_1, v0_1, &u0_1);
				cc = _addcarry_u64(cc, u0_2, v0_2, &u0_2);
				(void)_addcarry_u64(cc, u0_3, v0_3, &u0_3);

				cc = _addcarry_u64(0, u1_0, v1_0, &u1_0);
				cc = _addcarry_u64(cc, u1_1, v1_1, &u1_1);
				cc = _addcarry_u64(cc, u1_2, v1_2, &u1_2);
				(void)_addcarry_u64(cc, u1_3, v1_3, &u1_3);

				cc = _addcarry_u64(0, nu_0, nv_0, &nu_0);
				cc = _addcarry_u64(cc, nu_1, nv_1, &nu_1);
				cc = _addcarry_u64(cc, nu_2, nv_2, &nu_2);
				cc = _addcarry_u64(cc, nu_3, nv_3, &nu_3);
				cc = _addcarry_u64(cc, nu_4, nv_4, &nu_4);
				cc = _addcarry_u64(cc, nu_5, nv_5, &nu_5);
				cc = _addcarry_u64(cc, nu_6, nv_6, &nu_6);
				cc = _addcarry_u64(cc, nu_7, nv_7, &nu_7);
				cc = _addcarry_u64(cc, nu_8, nv_8, &nu_8);
				(void)_addcarry_u64(cc, nu_9, nv_9, &nu_9);
				// cc = _addcarry_u64(cc, nu_10, nv_10, &nu_10);
				// (void)_addcarry_u64(cc, nu_11, nv_11, &nu_11);
				// cc = _addcarry_u64(cc, nu_12, nv_12, &nu_12);
				// (void)_addcarry_u64(cc, nu_13, nv_13, &nu_13);

				cc = _addcarry_u64(0, nu_0,
					(sp_0 << 1), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(sp_1 << 1) | (sp_0 >> 63), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(sp_2 << 1) | (sp_1 >> 63), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(sp_3 << 1) | (sp_2 >> 63), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(sp_4 << 1) | (sp_3 >> 63), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(sp_5 << 1) | (sp_4 >> 63), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(sp_6 << 1) | (sp_5 >> 63), &nu_6);
				cc = _addcarry_u64(cc, nu_7,
					(sp_7 << 1) | (sp_6 >> 63), &nu_7);
				cc = _addcarry_u64(cc, nu_8,
					(sp_8 << 1) | (sp_7 >> 63), &nu_8);
				(void)_addcarry_u64(cc, nu_9,
					(sp_9 << 1) | (sp_8 >> 63), &nu_9);
				// cc = _addcarry_u64(cc, nu_10,
				// 	(sp_10 << 1) | (sp_9 >> 63), &nu_10);
				// (void)_addcarry_u64(cc, nu_11,
				// 	(sp_11 << 1) | (sp_10 >> 63), &nu_11);
				// cc = _addcarry_u64(cc, nu_12,
				// 	(sp_12 << 1) | (sp_11 >> 63), &nu_12);
				// (void)_addcarry_u64(cc, nu_13,
				// 	(sp_13 << 1) | (sp_12 >> 63), &nu_13);


				cc = _addcarry_u64(0, sp_0, nv_0, &sp_0);
				cc = _addcarry_u64(cc, sp_1, nv_1, &sp_1);
				cc = _addcarry_u64(cc, sp_2, nv_2, &sp_2);
				cc = _addcarry_u64(cc, sp_3, nv_3, &sp_3);
				cc = _addcarry_u64(cc, sp_4, nv_4, &sp_4);
				cc = _addcarry_u64(cc, sp_5, nv_5, &sp_5);
				cc = _addcarry_u64(cc, sp_6, nv_6, &sp_6);
				cc = _addcarry_u64(cc, sp_7, nv_7, &sp_7);
				cc = _addcarry_u64(cc, sp_8, nv_8, &sp_8);
				(void)_addcarry_u64(cc, sp_9, nv_9, &sp_9);
				// cc = _addcarry_u64(cc, sp_10, nv_10, &sp_10);
				// (void)_addcarry_u64(cc, sp_11, nv_11, &sp_11);
				// cc = _addcarry_u64(cc, sp_12, nv_12, &sp_12);
				// (void)_addcarry_u64(cc, sp_13, nv_13, &sp_13);

			} else {
				unsigned rs, s1, rs1, s2, rs2;

				rs = 64 - s;
				s1 = s + 1;
				rs1 = 64 - s1;
				s2 = s << 1;
				rs2 = 64 - s2;

				cc = _addcarry_u64(0, u0_0,
					(v0_0 << s), &u0_0);
				cc = _addcarry_u64(cc, u0_1,
					(v0_1 << s) | (v0_0 >> rs), &u0_1);
				cc = _addcarry_u64(cc, u0_2,
					(v0_2 << s) | (v0_1 >> rs), &u0_2);
				(void)_addcarry_u64(cc, u0_3,
					(v0_3 << s) | (v0_2 >> rs), &u0_3);


				cc = _addcarry_u64(0, u1_0,
					(v1_0 << s), &u1_0);
				cc = _addcarry_u64(cc, u1_1,
					(v1_1 << s) | (v1_0 >> rs), &u1_1);
				cc = _addcarry_u64(cc, u1_2,
					(v1_2 << s) | (v1_1 >> rs), &u1_2);
				(void)_addcarry_u64(cc, u1_3,
					(v1_3 << s) | (v1_2 >> rs), &u1_3);

				cc = _addcarry_u64(0, nu_0,
					(nv_0 << s2), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(nv_1 << s2) | (nv_0 >> rs2), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(nv_2 << s2) | (nv_1 >> rs2), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(nv_3 << s2) | (nv_2 >> rs2), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(nv_4 << s2) | (nv_3 >> rs2), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(nv_5 << s2) | (nv_4 >> rs2), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(nv_6 << s2) | (nv_5 >> rs2), &nu_6);
				cc = _addcarry_u64(cc, nu_7,
					(nv_7 << s2) | (nv_6 >> rs2), &nu_7);
				cc = _addcarry_u64(cc, nu_8,
					(nv_8 << s2) | (nv_7 >> rs2), &nu_8);
				(void)_addcarry_u64(cc, nu_9,
					(nv_9 << s2) | (nv_8 >> rs2), &nu_9);
				// cc = _addcarry_u64(cc, nu_10,
				// 	(nv_10 << s2) | (nv_9 >> rs2), &nu_10);
				// (void)_addcarry_u64(cc, nu_11,
				// 	(nv_11 << s2) | (nv_10 >> rs2), &nu_11);
				// cc = _addcarry_u64(cc, nu_12,
				// 	(nv_12 << s2) | (nv_11 >> rs2), &nu_12);
				// (void)_addcarry_u64(cc, nu_13,
				// 	(nv_13 << s2) | (nv_12 >> rs2), &nu_13);

				cc = _addcarry_u64(0, nu_0,
					(sp_0 << s1), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(sp_1 << s1) | (sp_0 >> rs1), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(sp_2 << s1) | (sp_1 >> rs1), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(sp_3 << s1) | (sp_2 >> rs1), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(sp_4 << s1) | (sp_3 >> rs1), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(sp_5 << s1) | (sp_4 >> rs1), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(sp_6 << s1) | (sp_5 >> rs1), &nu_6);
				cc = _addcarry_u64(cc, nu_7,
					(sp_7 << s1) | (sp_6 >> rs1), &nu_7);
				cc = _addcarry_u64(cc, nu_8,
					(sp_8 << s1) | (sp_7 >> rs1), &nu_8);
				(void)_addcarry_u64(cc, nu_9,
					(sp_9 << s1) | (sp_8 >> rs1), &nu_9);
				// cc = _addcarry_u64(cc, nu_10,
				// 	(sp_10 << s1) | (sp_9 >> rs1), &nu_10);
				// (void)_addcarry_u64(cc, nu_11,
				// 	(sp_11 << s1) | (sp_10 >> rs1), &nu_11);
				// cc = _addcarry_u64(cc, nu_12,
				// 	(sp_12 << s1) | (sp_11 >> rs1), &nu_12);
				// (void)_addcarry_u64(cc, nu_13,
				// 	(sp_13 << s1) | (sp_12 >> rs1), &nu_13);

				cc = _addcarry_u64(0, sp_0,
					(nv_0 << s), &sp_0);
				cc = _addcarry_u64(cc, sp_1,
					(nv_1 << s) | (nv_0 >> rs), &sp_1);
				cc = _addcarry_u64(cc, sp_2,
					(nv_2 << s) | (nv_1 >> rs), &sp_2);
				cc = _addcarry_u64(cc, sp_3,
					(nv_3 << s) | (nv_2 >> rs), &sp_3);
				cc = _addcarry_u64(cc, sp_4,
					(nv_4 << s) | (nv_3 >> rs), &sp_4);
				cc = _addcarry_u64(cc, sp_5,
					(nv_5 << s) | (nv_4 >> rs), &sp_5);
				cc = _addcarry_u64(cc, sp_6,
					(nv_6 << s) | (nv_5 >> rs), &sp_6);
				cc = _addcarry_u64(cc, sp_7,
					(nv_7 << s) | (nv_6 >> rs), &sp_7);
				cc = _addcarry_u64(cc, sp_8,
					(nv_8 << s) | (nv_7 >> rs), &sp_8);
				(void)_addcarry_u64(cc, sp_9,
					(nv_9 << s) | (nv_8 >> rs), &sp_9);
				// cc = _addcarry_u64(cc, sp_10,
				// 	(nv_10 << s) | (nv_9 >> rs), &sp_10);
				// (void)_addcarry_u64(cc, sp_11,
				// 	(nv_11 << s) | (nv_10 >> rs), &sp_11);
				// cc = _addcarry_u64(cc, sp_12,
				// 	(nv_12 << s) | (nv_11 >> rs), &sp_12);
				// (void)_addcarry_u64(cc, sp_13,
				// 	(nv_13 << s) | (nv_12 >> rs), &sp_13);
			}
		}
	}


	/*
	 * The part below is reached when large values have shrunk to 8
	 * words.
	 */

	for (;;) {
		unsigned s;
		unsigned long long m;
		unsigned bl_nv, bl_sp;

		/*
		 * If nu < nv, then swap(u,v) and swap(nu,nv).
		 */
		cc = _subborrow_u64(0, nu_0, nv_0, &w);
		cc = _subborrow_u64(cc, nu_1, nv_1, &w);
		cc = _subborrow_u64(cc, nu_2, nv_2, &w);
		cc = _subborrow_u64(cc, nu_3, nv_3, &w);
		cc = _subborrow_u64(cc, nu_4, nv_4, &w);
		cc = _subborrow_u64(cc, nu_5, nv_5, &w);
		cc = _subborrow_u64(cc, nu_6, nv_6, &w);
		cc = _subborrow_u64(cc, nu_7, nv_7, &w);
		// cc = _subborrow_u64(cc, nu_8, nv_8, &w);
		// cc = _subborrow_u64(cc, nu_9, nv_9, &w);
		// cc = _subborrow_u64(cc, nu_10, nv_10, &w);
		// cc = _subborrow_u64(cc, nu_11, nv_11, &w);
		// cc = _subborrow_u64(cc, nu_12, nv_12, &w);
		// cc = _subborrow_u64(cc, nu_13, nv_13, &w);
		m = -(unsigned long long)cc;


		SWAPCOND(u0_0, v0_0, m);
		SWAPCOND(u0_1, v0_1, m);
		SWAPCOND(u0_2, v0_2, m);
		SWAPCOND(u0_3, v0_3, m);
		
		SWAPCOND(u1_0, v1_0, m);
		SWAPCOND(u1_1, v1_1, m);
		SWAPCOND(u1_2, v1_2, m);
		SWAPCOND(u1_3, v1_3, m);

		SWAPCOND(nu_0, nv_0, m);
		SWAPCOND(nu_1, nv_1, m);
		SWAPCOND(nu_2, nv_2, m);
		SWAPCOND(nu_3, nv_3, m);
		SWAPCOND(nu_4, nv_4, m);
		SWAPCOND(nu_5, nv_5, m);
		SWAPCOND(nu_6, nv_6, m);
		SWAPCOND(nu_7, nv_7, m);
		// SWAPCOND(nu_8, nv_8, m);
		// SWAPCOND(nu_9, nv_9, m);
		// SWAPCOND(nu_10, nv_10, m);
		// SWAPCOND(nu_11, nv_11, m);
		// SWAPCOND(nu_12, nv_12, m);
		// SWAPCOND(nu_13, nv_13, m);


#define BITLENGTH8(size, bb)   do { \
		unsigned bitlength_acc = 8; \
		int bitlength_flag = 1; \
		unsigned long long bitlength_mask = -(bb ## 7 >> 63); \
		unsigned long long bitlength_top = 0; \
		unsigned long long bitlength_word; \
		bitlength_word = bb ## 7 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 6 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 5 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 4 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 3 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 2 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 1 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		bitlength_flag &= (bitlength_word == 0); \
		bitlength_word = bb ## 0 ^ bitlength_mask; \
		bitlength_top |= \
			bitlength_word & -(unsigned long long)bitlength_flag; \
		bitlength_acc -= bitlength_flag; \
		(size) = (bitlength_acc << 6) - _lzcnt_u64(bitlength_top); \
	} while (0)


		BITLENGTH8(bl_nv, nv_);
		BITLENGTH8(bl_sp, sp_);

		/*
		 * If v is small enough, return.
		 */
		if (bl_nv <= 447) {
			c0[0] = v0_0;
			c0[1] = v0_1;
			c0[2] = v0_2;
			c0[3] = v0_3;
			
			c1[0] = v1_0;
			c1[1] = v1_1;
			c1[2] = v1_2;
			c1[3] = v1_3;

			return;
		}

		/*
		 * Compute shift amount.
		 */
		s = bl_sp - bl_nv;
		s &= ~-(s >> 31);

		/*
		 * It is very rare that s > 31. We handle it with some
		 * generic code; branch prediction will soon learn that
		 * this path is normally not taken.
		 */
		if (s > 31) {
			uint64_t tu0[4], tu1[4], tv0[4], tv1[4];
			uint64_t tnu[14], tnv[14], tsp[14];

			tu0[0] = u0_0;
			tu0[1] = u0_1;
			tu0[2] = u0_2;
			tu0[3] = u0_3;

			tu1[0] = u1_0;
			tu1[1] = u1_1;
			tu1[2] = u1_2;
			tu1[3] = u1_3;

			tv0[0] = v0_0;
			tv0[1] = v0_1;
			tv0[2] = v0_2;
			tv0[3] = v0_3;

			tv1[0] = v1_0;
			tv1[1] = v1_1;
			tv1[2] = v1_2;
			tv1[3] = v1_3;

			tnu[0] = nu_0;
			tnu[1] = nu_1;
			tnu[2] = nu_2;
			tnu[3] = nu_3;
			tnu[4] = nu_4;
			tnu[5] = nu_5;
			tnu[6] = nu_6;
			tnu[7] = nu_7;
			tnu[8] = 0;
			tnu[9] = 0;
			tnu[10] = 0;
			tnu[11] = 0;
			tnu[12] = 0;
			tnu[13] = 0;

			tnv[0] = nv_0;
			tnv[1] = nv_1;
			tnv[2] = nv_2;
			tnv[3] = nv_3;
			tnv[4] = nv_4;
			tnv[5] = nv_5;
			tnv[6] = nv_6;
			tnv[7] = nv_7;
			tnv[8] = 0;
			tnv[9] = 0;
			tnv[10] = 0;
			tnv[11] = 0;
			tnv[12] = 0;
			tnv[13] = 0;

			tsp[0] = sp_0;
			tsp[1] = sp_1;
			tsp[2] = sp_2;
			tsp[3] = sp_3;
			tsp[4] = sp_4;
			tsp[5] = sp_5;
			tsp[6] = sp_6;
			tsp[7] = sp_7;
			tsp[8] = -(sp_7 >> 63);
			tsp[9] = tsp[8];
			tsp[10] = tsp[8];
			tsp[11] = tsp[8];
			tsp[12] = tsp[8];
			tsp[13] = tsp[8];

			if ((sp_13 >> 63) == 0) {
				sub_lshift_small(tu0, tv0, s);
				sub_lshift_small(tu1, tv1, s);
				add_lshift_large(tnu, tnv, 2 * s);
				sub_lshift_large(tnu, tsp, s + 1);
				sub_lshift_large(tsp, tnv, s);
			} else {
				add_lshift_small(tu0, tv0, s);
				add_lshift_small(tu1, tv1, s);
				add_lshift_large(tnu, tnv, 2 * s);
				add_lshift_large(tnu, tsp, s + 1);
				add_lshift_large(tsp, tnv, s);
			}


			u0_0 = tu0[0];
			u0_1 = tu0[1];
			u0_2 = tu0[2];
			u0_3 = tu0[3];

			u1_0 = tu1[0];
			u1_1 = tu1[1];
			u1_2 = tu1[2];
			u1_3 = tu1[3];

			v0_0 = tv0[0];
			v0_1 = tv0[1];
			v0_2 = tv0[2];
			v0_3 = tv0[3];

			v1_0 = tv1[0];
			v1_1 = tv1[1];
			v1_2 = tv1[2];
			v1_3 = tv1[3];

			nu_0 = tnu[0];
			nu_1 = tnu[1];
			nu_2 = tnu[2];
			nu_3 = tnu[3];
			nu_4 = tnu[4];
			nu_5 = tnu[5];
			nu_6 = tnu[6];
			nu_7 = tnu[7];
			// nu_8 = tnu[8];
			// nu_9 = tnu[9];
			// nu_10 = tnu[10];
			// nu_11 = tnu[11];
			// nu_12 = tnu[12];
			// nu_13 = tnu[13];

			

			nv_0 = tnv[0];
			nv_1 = tnv[1];
			nv_2 = tnv[2];
			nv_3 = tnv[3];
			nv_4 = tnv[4];
			nv_5 = tnv[5];
			nv_6 = tnv[6];
			nv_7 = tnv[7];
			// nv_8 = tnv[8];
			// nv_9 = tnv[9];
			// nv_10 = tnv[10];
			// nv_11 = tnv[11];
			// nv_12 = tnv[12];
			// nv_13 = tnv[13];

			sp_0 = tsp[0];
			sp_1 = tsp[1];
			sp_2 = tsp[2];
			sp_3 = tsp[3];
			sp_4 = tsp[4];
			sp_5 = tsp[5];
			sp_6 = tsp[6];
			sp_7 = tsp[7];
			// sp_8 = tsp[8];
			// sp_9 = tsp[9];
			// sp_10 = tsp[10];
			// sp_11 = tsp[11];
			// sp_12 = tsp[12];
			// sp_13 = tsp[13];
			continue;
		}

		if ((sp_7 >> 63) == 0) {
			if (s == 0) {
				cc = _subborrow_u64(0, u0_0, v0_0, &u0_0);
				cc = _subborrow_u64(cc, u0_1, v0_1, &u0_1);
				cc = _subborrow_u64(cc, u0_2, v0_2, &u0_2);
				(void)_subborrow_u64(cc, u0_3, v0_3, &u0_3);

				cc = _subborrow_u64(0, u1_0, v1_0, &u1_0);
				cc = _subborrow_u64(cc, u1_1, v1_1, &u1_1);
				cc = _subborrow_u64(cc, u1_2, v1_2, &u1_2);
				(void)_subborrow_u64(cc, u1_3, v1_3, &u1_3);

				cc = _addcarry_u64(0, nu_0, nv_0, &nu_0);
				cc = _addcarry_u64(cc, nu_1, nv_1, &nu_1);
				cc = _addcarry_u64(cc, nu_2, nv_2, &nu_2);
				cc = _addcarry_u64(cc, nu_3, nv_3, &nu_3);
				cc = _addcarry_u64(cc, nu_4, nv_4, &nu_4);
				cc = _addcarry_u64(cc, nu_5, nv_5, &nu_5);
				cc = _addcarry_u64(cc, nu_6, nv_6, &nu_6);
				(void)_addcarry_u64(cc, nu_7, nv_7, &nu_7);
				// cc = _addcarry_u64(cc, nu_8, nv_8, &nu_8);
				// (void)_addcarry_u64(cc, nu_9, nv_9, &nu_9);
				// cc = _addcarry_u64(cc, nu_10, nv_10, &nu_10);
				// (void)_addcarry_u64(cc, nu_11, nv_11, &nu_11);
				// cc = _addcarry_u64(cc, nu_12, nv_12, &nu_12);
				// (void)_addcarry_u64(cc, nu_13, nv_13, &nu_13);

				cc = _subborrow_u64(0, nu_0,
					(sp_0 << 1), &nu_0);
				cc = _subborrow_u64(cc, nu_1,
					(sp_1 << 1) | (sp_0 >> 63), &nu_1);
				cc = _subborrow_u64(cc, nu_2,
					(sp_2 << 1) | (sp_1 >> 63), &nu_2);
				cc = _subborrow_u64(cc, nu_3,
					(sp_3 << 1) | (sp_2 >> 63), &nu_3);
				cc = _subborrow_u64(cc, nu_4,
					(sp_4 << 1) | (sp_3 >> 63), &nu_4);
				cc = _subborrow_u64(cc, nu_5,
					(sp_5 << 1) | (sp_4 >> 63), &nu_5);
				cc = _subborrow_u64(cc, nu_6,
					(sp_6 << 1) | (sp_5 >> 63), &nu_6);
				(void)_subborrow_u64(cc, nu_7,
					(sp_7 << 1) | (sp_6 >> 63), &nu_7);
				// cc = _subborrow_u64(cc, nu_8,
				// 	(sp_8 << 1) | (sp_7 >> 63), &nu_8);
				// (void)_subborrow_u64(cc, nu_9,
				// 	(sp_9 << 1) | (sp_8 >> 63), &nu_9);
				// cc = _subborrow_u64(cc, nu_10,
				// 	(sp_10 << 1) | (sp_9 >> 63), &nu_10);
				// (void)_subborrow_u64(cc, nu_11,
				// 	(sp_11 << 1) | (sp_10 >> 63), &nu_11);
				// cc = _subborrow_u64(cc, nu_12,
				// 	(sp_12 << 1) | (sp_11 >> 63), &nu_12);
				// (void)_subborrow_u64(cc, nu_13,
				// 	(sp_13 << 1) | (sp_12 >> 63), &nu_13);

				cc = _subborrow_u64(0, sp_0, nv_0, &sp_0);
				cc = _subborrow_u64(cc, sp_1, nv_1, &sp_1);
				cc = _subborrow_u64(cc, sp_2, nv_2, &sp_2);
				cc = _subborrow_u64(cc, sp_3, nv_3, &sp_3);
				cc = _subborrow_u64(cc, sp_4, nv_4, &sp_4);
				cc = _subborrow_u64(cc, sp_5, nv_5, &sp_5);
				cc = _subborrow_u64(cc, sp_6, nv_6, &sp_6);
				(void)_subborrow_u64(cc, sp_7, nv_7, &sp_7);
				// cc = _subborrow_u64(cc, sp_8, nv_8, &sp_8);
				// (void)_subborrow_u64(cc, sp_9, nv_9, &sp_9);
				// cc = _subborrow_u64(cc, sp_10, nv_10, &sp_10);
				// (void)_subborrow_u64(cc, sp_11, nv_11, &sp_11);
				// cc = _subborrow_u64(cc, sp_12, nv_12, &sp_12);
				// (void)_subborrow_u64(cc, sp_13, nv_13, &sp_13);

			} else {
				unsigned rs, s1, rs1, s2, rs2;

				rs = 64 - s;
				s1 = s + 1;
				rs1 = 64 - s1;
				s2 = s << 1;
				rs2 = 64 - s2;

				cc = _subborrow_u64(0, u0_0,
					(v0_0 << s), &u0_0);
				cc = _subborrow_u64(cc, u0_1,
					(v0_1 << s) | (v0_0 >> rs), &u0_1);
				cc = _subborrow_u64(cc, u0_2,
					(v0_2 << s) | (v0_1 >> rs), &u0_2);
				(void)_subborrow_u64(cc, u0_3,
					(v0_3 << s) | (v0_2 >> rs), &u0_3);


				cc = _subborrow_u64(0, u1_0,
					(v1_0 << s), &u1_0);
				cc = _subborrow_u64(cc, u1_1,
					(v1_1 << s) | (v1_0 >> rs), &u1_1);
				cc = _subborrow_u64(cc, u1_2,
					(v1_2 << s) | (v1_1 >> rs), &u1_2);
				(void)_subborrow_u64(cc, u1_3,
					(v1_3 << s) | (v1_2 >> rs), &u1_3);

				cc = _addcarry_u64(0, nu_0,
					(nv_0 << s2), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(nv_1 << s2) | (nv_0 >> rs2), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(nv_2 << s2) | (nv_1 >> rs2), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(nv_3 << s2) | (nv_2 >> rs2), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(nv_4 << s2) | (nv_3 >> rs2), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(nv_5 << s2) | (nv_4 >> rs2), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(nv_6 << s2) | (nv_5 >> rs2), &nu_6);
				(void)_addcarry_u64(cc, nu_7,
					(nv_7 << s2) | (nv_6 >> rs2), &nu_7);
				// cc = _addcarry_u64(cc, nu_8,
				// 	(nv_8 << s2) | (nv_7 >> rs2), &nu_8);
				// (void)_addcarry_u64(cc, nu_9,
				// 	(nv_9 << s2) | (nv_8 >> rs2), &nu_9);
				// cc = _addcarry_u64(cc, nu_10,
				// 	(nv_10 << s2) | (nv_9 >> rs2), &nu_10);
				// (void)_addcarry_u64(cc, nu_11,
				// 	(nv_11 << s2) | (nv_10 >> rs2), &nu_11);
				// cc = _addcarry_u64(cc, nu_12,
				// 	(nv_12 << s2) | (nv_11 >> rs2), &nu_12);
				// (void)_addcarry_u64(cc, nu_13,
				// 	(nv_13 << s2) | (nv_12 >> rs2), &nu_13);

				cc = _subborrow_u64(0, nu_0,
					(sp_0 << s1), &nu_0);
				cc = _subborrow_u64(cc, nu_1,
					(sp_1 << s1) | (sp_0 >> rs1), &nu_1);
				cc = _subborrow_u64(cc, nu_2,
					(sp_2 << s1) | (sp_1 >> rs1), &nu_2);
				cc = _subborrow_u64(cc, nu_3,
					(sp_3 << s1) | (sp_2 >> rs1), &nu_3);
				cc = _subborrow_u64(cc, nu_4,
					(sp_4 << s1) | (sp_3 >> rs1), &nu_4);
				cc = _subborrow_u64(cc, nu_5,
					(sp_5 << s1) | (sp_4 >> rs1), &nu_5);
				cc = _subborrow_u64(cc, nu_6,
					(sp_6 << s1) | (sp_5 >> rs1), &nu_6);
				(void)_subborrow_u64(cc, nu_7,
					(sp_7 << s1) | (sp_6 >> rs1), &nu_7);
				// cc = _subborrow_u64(cc, nu_8,
				// 	(sp_8 << s1) | (sp_7 >> rs1), &nu_8);
				// (void)_subborrow_u64(cc, nu_9,
				// 	(sp_9 << s1) | (sp_8 >> rs1), &nu_9);
				// cc = _subborrow_u64(cc, nu_10,
				// 	(sp_10 << s1) | (sp_9 >> rs1), &nu_10);
				// (void)_subborrow_u64(cc, nu_11,
				// 	(sp_11 << s1) | (sp_10 >> rs1), &nu_11);
				// cc = _subborrow_u64(cc, nu_12,
				// 	(sp_12 << s1) | (sp_11 >> rs1), &nu_12);
				// (void)_subborrow_u64(cc, nu_13,
				// 	(sp_13 << s1) | (sp_12 >> rs1), &nu_13);

				cc = _subborrow_u64(0, sp_0,
					(nv_0 << s), &sp_0);
				cc = _subborrow_u64(cc, sp_1,
					(nv_1 << s) | (nv_0 >> rs), &sp_1);
				cc = _subborrow_u64(cc, sp_2,
					(nv_2 << s) | (nv_1 >> rs), &sp_2);
				cc = _subborrow_u64(cc, sp_3,
					(nv_3 << s) | (nv_2 >> rs), &sp_3);
				cc = _subborrow_u64(cc, sp_4,
					(nv_4 << s) | (nv_3 >> rs), &sp_4);
				cc = _subborrow_u64(cc, sp_5,
					(nv_5 << s) | (nv_4 >> rs), &sp_5);
				cc = _subborrow_u64(cc, sp_6,
					(nv_6 << s) | (nv_5 >> rs), &sp_6);
				(void)_subborrow_u64(cc, sp_7,
					(nv_7 << s) | (nv_6 >> rs), &sp_7);
				// cc = _subborrow_u64(cc, sp_8,
				// 	(nv_8 << s) | (nv_7 >> rs), &sp_8);
				// (void)_subborrow_u64(cc, sp_9,
				// 	(nv_9 << s) | (nv_8 >> rs), &sp_9);
				// cc = _subborrow_u64(cc, sp_10,
				// 	(nv_10 << s) | (nv_9 >> rs), &sp_10);
				// (void)_subborrow_u64(cc, sp_11,
				// 	(nv_11 << s) | (nv_10 >> rs), &sp_11);
				// cc = _subborrow_u64(cc, sp_12,
				// 	(nv_12 << s) | (nv_11 >> rs), &sp_12);
				// (void)_subborrow_u64(cc, sp_13,
				// 	(nv_13 << s) | (nv_12 >> rs), &sp_13);
			}
		} else {
			if (s == 0) {
				cc = _addcarry_u64(0, u0_0, v0_0, &u0_0);
				cc = _addcarry_u64(cc, u0_1, v0_1, &u0_1);
				cc = _addcarry_u64(cc, u0_2, v0_2, &u0_2);
				(void)_addcarry_u64(cc, u0_3, v0_3, &u0_3);

				cc = _addcarry_u64(0, u1_0, v1_0, &u1_0);
				cc = _addcarry_u64(cc, u1_1, v1_1, &u1_1);
				cc = _addcarry_u64(cc, u1_2, v1_2, &u1_2);
				(void)_addcarry_u64(cc, u1_3, v1_3, &u1_3);

				cc = _addcarry_u64(0, nu_0, nv_0, &nu_0);
				cc = _addcarry_u64(cc, nu_1, nv_1, &nu_1);
				cc = _addcarry_u64(cc, nu_2, nv_2, &nu_2);
				cc = _addcarry_u64(cc, nu_3, nv_3, &nu_3);
				cc = _addcarry_u64(cc, nu_4, nv_4, &nu_4);
				cc = _addcarry_u64(cc, nu_5, nv_5, &nu_5);
				cc = _addcarry_u64(cc, nu_6, nv_6, &nu_6);
				(void)_addcarry_u64(cc, nu_7, nv_7, &nu_7);
				// cc = _addcarry_u64(cc, nu_8, nv_8, &nu_8);
				// (void)_addcarry_u64(cc, nu_9, nv_9, &nu_9);
				// cc = _addcarry_u64(cc, nu_10, nv_10, &nu_10);
				// (void)_addcarry_u64(cc, nu_11, nv_11, &nu_11);
				// cc = _addcarry_u64(cc, nu_12, nv_12, &nu_12);
				// (void)_addcarry_u64(cc, nu_13, nv_13, &nu_13);

				cc = _addcarry_u64(0, nu_0,
					(sp_0 << 1), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(sp_1 << 1) | (sp_0 >> 63), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(sp_2 << 1) | (sp_1 >> 63), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(sp_3 << 1) | (sp_2 >> 63), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(sp_4 << 1) | (sp_3 >> 63), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(sp_5 << 1) | (sp_4 >> 63), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(sp_6 << 1) | (sp_5 >> 63), &nu_6);
				(void)_addcarry_u64(cc, nu_7,
					(sp_7 << 1) | (sp_6 >> 63), &nu_7);
				// cc = _addcarry_u64(cc, nu_8,
				// 	(sp_8 << 1) | (sp_7 >> 63), &nu_8);
				// (void)_addcarry_u64(cc, nu_9,
				// 	(sp_9 << 1) | (sp_8 >> 63), &nu_9);
				// cc = _addcarry_u64(cc, nu_10,
				// 	(sp_10 << 1) | (sp_9 >> 63), &nu_10);
				// (void)_addcarry_u64(cc, nu_11,
				// 	(sp_11 << 1) | (sp_10 >> 63), &nu_11);
				// cc = _addcarry_u64(cc, nu_12,
				// 	(sp_12 << 1) | (sp_11 >> 63), &nu_12);
				// (void)_addcarry_u64(cc, nu_13,
				// 	(sp_13 << 1) | (sp_12 >> 63), &nu_13);


				cc = _addcarry_u64(0, sp_0, nv_0, &sp_0);
				cc = _addcarry_u64(cc, sp_1, nv_1, &sp_1);
				cc = _addcarry_u64(cc, sp_2, nv_2, &sp_2);
				cc = _addcarry_u64(cc, sp_3, nv_3, &sp_3);
				cc = _addcarry_u64(cc, sp_4, nv_4, &sp_4);
				cc = _addcarry_u64(cc, sp_5, nv_5, &sp_5);
				cc = _addcarry_u64(cc, sp_6, nv_6, &sp_6);
				(void)_addcarry_u64(cc, sp_7, nv_7, &sp_7);
				// cc = _addcarry_u64(cc, sp_8, nv_8, &sp_8);
				// (void)_addcarry_u64(cc, sp_9, nv_9, &sp_9);
				// cc = _addcarry_u64(cc, sp_10, nv_10, &sp_10);
				// (void)_addcarry_u64(cc, sp_11, nv_11, &sp_11);
				// cc = _addcarry_u64(cc, sp_12, nv_12, &sp_12);
				// (void)_addcarry_u64(cc, sp_13, nv_13, &sp_13);

			} else {
				unsigned rs, s1, rs1, s2, rs2;

				rs = 64 - s;
				s1 = s + 1;
				rs1 = 64 - s1;
				s2 = s << 1;
				rs2 = 64 - s2;

				cc = _addcarry_u64(0, u0_0,
					(v0_0 << s), &u0_0);
				cc = _addcarry_u64(cc, u0_1,
					(v0_1 << s) | (v0_0 >> rs), &u0_1);
				cc = _addcarry_u64(cc, u0_2,
					(v0_2 << s) | (v0_1 >> rs), &u0_2);
				(void)_addcarry_u64(cc, u0_3,
					(v0_3 << s) | (v0_2 >> rs), &u0_3);


				cc = _addcarry_u64(0, u1_0,
					(v1_0 << s), &u1_0);
				cc = _addcarry_u64(cc, u1_1,
					(v1_1 << s) | (v1_0 >> rs), &u1_1);
				cc = _addcarry_u64(cc, u1_2,
					(v1_2 << s) | (v1_1 >> rs), &u1_2);
				(void)_addcarry_u64(cc, u1_3,
					(v1_3 << s) | (v1_2 >> rs), &u1_3);

				cc = _addcarry_u64(0, nu_0,
					(nv_0 << s2), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(nv_1 << s2) | (nv_0 >> rs2), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(nv_2 << s2) | (nv_1 >> rs2), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(nv_3 << s2) | (nv_2 >> rs2), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(nv_4 << s2) | (nv_3 >> rs2), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(nv_5 << s2) | (nv_4 >> rs2), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(nv_6 << s2) | (nv_5 >> rs2), &nu_6);
				(void)_addcarry_u64(cc, nu_7,
					(nv_7 << s2) | (nv_6 >> rs2), &nu_7);
				// cc = _addcarry_u64(cc, nu_8,
				// 	(nv_8 << s2) | (nv_7 >> rs2), &nu_8);
				// (void)_addcarry_u64(cc, nu_9,
				// 	(nv_9 << s2) | (nv_8 >> rs2), &nu_9);
				// cc = _addcarry_u64(cc, nu_10,
				// 	(nv_10 << s2) | (nv_9 >> rs2), &nu_10);
				// (void)_addcarry_u64(cc, nu_11,
				// 	(nv_11 << s2) | (nv_10 >> rs2), &nu_11);
				// cc = _addcarry_u64(cc, nu_12,
				// 	(nv_12 << s2) | (nv_11 >> rs2), &nu_12);
				// (void)_addcarry_u64(cc, nu_13,
				// 	(nv_13 << s2) | (nv_12 >> rs2), &nu_13);

				cc = _addcarry_u64(0, nu_0,
					(sp_0 << s1), &nu_0);
				cc = _addcarry_u64(cc, nu_1,
					(sp_1 << s1) | (sp_0 >> rs1), &nu_1);
				cc = _addcarry_u64(cc, nu_2,
					(sp_2 << s1) | (sp_1 >> rs1), &nu_2);
				cc = _addcarry_u64(cc, nu_3,
					(sp_3 << s1) | (sp_2 >> rs1), &nu_3);
				cc = _addcarry_u64(cc, nu_4,
					(sp_4 << s1) | (sp_3 >> rs1), &nu_4);
				cc = _addcarry_u64(cc, nu_5,
					(sp_5 << s1) | (sp_4 >> rs1), &nu_5);
				cc = _addcarry_u64(cc, nu_6,
					(sp_6 << s1) | (sp_5 >> rs1), &nu_6);
				(void)_addcarry_u64(cc, nu_7,
					(sp_7 << s1) | (sp_6 >> rs1), &nu_7);
				// cc = _addcarry_u64(cc, nu_8,
				// 	(sp_8 << s1) | (sp_7 >> rs1), &nu_8);
				// (void)_addcarry_u64(cc, nu_9,
				// 	(sp_9 << s1) | (sp_8 >> rs1), &nu_9);
				// cc = _addcarry_u64(cc, nu_10,
				// 	(sp_10 << s1) | (sp_9 >> rs1), &nu_10);
				// (void)_addcarry_u64(cc, nu_11,
				// 	(sp_11 << s1) | (sp_10 >> rs1), &nu_11);
				// cc = _addcarry_u64(cc, nu_12,
				// 	(sp_12 << s1) | (sp_11 >> rs1), &nu_12);
				// (void)_addcarry_u64(cc, nu_13,
				// 	(sp_13 << s1) | (sp_12 >> rs1), &nu_13);

				cc = _addcarry_u64(0, sp_0,
					(nv_0 << s), &sp_0);
				cc = _addcarry_u64(cc, sp_1,
					(nv_1 << s) | (nv_0 >> rs), &sp_1);
				cc = _addcarry_u64(cc, sp_2,
					(nv_2 << s) | (nv_1 >> rs), &sp_2);
				cc = _addcarry_u64(cc, sp_3,
					(nv_3 << s) | (nv_2 >> rs), &sp_3);
				cc = _addcarry_u64(cc, sp_4,
					(nv_4 << s) | (nv_3 >> rs), &sp_4);
				cc = _addcarry_u64(cc, sp_5,
					(nv_5 << s) | (nv_4 >> rs), &sp_5);
				cc = _addcarry_u64(cc, sp_6,
					(nv_6 << s) | (nv_5 >> rs), &sp_6);
				(void)_addcarry_u64(cc, sp_7,
					(nv_7 << s) | (nv_6 >> rs), &sp_7);
				// cc = _addcarry_u64(cc, sp_8,
				// 	(nv_8 << s) | (nv_7 >> rs), &sp_8);
				// (void)_addcarry_u64(cc, sp_9,
				// 	(nv_9 << s) | (nv_8 >> rs), &sp_9);
				// cc = _addcarry_u64(cc, sp_10,
				// 	(nv_10 << s) | (nv_9 >> rs), &sp_10);
				// (void)_addcarry_u64(cc, sp_11,
				// 	(nv_11 << s) | (nv_10 >> rs), &sp_11);
				// cc = _addcarry_u64(cc, sp_12,
				// 	(nv_12 << s) | (nv_11 >> rs), &sp_12);
				// (void)_addcarry_u64(cc, sp_13,
				// 	(nv_13 << s) | (nv_12 >> rs), &sp_13);
			}
		}
	}


#undef SWAPCOND
#undef BITLENGTH
#undef BITLENGTH12
#undef BITLENGTH10
#undef BITLENGTH8
}