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

#include "curve25519_hgcd_vartime.h"
#include "gmp-impl_custom.h"


/* custom hgcd */
inline mp_size_t
mpn_hgcd_custom (mp_ptr ap, mp_ptr bp, mp_size_t n, mp_size_t s,
	  struct hgcd_matrix *M, mp_ptr tp)
{
    // mp_size_t s = n/2 + 1;
    mp_size_t nn;
    int success = 0;

    for (;;)
    {
      nn = mpn_hgcd_step (n, ap, bp, s, M, tp);
      if (!nn)
	    return success ? n : 0;

      n = nn;
      success = 1;
    }
}


/* return rp = abs(ap - bp) and t = sign(ap - bp) * -(u00 + u01) */
mp_size_t final_step(mp_ptr rp, mp_ptr ap, mp_ptr bp, mp_size_t n, struct hgcd_matrix *M, mp_ptr tp) 
{
    int c;
    int swapped = 0;
    mp_size_t tn = M->n;
    c = mpn_cmp(ap, bp, n);
    if(c < 0){
        MP_PTR_SWAP(ap, bp);
        swapped = 1;
    }
    mpn_sub(rp, ap, n, bp, n);
    mpn_add(tp, M->p[0][0], tn, M->p[0][1], tn);
    MPN_NORMALIZE(tp, tn);
    /* convert tp from unsigned to signed */
    if (tp[tn-1] >> 63){
      tn += 1;
      tp[tn-1] = 0;
    }
    if(!swapped)
        mpn_neg(tp, tp, tn);
    
    return tn;

}



/* 
 * Return signed integers rho and tau s.t. len(rho) and len(tau) ~ len(el)/2, and
 * rho = tau * v mod el 
*/
void
curve25519_hgcd_vartime(
	uint64_t *rho, uint64_t *tau, const uint64_t *v, int *rho_size, int *tau_size)
{
    mp_size_t n = 4; // Number of limbs
    
    /* ap = el */
    mp_limb_t ap[] = {
        0x5812631a5cf5d3ed,
		0x14def9dea2f79cd6,
		0x0000000000000000,
		0x1000000000000000
    };
    
    /* bp = v */
    mp_limb_t bp[n];
    bp[0] = v[0];
    bp[1] = v[1];
    bp[2] = v[2];
    bp[3] = v[3];

    
    mp_size_t s_matrix = (n+1)/2 + 1;
    struct hgcd_matrix M;        
    mp_limb_t tmp[4 * n];         
    mp_limb_t tmp_matrix[4 * s_matrix];
    mp_size_t n_reduced;
    
    mpn_hgcd_matrix_init(&M, n, tmp_matrix);
   
    /*
     (c;d) = M^{-1}(a;b)
     M = [u00 u01]
         [u10 u11]
     size(u**) ~< N/2
     c =  u11 * a - u01 * b
     d = -u10 * a + u00 * b
     size(c) = size(d) > s
     size(abs(c-d)) = s
     */
    n_reduced = mpn_hgcd(ap, bp, n, &M, tmp);

    mpz_t r2, r1, t2, t1, q;
    mpz_inits(r2, r1, t2, t1, q, NULL);

    int c;
    c = mpn_cmp(ap, bp, n_reduced);
    if(c>0){
        mp_limb_t_2_mpz(r2, ap, n);
        mp_limb_t_2_mpz(r1, bp, n);
        mp_limb_t_2_mpz(t2, M.p[0][1], M.n);
        mpz_neg(t2, t2);
        mp_limb_t_2_mpz(t1, M.p[0][0], M.n);
    }else{
        mp_limb_t_2_mpz(r1, ap, n);
        mp_limb_t_2_mpz(r2, bp, n);
        mp_limb_t_2_mpz(t1, M.p[0][1], M.n);
        mpz_neg(t1, t1);
        mp_limb_t_2_mpz(t2, M.p[0][0], M.n);
    }

    for(;;){
    	if (mpz_sizeinbase(r1, 2) <= 127) {
            for(int i=0; i<r1->_mp_size;i++){
                rho[i] = r1->_mp_d[i];
            }
            *rho_size = r1->_mp_size;

            if (t1->_mp_size < 0 ){
                t1->_mp_size *= -1;
                for(int i=0; i< t1->_mp_size;i++){
                    t1->_mp_d[i] = ~ t1->_mp_d[i];
                }
                mpz_add_ui(t1, t1, 1);
            }
            *tau_size = t1->_mp_size;
            for(int i=0; i< *tau_size;i++){
                tau[i] = t1->_mp_d[i];
            }
        
			mpz_clears(q, r1, r2, t1, t2, NULL);
			return;
		}
        mpz_tdiv_qr(q, r2, r2, r1);
        mpz_swap(r2, r1);
        mpz_submul(t2, q, t1);
        mpz_swap(t2, t1);
    }
}


/* 
 * Return signed integers rho and tau s.t. len(rho) and len(tau) ~ len(el)/2, and
 * rho = tau * v mod el 
*/
void
curve25519_hgcd_vartime_enhance1(
	uint64_t *rho, uint64_t *tau, const uint64_t *v, int *rho_size, int *tau_size)
{
    mp_size_t n = 4; // Number of limbs
    
    /* ap = el */
    mp_limb_t ap[] = {
        0x5812631a5cf5d3ed,
		0x14def9dea2f79cd6,
		0x0000000000000000,
		0x1000000000000000
    };
    
    /* bp = v */
    mp_limb_t bp[n];
    bp[0] = v[0];
    bp[1] = v[1];
    bp[2] = v[2];
    bp[3] = v[3];

    
    mp_size_t s = 2; // s = floor(n/2)
    mp_size_t s_matrix = (n+1)/2 + 1;
    struct hgcd_matrix M;        
    mp_limb_t tmp[4 * n];         
    mp_limb_t tmp_matrix[4 * s_matrix];
    mp_size_t n_reduced;
    
    mpn_hgcd_matrix_init(&M, n, tmp_matrix);
   
    /*
     (c;d) = M^{-1}(a;b)
     M = [u00 u01]
         [u10 u11]
     size(u**) ~< N/2
     c =  u11 * a - u01 * b
     d = -u10 * a + u00 * b
     size(c) = size(d) > s
     size(abs(c-d)) = s
     */
    n_reduced = mpn_hgcd_custom(ap, bp, n, s, &M, tmp);

    mpz_t r2, r1, t2, t1, q;
    mpz_inits(r2, r1, t2, t1, q, NULL);

    int c;
    c = mpn_cmp(ap, bp, n_reduced);
    if(c>0){
        mp_limb_t_2_mpz(r2, ap, n);
        mp_limb_t_2_mpz(r1, bp, n);
        mp_limb_t_2_mpz(t2, M.p[0][1], M.n);
        mpz_neg(t2, t2);
        mp_limb_t_2_mpz(t1, M.p[0][0], M.n);
    }else{
        mp_limb_t_2_mpz(r1, ap, n);
        mp_limb_t_2_mpz(r2, bp, n);
        mp_limb_t_2_mpz(t1, M.p[0][1], M.n);
        mpz_neg(t1, t1);
        mp_limb_t_2_mpz(t2, M.p[0][0], M.n);
    }

    for(;;){
    	if (mpz_sizeinbase(r1, 2) <= 127) {
            for(int i=0; i<r1->_mp_size;i++){
                rho[i] = r1->_mp_d[i];
            }
            *rho_size = r1->_mp_size;

            if (t1->_mp_size < 0 ){
                t1->_mp_size *= -1;
                for(int i=0; i< t1->_mp_size;i++){
                    t1->_mp_d[i] = ~ t1->_mp_d[i];
                }
                mpz_add_ui(t1, t1, 1);
            }
            *tau_size = t1->_mp_size;
            for(int i=0; i< *tau_size;i++){
                tau[i] = t1->_mp_d[i];
            }
			mpz_clears(q, r1, r2, t1, t2, NULL);
			return;
		}
        mpz_tdiv_qr(q, r2, r2, r1);
        mpz_swap(r2, r1);
        mpz_submul(t2, q, t1);
        mpz_swap(t2, t1);
    }
}


/* 
 * Return signed integers rho and tau s.t. len(rho) and len(tau) ~ len(el)/2, and
 * rho = tau * v mod el 
*/
void
curve25519_hgcd_vartime_enhance2(
	uint64_t *rho, uint64_t *tau, const uint64_t *v, int *rho_size, int *tau_size)
{
    mp_size_t n = 4; // Number of limbs
    
    /* ap = el */
    mp_limb_t ap[] = {
        0x5812631a5cf5d3ed,
		0x14def9dea2f79cd6,
		0x0000000000000000,
		0x1000000000000000
    };
    
    /* bp = v */
    mp_limb_t bp[n];
    bp[0] = v[0];
    bp[1] = v[1];
    bp[2] = v[2];
    bp[3] = v[3];

    
    mp_size_t s = 2; // s = floor(n/2)
    mp_size_t s_matrix = (n+1)/2 + 1;
    struct hgcd_matrix M;        
    mp_limb_t tmp[4 * n];         
    mp_limb_t tmp_matrix[4 * s_matrix];
    mp_size_t n_reduced;

    mp_limb_t rp[n];
    mp_limb_t tp[s];
    mp_size_t tn;
    
    mpn_hgcd_matrix_init(&M, n, tmp_matrix);
   
    /*
     (c;d) = M^{-1}(a;b)
     M = [u00 u01]
         [u10 u11]
     size(u**) ~< N/2
     c =  u11 * a - u01 * b
     d = -u10 * a + u00 * b
     size(c) = size(d) > s
     size(abs(c-d)) = s
     */
    n_reduced = mpn_hgcd_custom(ap, bp, n, s, &M, tmp);

    /*
    => c - b = (u11 + u10) * a - (u01 + u00) * b
     r = abs(c - d); size(r) = s
     t = sign(c-b) * -(u01 + u00)
    => r == t * b mod a
    */
    tn = final_step(rp, ap, bp, n_reduced, &M, tp);
    MPN_NORMALIZE(rp, n_reduced);

    /* convert rp from unsigned to signed */
    if (rp[n_reduced-1] >> 63){
      n_reduced += 1;
      rp[n_reduced-1] = 0;
    }
    
    /* set rho = r */
    for (mp_size_t i=0; i<n_reduced; i++){
        rho[i] = rp[i];
    }
    *rho_size = n_reduced;
    
    /* set tau = t */
    for (mp_size_t i=0; i<tn; i++){
        tau[i] = tp[i];
    }
    *tau_size = tn;
    return;
}
