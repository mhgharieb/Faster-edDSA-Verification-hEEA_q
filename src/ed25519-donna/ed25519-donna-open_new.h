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

/* New single verification using QSM_B_B' with hEEA */
int
ED25519_FN(ed25519_sign_open_hEEA) (const unsigned char *m, size_t mlen, const ed25519_public_key pk, const ed25519_signature RS) {
	ge25519 ALIGN(16) R, A, sumBRA;
	hash_512bits hash;
	bignum256modm hram, S1, S2={0}, rho, tau;
    int tau_isNegative;
	int rho_isNegative;

	if ((RS[63] & 224))
		return -1;

	/* hram <-- H(R,A,m) */
	ed25519_hram(hash, RS, pk, m, mlen);
	expand256_modm(hram, hash, 64);

	/* compute tau and rho s.t. tau * h = rho mod el */
    curve25519_half_size_scalar_vartime_hEEA(rho, tau, hram, &rho_isNegative, &tau_isNegative);
    
    /* unpacking (-R) */
    if (!ge25519_unpack_negative_vartime(&R, RS))
		return -1;

	/*
		[S]B + (-R) + h(-A) =? 0 ==> [tau * S]B + [tau](-R) + [rho](-A) =? 0
		----------------------------------------------------------------------------------
		tau_isNegative | rho_isNegative | [tau * s]B + [tau](-R) + [rho](-A) =? 0
		----------------------------------------------------------------------------------
		       0   	   |        0       | [|tau| * s]B + [|tau|](-R) + [|rho|](-A)
		----------------------------------------------------------------------------------
		       0   	   |        1       | [|tau| * s]B + [|tau|](-R) - [|rho|](-A) ==>
		               |                | [|tau| * s]B + [|tau|](-R) + [|rho|](A)
		----------------------------------------------------------------------------------
		       1       |        0       | [-|tau| * s]B + [-|tau|](-R) + [|rho|](-A) ==> 
		               |                | [|tau| * s] B + [|tau|](-R) + [|rho|](A)
		----------------------------------------------------------------------------------
		       1       |        1       | [-|tau| * s]B + [-|tau|](-R) + [-|rho|](-A) ==> 
		               |                | [|tau| * s] B + [|tau|](-R) + [|rho|](-A)
		----------------------------------------------------------------------------------
	*/
	/* unpacking (-A) and adjust the sign based on the previous table */
    if (tau_isNegative == rho_isNegative){
		if (!ge25519_unpack_negative_vartime(&A, pk))
			return -1;
	}else{
		if (!ge25519_unpack_positive_vartime(&A, pk))
			return -1;
	}
    
    /* S */
	expand256_modm(S1, RS + 32, 32);

    /* S <-- tau * S */
    mul256_modm(S1, S1, tau);
    
    /* split S to S1 and S2, s.t. S = (S2<<126 + S1), bl(S1) = 126, bl(S2) = 127 */
    S2[0] = (S1[2] >> 14) | ((S1[3] & 0x3FFF) << 42);
    S2[1] = (S1[3] >> 14) | ((S1[4] & 0x3FFF) << 42);
    S2[2] = S1[4] >> 14;
    S1[2] &= 0x3FFF;
    S1[3] = 0;
    S1[4] = 0;

	/* 
	 * [S]B + (-R) + h(-A) =? 0 ==>
	 * [tau * S]B + [tau](-R) + [rho](-A) =? 0 ==>
	 * [S1]B + [S2]([2^126]B) + [tau](-R) + [rho](-A) =? 0 
	*/
	ge25519_quadruple_scalarmult_vartime(&sumBRA, &R, &A, tau, rho, S1, S2);

    return ge25519_is_neutral_vartime(&sumBRA) ? 0: -1;
}

/* New single verification using QSM_B with hEEA */
int
ED25519_FN(ed25519_sign_open_hEEA_samePre) (const unsigned char *m, size_t mlen, const ed25519_public_key pk, const ed25519_signature RS) {
	ge25519 ALIGN(16) R, A, sumBRA;
	hash_512bits hash;
	bignum256modm hram, S1, S2={0}, rho, tau;
    int tau_isNegative;
	int rho_isNegative;

	if ((RS[63] & 224))
		return -1;

	/* hram <-- H(R,A,m) */
	ed25519_hram(hash, RS, pk, m, mlen);
	expand256_modm(hram, hash, 64);

	/* compute tau and rho s.t. tau * h = rho mod el */
    curve25519_half_size_scalar_vartime_hEEA(rho, tau, hram, &rho_isNegative, &tau_isNegative);
    
    /* unpacking (-R) */
    if (!ge25519_unpack_negative_vartime(&R, RS))
		return -1;

	/*
		[S]B + (-R) + h(-A) =? 0 ==> [tau * S]B + [tau](-R) + [rho](-A) =? 0
		----------------------------------------------------------------------------------
		tau_isNegative | rho_isNegative | [tau * s]B + [tau](-R) + [rho](-A) =? 0
		----------------------------------------------------------------------------------
		       0   	   |        0       | [|tau| * s]B + [|tau|](-R) + [|rho|](-A)
		----------------------------------------------------------------------------------
		       0   	   |        1       | [|tau| * s]B + [|tau|](-R) - [|rho|](-A) ==>
		               |                | [|tau| * s]B + [|tau|](-R) + [|rho|](A)
		----------------------------------------------------------------------------------
		       1       |        0       | [-|tau| * s]B + [-|tau|](-R) + [|rho|](-A) ==> 
		               |                | [|tau| * s] B + [|tau|](-R) + [|rho|](A)
		----------------------------------------------------------------------------------
		       1       |        1       | [-|tau| * s]B + [-|tau|](-R) + [-|rho|](-A) ==> 
		               |                | [|tau| * s] B + [|tau|](-R) + [|rho|](-A)
		----------------------------------------------------------------------------------
	*/
	/* unpacking (-A) and adjust the sign based on the previous table */
    if (tau_isNegative == rho_isNegative){
		if (!ge25519_unpack_negative_vartime(&A, pk))
			return -1;
	}else{
		if (!ge25519_unpack_positive_vartime(&A, pk))
			return -1;
	}
    
    /* S */
	expand256_modm(S1, RS + 32, 32);

    /* S <-- rS */
    mul256_modm(S1, S1, tau);
    
    /* split S to S1 and S2, s.t. S = (S2<<126 + S1), bl(S1) = 126, bl(S2) = 127 */
    S2[0] = (S1[2] >> 14) | ((S1[3] & 0x3FFF) << 42);
    S2[1] = (S1[3] >> 14) | ((S1[4] & 0x3FFF) << 42);
    S2[2] = S1[4] >> 14;
    S1[2] &= 0x3FFF;
    S1[3] = 0;
    S1[4] = 0;

	/* 
	 * [S]B + (-R) + h(-A) =? 0 ==>
	 * [tau * S]B + [tau](-R) + [rho](-A) =? 0 ==>
	 * [S1]B + [S2]([2^126]B) + [tau](-R) + [rho](-A) =? 0 
	*/
	// ge25519_quadruple_scalarmult_vartime(&sumBRA, &R, &A, tau, rho, S1, S2);
	ge25519_quadruple_scalarmult_samePre_vartime(&sumBRA, &R, &A, tau, rho, S1, S2);

    return ge25519_is_neutral_vartime(&sumBRA) ? 0: -1;
}


/* New single verification using QSM_B_B' with hgcd_enhance2 */
int
ED25519_FN(ed25519_sign_open_hgcd) (const unsigned char *m, size_t mlen, const ed25519_public_key pk, const ed25519_signature RS) {
	ge25519 ALIGN(16) R, A, sumBRA;
	hash_512bits hash;
	bignum256modm hram, S1, S2={0}, rho, tau;
    int tau_isNegative;
	int rho_isNegative;

	if ((RS[63] & 224))
		return -1;

	/* hram <-- H(R,A,m) */
	ed25519_hram(hash, RS, pk, m, mlen);
	expand256_modm(hram, hash, 64);

	/* compute tau and rho s.t. tau * h = rho mod el */
    // curve25519_half_size_scalar_vartime(rho, tau, hram, &rho_isNegative, &tau_isNegative);
    curve25519_half_size_scalar_vartime_hgcd(rho, tau, hram, &rho_isNegative, &tau_isNegative);
    
    /* unpacking (-R) */
    if (!ge25519_unpack_negative_vartime(&R, RS))
		return -1;

	/*
		[S]B + (-R) + h(-A) =? 0 ==> [tau * S]B + [tau](-R) + [rho](-A) =? 0
		----------------------------------------------------------------------------------
		tau_isNegative | rho_isNegative | [tau * s]B + [tau](-R) + [rho](-A) =? 0
		----------------------------------------------------------------------------------
		       0   	   |        0       | [|tau| * s]B + [|tau|](-R) + [|rho|](-A)
		----------------------------------------------------------------------------------
		       0   	   |        1       | [|tau| * s]B + [|tau|](-R) - [|rho|](-A) ==>
		               |                | [|tau| * s]B + [|tau|](-R) + [|rho|](A)
		----------------------------------------------------------------------------------
		       1       |        0       | [-|tau| * s]B + [-|tau|](-R) + [|rho|](-A) ==> 
		               |                | [|tau| * s] B + [|tau|](-R) + [|rho|](A)
		----------------------------------------------------------------------------------
		       1       |        1       | [-|tau| * s]B + [-|tau|](-R) + [-|rho|](-A) ==> 
		               |                | [|tau| * s] B + [|tau|](-R) + [|rho|](-A)
		----------------------------------------------------------------------------------
	*/
	/* unpacking (-A) and adjust the sign based on the previous table */
    if (tau_isNegative == rho_isNegative){
		if (!ge25519_unpack_negative_vartime(&A, pk))
			return -1;
	}else{
		if (!ge25519_unpack_positive_vartime(&A, pk))
			return -1;
	}
    
    /* S */
	expand256_modm(S1, RS + 32, 32);

    /* S <-- tau * S */
    mul256_modm(S1, S1, tau);
    
    /* split S to S1 and S2, s.t. S = (S2<<126 + S1), bl(S1) = 126, bl(S2) = 127 */
    S2[0] = (S1[2] >> 14) | ((S1[3] & 0x3FFF) << 42);
    S2[1] = (S1[3] >> 14) | ((S1[4] & 0x3FFF) << 42);
    S2[2] = S1[4] >> 14;
    S1[2] &= 0x3FFF;
    S1[3] = 0;
    S1[4] = 0;

	/* 
	 * [S]B + (-R) + h(-A) =? 0 ==>
	 * [tau * S]B + [tau](-R) + [rho](-A) =? 0 ==>
	 * [S1]B + [S2]([2^126]B) + [tau](-R) + [rho](-A) =? 0 
	*/
	ge25519_quadruple_scalarmult_vartime(&sumBRA, &R, &A, tau, rho, S1, S2);

    return ge25519_is_neutral_vartime(&sumBRA) ? 0: -1;
}

/* Old single verification using DSM_B_doublePre in which the precomputed table for B is doubled */
int
ED25519_FN(ed25519_sign_open_2B) (const unsigned char *m, size_t mlen, const ed25519_public_key pk, const ed25519_signature RS) {
	ge25519 ALIGN(16) R, A;
	hash_512bits hash;
	bignum256modm hram, S;
	unsigned char checkR[32];

	if ((RS[63] & 224) || !ge25519_unpack_negative_vartime(&A, pk))
		return -1;

	/* hram = H(R,A,m) */
	ed25519_hram(hash, RS, pk, m, mlen);
	expand256_modm(hram, hash, 64);

	/* S */
	expand256_modm(S, RS + 32, 32);

	/* SB - H(R,A,m)A */
	ge25519_double_scalarmult_vartime2B(&R, &A, hram, S);
	ge25519_pack(checkR, &R);

	/* check that R = SB - H(R,A,m)A */
	return ed25519_verify(RS, checkR, 32) ? 0 : -1;
}