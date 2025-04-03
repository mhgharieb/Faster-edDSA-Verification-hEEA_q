#ifndef ED25519_H
#define ED25519_H

#include <stdlib.h>

#if defined(__cplusplus)
extern "C" {
#endif

typedef unsigned char ed25519_signature[64];
typedef unsigned char ed25519_public_key[32];
typedef unsigned char ed25519_secret_key[32];

typedef unsigned char curved25519_key[32];

void ed25519_publickey(const ed25519_secret_key sk, ed25519_public_key pk);
int ed25519_sign_open(const unsigned char *m, size_t mlen, const ed25519_public_key pk, const ed25519_signature RS);
void ed25519_sign(const unsigned char *m, size_t mlen, const ed25519_secret_key sk, const ed25519_public_key pk, ed25519_signature RS);
int ed25519_sign_open_batch(const unsigned char **m, size_t *mlen, const unsigned char **pk, const unsigned char **RS, size_t num, int *valid);

/* New single verification using QSM_B_B' with hEEA */
int ed25519_sign_open_hEEA(const unsigned char *m, size_t mlen, const ed25519_public_key pk, const ed25519_signature RS);
/* New single verification using QSM_B with hEEA */
int ed25519_sign_open_hEEA_samePre(const unsigned char *m, size_t mlen, const ed25519_public_key pk, const ed25519_signature RS);
/* New single verification using QSM_B_B' with hgcd_enhance2 */
int ed25519_sign_open_hgcd(const unsigned char *m, size_t mlen, const ed25519_public_key pk, const ed25519_signature RS);
/* Old single verification using DSM_B_doublePre in which the precomputed table for B is doubled */
int ed25519_sign_open_2B(const unsigned char *m, size_t mlen, const ed25519_public_key pk, const ed25519_signature RS);

/* New batch verification using the new randamization with hEEA */
int ed25519_sign_open_batch_hEEA(const unsigned char **m, size_t *mlen, const unsigned char **pk, const unsigned char **RS, size_t num, int *valid);
/* New batch verification using the new randamization with hgcd_enhance2 */
int ed25519_sign_open_batch_hgcd(const unsigned char **m, size_t *mlen, const unsigned char **pk, const unsigned char **RS, size_t num, int *valid);

void ed25519_randombytes_unsafe(void *out, size_t count);

void curved25519_scalarmult_basepoint(curved25519_key pk, const curved25519_key e);

#if defined(__cplusplus)
}
#endif

#endif // ED25519_H
