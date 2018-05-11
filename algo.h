#ifndef _H_ALGO
#define _H_ALGO

#include "util.h"

#include "flint/fmpz.h"
#include <openssl/conf.h>
#include <openssl/evp.h>
#include <openssl/err.h>
#include <openssl/pem.h>
#include <openssl/rand.h>
#include <openssl/engine.h>

#include <sys/random.h>

void DLEQ_prove(fmpz_t g1, fmpz_t h1, fmpz_t g2, fmpz_t h2, fmpz_t alpha, fmpz_t c, fmpz_t w, fmpz_t r);
int DLEQ_verify(fmpz_t g1, fmpz_t h1, fmpz_t g2, fmpz_t h2, fmpz_t c, fmpz_t a1, fmpz_t a2, fmpz_t r);
                 
                 

void sharedRandGenerator(fmpz_t x);
void sharedRand(fmpz_t x);

#endif
