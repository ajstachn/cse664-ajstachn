#ifndef _H_UTIL
#define _H_UTIL

#include "prime.h"

#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"
#include <openssl/conf.h>
#include <openssl/evp.h>
#include <openssl/err.h>
#include <openssl/pem.h>
#include <openssl/rand.h>
#include <openssl/engine.h>

#include <sys/random.h>

#define SEED_SIZE 256
#define GROUP_ELT_SIZE 1028

extern int g_rank, g_np;
extern fmpz_t g_groupOrder,g_groupMod, g_generator1, g_generator2;

void init(int argc, char** argv);
void cleanup();
void randInt(fmpz_t n, fmpz_t maxVal, int maxValLogCeil);
void seedRNG();
void print_hex_memory(void *mem);
void randModQ(fmpz_t n);
void randUnitModQ(fmpz_t n);
void randGroupElt(fmpz_t n);
void evalPoly(fmpz_t y,fmpz_poly_t p, ulong x);
void randPoly(fmpz_poly_t p, int degree);
void lagrangeCoeff(fmpz_t c,int n,int j);
void hashElts(fmpz_t result,fmpz* elts,int n);
void hashToFmpz(fmpz_t n,char* hash,int md_len);
void exportGroupElt(fmpz_t n,char* buf);
void importGroupElt(fmpz_t n,char* buf);
void read_private_key(fmpz_t key);
void read_public_key(fmpz_t key,int pid);
void initKeys(fmpz_t privateKey,fmpz **pubKeys);
void freeKeys(fmpz_t privateKey,fmpz *pubKeys);
void bcast_fmpz(fmpz_t n, int root);
void bcast_n_fmpz(fmpz *arr, int n, int root);
void bcast_all_n_fmpz(fmpz *arr, int n);
void bcast_fmpz_vec(fmpz *arr);
void setFmpz(fmpz_t a, const fmpz_t b);

unsigned int debugHash(fmpz_t n);

void init_fmpz_array(fmpz **arr, int n);
void free_fmpz_array(fmpz *arr, int n);

#endif

