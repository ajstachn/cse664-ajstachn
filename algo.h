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

struct LocalVars {
	int myVote;
	fmpz_t privateKey;
	fmpz *pubKeys;
	fmpz_poly_t myPoly;

	fmpz *uVotes;
	fmpz *myPolyPoints;
	fmpz *myPolyPointsExp;
	fmpz *allPolyCoeffExp;
	fmpz *allPolyPointsEnc;
	fmpz *challenges;

	fmpz *wVec;

	fmpz *allResponses;

	fmpz *shares;
	fmpz *shareProofs;
	fmpz *shareA1s, *shareA2s;

	fmpz_t sharedSecret;

	FILE *log;
};

void initVars(struct LocalVars *this);
void freeVars(struct LocalVars *this);
void computeCYU(struct LocalVars *this);
void publishCYUR(struct LocalVars *this);
void generateChallenge(struct LocalVars *this);
void computeDLEQ(struct LocalVars *this);
void verifySecret(struct LocalVars *this);
void DLEQ_prove(fmpz_t g1, fmpz_t h1, fmpz_t g2, fmpz_t h2, fmpz_t alpha, fmpz_t c, fmpz_t w, fmpz_t r);
void newRand(struct LocalVars *this);
void computeShare(struct LocalVars *this);
void publishShares(struct LocalVars *this);
void verifyShares(struct LocalVars *this);
void reconstructSecret(struct LocalVars *this);
int tallyVotes(struct LocalVars *this);

int discreteLog(fmpz_t g,fmpz_t x);

void logShared(struct LocalVars *this);
void logHashArray(struct LocalVars *this,fmpz *hashArray);

void sharedRandGenerator(fmpz_t x);
void sharedRand(fmpz_t x);

#endif
