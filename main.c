#include "prime.h"
#include "util.h"
#include "algo.h"

#include <stdio.h>
#include <mpi.h>
#include <string.h>

#include "flint/fmpz.h"

#include <openssl/conf.h>
#include <openssl/evp.h>
#include <openssl/err.h>
#include <openssl/pem.h>
#include <openssl/rand.h>
#include <openssl/engine.h>



const char votes[]={1,0,1,0,1};


int main(int argc, char** argv)
{
	int i,j,tally;
	struct LocalVars this;
	
	init(argc,argv);
	initVars(&this);
	this.myVote=votes[g_rank];
	computeCYU(&this);

	fprintf(this.log,"Secret:\n");
	fmpz_fprint(this.log,fmpz_poly_get_coeff_ptr(this.myPoly,0));
	fprintf(this.log,"\n\n");

	generateChallenge(&this);
	computeDLEQ(&this);

	logShared(&this);
	publishCYUR(&this);

	logShared(&this);

	verifySecret(&this);
	if (g_rank == 0) {
		printf("Shared secret established and verified\n");
	}

	newRand(&this);

	computeShare(&this);
	publishShares(&this);
	verifyShares(&this);
	if (g_rank == 0) {
		printf("Decrypted votes verified\n");
	}

	reconstructSecret(&this);
	//tally=tallyVotes(&this);
	//printf("%d\n");

	freeVars(&this);
	cleanup();
	return 0;
}
