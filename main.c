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

	fprintf(this.log,"Polynomial:\n");
	for (i=0;i<g_np;++i) {
		fmpz_fprint(this.log,fmpz_poly_get_coeff_ptr(this.myPoly,i));
		fprintf(this.log,"\n\n");
	}

	fprintf(this.log,"Key:\n");
	fmpz_fprint(this.log,this.privateKey);
	fprintf(this.log,"\n\n");

	generateChallenge(&this);
	computeDLEQ(&this);

	logShared(&this);
	publishCYUR(&this);
	if (g_rank == 0) {
		printf("Secret shares distributed\n");
	}

	proveVoteSoundness(&this);
	publishVoteSoundnessProof(&this);
	verifyVoteSoundness(&this);
	if (g_rank == 0) {
		printf("Vote soundness verified\n");
	}

	logShared(&this);

	verifySecret(&this);
	if (g_rank == 0) {
		printf("Shared secrets verified\n");
	}

	reducePolyPoints(&this);
	computeShare(&this);
	
	if (g_rank == 0) {
		printf("Shares decrypted\n");
	}
	
	publishShares(&this);
	verifyShares(&this);
	
	if (g_rank == 0) {
		printf("Decrypted votes verified\n");
	}

	reconstructSecret(&this);
	tally=tallyVotes(&this);

	if (g_rank == 0) {
		printf("The vote is Yay:%d Nay:%d\n",tally,g_np-tally);
	}

	freeVars(&this);
	cleanup();
	return 0;
}
