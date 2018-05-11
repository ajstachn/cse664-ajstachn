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

void initKeys(fmpz_t privateKey, fmpz **pubKeys)
{
	int i;
	fmpz_init(privateKey);
	read_private_key(privateKey);

	*pubKeys=malloc(sizeof(fmpz)*g_np);
	for (i=0;i<g_np;++i) {
		read_public_key((*pubKeys)+i,i);
	}
}

void freeKeys(fmpz_t privateKey, fmpz *pubKeys)
{
	int i;
	fmpz_clear(privateKey);
	for (i=0;i<g_np;++i) {
		fmpz_clear(pubKeys+i);
	}
	free(pubKeys);
}

void secret_consistency_proof(fmpz *proofArr,fmpz *xArr, fmpz *pubKeys,
                              fmpz *pointVals,fmpz *pointShares, fmpz_t hash)
{
	int i;
	fmpz_t c;
	fmpz *hashInVals;
	fmpz *ws;
	
	fmpz_init(c);
	init_fmpz_array(&hashInVals,g_np*4);
	init_fmpz_array(&ws,g_np);
	for (i=0;i<g_np;++i) {
		fmpz_add_ui(hashInVals+(4*i),xArr+i,0);
		fmpz_add_ui(hashInVals+(4*i+1),pointShares+i,0);
		randGroupElt(ws+i);
		fmpz_powm(hashInVals+(4*i+2),g_generator1,ws+i,g_groupOrder);
		fmpz_powm(hashInVals+(4*i+3),pubKeys+i,ws+i,g_groupOrder);
	}
	if (g_rank == 0) {
		for (i=0;i<g_np;++i) {
			printf("%u %u %u %u\n",
			       debugHash(hashInVals+(4*i)),
			       debugHash(hashInVals+(4*i+1)),
			       debugHash(hashInVals+(4*i+2)),
			       debugHash(hashInVals+(4*i+3)));
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	hashElts(hash,hashInVals,g_np*4);
	for (i=0;i<g_np;++i) {
		DLEQ_prove(g_generator1,xArr+i,pubKeys+i,pointShares+i+g_np*g_rank,
		           pointVals+i,hash,ws+i,
		           proofArr+i+g_np*g_rank);
	}

	free_fmpz_array(ws,g_np);
	free_fmpz_array(hashInVals,g_np*4);
	fmpz_clear(c);
}

int verify_secret_consistency(fmpz *proofArr, fmpz *polyCommitments, fmpz *pubKeys, fmpz *pointShares, fmpz *hashArray)
{
	int i,j,pid;
	fmpz_t t1,t2,Xi;
	fmpz_t hash;
	fmpz *hashInVals;
	init_fmpz_array(&hashInVals,g_np*4);
	fmpz_init(t1);
	fmpz_init(t2);
	fmpz_init(Xi);
	fmpz_init(hash);

	for (pid=0;pid<g_np;++pid) {
		for (i=0;i<g_np;++i) {
			fmpz_set_ui(Xi,1);
			for (j=0;j<g_np;++j) {
				fmpz_set_ui(t2,i+1);
				fmpz_pow_ui(t2,t2,j);
				if (g_rank == 1) {
					fmpz_print(t2);
					printf("\n");
				}
				fmpz_powm(t1,polyCommitments+(g_np*pid+j),t2,g_groupOrder);
				fmpz_mul(Xi,Xi,t1);
			}
			fmpz_add_ui(hashInVals+4*i,Xi,0);
			if (g_rank == 1) {
				printf("\ninit: %u\n",debugHash(hashInVals+4*i));
			}
			if (g_rank == 1) {
				printf("init1.5: %u\n",debugHash(hashInVals+4*i));
			}
			if (g_rank == 1) {
				printf("init1.7: %u\n",debugHash(hashInVals+4*i));
			}
			fmpz_add_ui(hashInVals+(4*i+1),pointShares+i,0);
			fmpz_powm(t1,g_generator1,proofArr+(pid*g_np+i),g_groupOrder);
			fmpz_powm(t2,Xi,hashArray+pid,g_groupOrder);
			fmpz_mul(hashInVals+(4*i+2),t1,t2);
			fmpz_mod(hashInVals+(4*i+2),hashInVals+(4*i+2),g_groupOrder);
			fmpz_powm(t1,pubKeys+i,proofArr+(pid*g_np+i),g_groupOrder);
			fmpz_powm(t2,pointShares+(i+g_np*pid),hashArray+pid,g_groupOrder);
			fmpz_mul(hashInVals+(4*i+3),t1,t2);
			fmpz_mod(hashInVals+(4*i+3),hashInVals+(4*i+3),g_groupOrder);
			if (g_rank == 1) {
				printf("init2: %u\n",debugHash(hashInVals+(4*i)));
			}
		}
		if (g_rank == 1) {
			for (i=0;i<g_np;++i) {
				printf("final: %u %u %u %u\n",
				       debugHash(hashInVals+(4*i)),
				       debugHash(hashInVals+(4*i+1)),
				       debugHash(hashInVals+(4*i+2)),
				       debugHash(hashInVals+(4*i+3)));
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
						


		hashElts(hash,hashInVals,g_np*4);
		
		if (fmpz_cmp(hash,hashArray+pid) != 0) {
			printf("Proc %d: Hashes from %d do not match\n",g_rank,pid);
			fmpz_print(hash);
			printf("\n");
			fmpz_print(hashArray+pid);
			printf("\n");
			exit(-1);
		} else {
			printf("Proc %d: Hashes from %d match\n",g_rank,pid);
		}

	}

	fmpz_clear(hash);
	fmpz_clear(Xi);
	fmpz_clear(t2);
	fmpz_clear(t1);
	free_fmpz_array(hashInVals,g_np*4);
}

int main(int argc, char** argv)
{
	fmpz_t privateKey;
	fmpz *pubKeys;


	
	init(argc,argv);
	initKeys(privateKey,&pubKeys);

	int i;
	fmpz_poly_t p;
	randPoly(p,g_np);
	for (i=0;i<g_np;++i) {
		fmpz_poly_set_coeff_ui(p,i,i+1);
	}

	fmpz *polyCommitments, *pointShares, *xArr, *pointVals;
	fmpz *proofArr, *hashArr;
	init_fmpz_array(&polyCommitments,g_np*g_np);
	init_fmpz_array(&pointShares,g_np*g_np);
	init_fmpz_array(&xArr,g_np);
	init_fmpz_array(&proofArr,g_np*g_np);
	init_fmpz_array(&pointVals,g_np*g_np);
	init_fmpz_array(&hashArr,g_np);

	for (i=0;i<g_np;++i) {
		int k=i+(g_rank*g_np);
		fmpz_poly_get_coeff_fmpz(polyCommitments+k,p,i);
		fmpz_powm(polyCommitments+k,g_generator1,polyCommitments+k,g_groupOrder);


		evalPoly(pointVals+i,p,i+1);
		if (g_rank == 0) {
			fmpz_print(pointVals+k);
			printf(" ");
		}
		fmpz_powm(xArr+i,g_generator1,pointVals+i,g_groupOrder);
		fmpz_powm(pointShares+i,pubKeys+i,pointVals+i,g_groupOrder);
	}
	if (g_rank == 0) {
		printf("\n");
	}
	secret_consistency_proof(proofArr,xArr,pubKeys,pointVals,pointShares,hashArr+g_rank);
	bcast_all_n_fmpz(hashArr,1);
	bcast_all_n_fmpz(proofArr,g_np);
	bcast_all_n_fmpz(polyCommitments,g_np);
	bcast_all_n_fmpz(pointShares,g_np);

	verify_secret_consistency(proofArr,polyCommitments,pubKeys,pointShares,hashArr);

	free_fmpz_array(hashArr,g_np);
	free_fmpz_array(pointVals,g_np*g_np);
	free_fmpz_array(proofArr,g_np*g_np);
	free_fmpz_array(xArr,g_np);
	free_fmpz_array(pointShares,g_np*g_np);
	free_fmpz_array(polyCommitments,g_np*g_np);

	fmpz_poly_clear(p);

	freeKeys(privateKey,pubKeys);
	cleanup();
	return 0;
}
