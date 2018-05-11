#include "algo.h"

#include <mpi.h>
#include "util.h"

void DLEQ_prove(fmpz_t g1, fmpz_t h1, fmpz_t g2, fmpz_t h2, fmpz_t alpha, fmpz_t c, fmpz_t w, fmpz_t r)
{
	fmpz_t a1,a2;
	fmpz_init(a1);
	fmpz_init(a2);
	
	fmpz_powm(a1,g1,w,g_groupOrder);
	fmpz_powm(a2,g2,w,g_groupOrder);
	fmpz_mul(r,alpha,c);
	fmpz_mod(r,r,g_groupOrder);
	fmpz_sub(r,g_groupOrder,r);
	fmpz_add(r,w,r);
	fmpz_mod(r,r,g_groupOrder);
}

int DLEQ_verify(fmpz_t g1, fmpz_t h1, fmpz_t g2, fmpz_t h2, fmpz_t c, fmpz_t a1, fmpz_t a2, fmpz_t r)
{
	fmpz_t t1,t2;
	fmpz_init(t1);
	fmpz_init(t2);

	fmpz_powm(t1,g1,r,g_groupOrder);
	fmpz_powm(t2,h1,c,g_groupOrder);
	fmpz_mul(t1,t1,t2);
	fmpz_mod(t1,t1,g_groupOrder);
	if (fmpz_cmp(t1,a1) != 0) {
		return -1;
	}

	fmpz_powm(t1,g2,r,g_groupOrder);
	fmpz_powm(t2,h2,c,g_groupOrder);
	fmpz_mul(t1,t1,t2);
	fmpz_mod(t1,t1,g_groupOrder);
	if (fmpz_cmp(t1,a2) != 0) {
		return -1;
	}
	
	return 0;
}

void sharedRandGenerator(fmpz_t x)
{
	fmpz_t y;
	fmpz_init(y);

	while (1) {
		if (g_rank == 1) {
			printf("loop\n");
		}
		sharedRand(x);
		fmpz_sub_ui(y,g_groupOrder,1);
		fmpz_divexact_ui(y,y,2);
		fmpz_powm(x,x,y,g_groupOrder);
		if (!fmpz_cmp_ui(x,1)) {
			break;
		}
	}
	
	fmpz_clear(y);
}

void sharedRand(fmpz_t x)
{
	fmpz_t y,gX,two;
	void *sendbuf, *recvbuf;
	int i;
	fmpz_t *commitments;
	
	fmpz_init(y);
	fmpz_init(gX);
	fmpz_init_set_ui(two,2);
	randGroupElt(x);
	fmpz_powm(gX,two,x,g_groupOrder);

	sendbuf=malloc(GROUP_ELT_SIZE);
	recvbuf=malloc(GROUP_ELT_SIZE*g_np);
	commitments=malloc(sizeof(fmpz_t)*g_np);

	exportGroupElt(gX,sendbuf);

	MPI_Allgather(sendbuf,GROUP_ELT_SIZE,MPI_CHAR,
	              recvbuf,GROUP_ELT_SIZE,MPI_CHAR,
	              MPI_COMM_WORLD);

	for (i=0;i<g_np;++i) {
		fmpz_init(commitments[i]);
		importGroupElt(commitments[i],recvbuf+GROUP_ELT_SIZE*i);
	}

	exportGroupElt(x,sendbuf);

	MPI_Allgather(sendbuf,GROUP_ELT_SIZE,MPI_CHAR,
	              recvbuf,GROUP_ELT_SIZE,MPI_CHAR,
	              MPI_COMM_WORLD);
	
	for (i=0;i<g_np;++i) {
		if (i != g_rank) {
			importGroupElt(y,recvbuf+GROUP_ELT_SIZE*i);
			fmpz_powm(gX,two,y,g_groupOrder);
			if (fmpz_cmp(gX,commitments[i]) != 0) {
				fprintf(stderr,"Error, commitments do not agree\n");
				exit(-1);
			}
			fmpz_add(x,x,y);
			fmpz_mod(x,x,g_groupOrder);
		}
	}

	fmpz_clear(y);
	fmpz_clear(gX);
	fmpz_clear(two);
	for (i=0;i<g_np;++i) {
		fmpz_clear(commitments[i]);
	}
	free(commitments);
	free(recvbuf);
	free(sendbuf);
}

