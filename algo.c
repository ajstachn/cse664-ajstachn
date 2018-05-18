#include "algo.h"

#include <mpi.h>
#include "util.h"

int discreteLog(fmpz_t g,fmpz_t x)
{
	fmpz_t temp;
	int i=0;

	fmpz_init(temp);
	while(1) {
		if (g_rank == 0) {
			printf("discretelog\n");
		}
		fmpz_set_ui(temp,i);
		fmpz_powm(temp,g,temp,g_groupMod);
		if (fmpz_cmp(temp,x) == 0) {
			break;
		}
	}
	return i;
}

int tallyVotes(struct LocalVars *this)
{
	int i;
	fmpz_t uStar,gv;

	fmpz_init(uStar);
	fmpz_init(gv);
	
	fmpz_set_ui(uStar,1);
	for (i=0;i<g_np;++i) {
		fmpz_mul(uStar,uStar,(this->uVotes)+i);
		fmpz_mod(uStar,uStar,g_groupMod);
	}

	fmpz_invmod(gv,this->sharedSecret,g_groupMod);
	fmpz_mul(gv,gv,uStar);
	fmpz_mod(gv,gv,g_groupMod);
	return discreteLog(g_generator2,gv);
}

void reconstructSecret(struct LocalVars *this)
{
	fmpz *yStar;
	fmpz_t lambda,temp;
	int voter,i;

	init_fmpz_array(&yStar,g_np);
	fmpz_init(lambda);
	fmpz_init(temp);
	
	for (i=0;i<g_np;++i) {
		fmpz_set_ui(yStar+i,1);
		for (voter=0;voter<g_np;++voter) {
			fmpz_mul(yStar+i,yStar+i,(this->shares)+(voter*g_np)+i);
			fmpz_mod(yStar+i,yStar+i,g_groupMod);
		}
	}

	fmpz_set_ui(this->sharedSecret,1);
	for (i=0;i<g_np;++i) {
		lagrangeCoeff(lambda,g_np,i);
		fmpz_powm(temp,yStar+i,lambda,g_groupMod);
		fmpz_mul(this->sharedSecret,this->sharedSecret,temp);
		fmpz_mod(this->sharedSecret,this->sharedSecret,g_groupMod);
	}
	fprintf(this->log,"Shared secret:\n");
	fmpz_fprint(this->log,this->sharedSecret);
	fprintf(this->log,"\n\n");

	free_fmpz_array(yStar,g_np);
	fmpz_clear(lambda);
	fmpz_clear(temp);
}

void verifyShares(struct LocalVars *this)
{
	int tallier,i;
	fmpz_t temp1,temp2,ri,ci;

	fmpz_init(temp1);
	fmpz_init(temp2);
	fmpz_init(ri);
	fmpz_init(ci);

	for (tallier=0;tallier<g_np;++tallier) {
		for (i=0;i<g_np;++i) {
			setFmpz(ri,(this->shareProofs)+(g_np*tallier)+i);
			setFmpz(ci,(this->challenges)+tallier);
			
			fmpz_powm(temp1,g_generator2,ri,g_groupMod);
			fmpz_powm(temp2,(this->pubKeys)+tallier,ci,g_groupMod);
			fmpz_mul(temp1,temp1,temp2);
			fmpz_mod(temp1,temp1,g_groupMod);
			fprintf(this->log,"%u %u\n",debugHash(temp1),debugHash((this->shareA1s)+(g_np*tallier)+i));
			if (fmpz_cmp(temp1,(this->shareA1s)+(g_np*tallier)+i) != 0) {
				printf("Error a1s do not match, Tallier: %d, Verifier %d\n",tallier,g_rank);
				cleanup();
				exit(0);
			}

			fmpz_powm(temp1,(this->shares)+(g_np*tallier)+i,ri,g_groupMod);
			fmpz_powm(temp2,(this->allPolyPointsEnc)+(i*g_np)+tallier,ci,g_groupMod);
			fmpz_mul(temp1,temp1,temp2);
			fmpz_mod(temp1,temp1,g_groupMod);
			fprintf(this->log,"%u %u\n",debugHash(temp1),debugHash((this->shareA2s)+(g_np*tallier)+i));
			if (fmpz_cmp(temp1,(this->shareA2s)+(g_np*tallier)+i) != 0) {
				printf("Error a2s do not match, Tallier: %d, Verifier %d\n",tallier,g_rank);
				cleanup();
				exit(0);
			}
		}
	}

	fmpz_clear(temp1);
	fmpz_clear(temp2);
	fmpz_clear(ri);
	fmpz_clear(ci);
}

void publishShares(struct LocalVars *this)
{
	bcast_all_n_fmpz(this->shares,g_np);
	bcast_all_n_fmpz(this->shareProofs,g_np);

	bcast_all_n_fmpz(this->shareA1s,g_np);
	bcast_all_n_fmpz(this->shareA2s,g_np);
}

void newRand(struct LocalVars *this)
{
	int i;
	for (i=0;i<g_np;++i) {
		randModQ((this->wVec)+i);
	}
}

void computeShare(struct LocalVars *this)
{
	int voter;
	fmpz_t x,pointEnc;

	fmpz_t y,z;
	fmpz_init(y);
	fmpz_init(z);
	
	fmpz_init(x);
	fmpz_init(pointEnc);
	for (voter=0;voter<g_np;++voter) {
		setFmpz(pointEnc,(this->allPolyPointsEnc)+(voter*g_np)+g_rank);
		fmpz_invmod(x,this->privateKey,g_groupOrder);
		fmpz_powm(y,pointEnc,x,g_groupMod);
		fmpz_powm(z,y,this->privateKey,g_groupMod);
		setFmpz((this->shares)+(g_np*g_rank)+voter,y);

		fmpz_powm((this->shareA1s)+(g_np*g_rank)+voter,g_generator2,(this->wVec)+voter,g_groupMod);
		fmpz_powm((this->shareA2s)+(g_np*g_rank)+voter,y,(this->wVec)+voter,g_groupMod);

		DLEQ_prove(g_generator2,
		           (this->pubKeys)+g_rank,
		           y,
		           pointEnc,
		           this->privateKey,
		           (this->challenges)+g_rank,
		           (this->wVec)+voter,
		           (this->shareProofs)+(g_np*g_rank)+voter);
	}
	fmpz_clear(x);
}

void computeDLEQ(struct LocalVars *this)
{
	int i;
	fmpz_t r;

	fmpz_init(r);

	for (i=0;i<g_np;++i) {
		DLEQ_prove(g_generator1,
		           (this->myPolyPointsExp)+i,
		           (this->pubKeys)+i,
		           (this->allPolyPointsEnc)+(g_np*g_rank)+i,
		           (this->myPolyPoints)+i,
		           (this->challenges)+g_rank,
		           (this->wVec)+i,
		           r);
		setFmpz((this->allResponses)+(g_np*g_rank)+i,r);
	}

	fmpz_clear(r);
}

void verifySecret(struct LocalVars *this)
{
	int i,j,voter;
	fmpz_t xi,yi,ri,challenge,pt,temp1,temp2,a1i,a2i,hash;
	fmpz *hashArr;

	fmpz_init(xi);
	fmpz_init(yi);
	fmpz_init(ri);
	fmpz_init(a1i);
	fmpz_init(a2i);
	fmpz_init(hash);
	fmpz_init(challenge);
	fmpz_init(pt);
	fmpz_init(temp1);
	fmpz_init(temp2);
	init_fmpz_array(&hashArr,g_np*4);

	for (voter=0;voter<g_np;++voter) {
		for (i=0;i<g_np;++i) {
			fmpz_set_ui(xi,1);
			for (j=0;j<g_np;++j) {
				fmpz_set_ui(pt,i+1);
				fmpz_pow_ui(temp1,pt,j);
				fmpz_powm(temp1,this->allPolyCoeffExp+(voter*g_np)+j,temp1,g_groupMod);
				fmpz_mul(xi,xi,temp1);
				fmpz_mod(xi,xi,g_groupMod);
			}
			setFmpz(yi,(this->allPolyPointsEnc)+(voter*g_np)+i);
			setFmpz(ri,(this->allResponses)+(voter*g_np)+i);
			setFmpz(challenge,(this->challenges)+voter);
			
			fmpz_powm(temp1,g_generator1,ri,g_groupMod);
			fmpz_powm(temp2,xi,challenge,g_groupMod);
			fmpz_mul(temp1,temp1,temp2);
			fmpz_mod(a1i,temp1,g_groupMod);

			fmpz_powm(temp1,(this->pubKeys)+i,ri,g_groupMod);
			fmpz_powm(temp2,yi,challenge,g_groupMod);
			fmpz_mul(temp1,temp1,temp2);
			fmpz_mod(a2i,temp1,g_groupMod);
			
			
			setFmpz(hashArr+(i*4),xi);
			setFmpz(hashArr+(i*4)+1,yi);
			setFmpz(hashArr+(i*4)+2,a1i);
			setFmpz(hashArr+(i*4)+3,a2i);
		}
		logHashArray(this,hashArr);
		hashElts(hash,hashArr,g_np*4);
		if (fmpz_cmp(hash,challenge) != 0) {
			printf("Error, hashes do not match, proc %d, challenge %d\n",g_rank,voter);
			cleanup();
			exit(0);
		}
	}

	free_fmpz_array(hashArr,g_np*4);
	fmpz_clear(xi);
	fmpz_clear(yi);
	fmpz_clear(ri);
	fmpz_clear(a1i);
	fmpz_clear(a2i);
	fmpz_clear(hash);
	fmpz_clear(challenge);
	fmpz_clear(pt);
	fmpz_clear(temp1);
	fmpz_clear(temp2);
}

void logHashArray(struct LocalVars *this,fmpz *hashArray)
{
	int i;
	fprintf(this->log,"BEGIN hash array\n");
	for (i=0;i<g_np;++i) {
		fprintf(this->log,"%u %u %u %u\n",
		        debugHash(hashArray+(i*4)),
		        debugHash(hashArray+(i*4)+1),
		        debugHash(hashArray+(i*4)+2),
		        debugHash(hashArray+(i*4)+3));
	}
	fprintf(this->log,"END hash array\n\n");
}

void logShared(struct LocalVars *this)
{
	int i,j;
	fprintf(this->log,"BEGIN allPolyCoeffExp\n");
	for(i=0;i<g_np;++i) {
		for (j=0;j<g_np;++j) {
			fprintf(this->log,"%u ",debugHash(this->allPolyCoeffExp+(i*g_np)+j));
		}
		fprintf(this->log,"\n");
	}
	fprintf(this->log,"END allPolyCoeffExp\n");

	fprintf(this->log,"BEGIN allPolyPointsEnc\n");
	for(i=0;i<g_np;++i) {
		for (j=0;j<g_np;++j) {
			fprintf(this->log,"%u ",debugHash(this->allPolyPointsEnc+(i*g_np)+j));
		}
		fprintf(this->log,"\n");
	}
	fprintf(this->log,"END allPolyPointsEnc\n");


	fprintf(this->log,"BEGIN challenges\n");
	for(i=0;i<g_np;++i) {
		fprintf(this->log,"%u ",debugHash((this->challenges)+i));
	}
	fprintf(this->log,"\nEND challenges\n");

	fprintf(this->log,"BEGIN responses\n");
	for(i=0;i<g_np;++i) {
		for (j=0;j<g_np;++j) {
			fprintf(this->log,"%u ",debugHash(this->allResponses+(i*g_np)+j));
		}
		fprintf(this->log,"\n");
	}
	fprintf(this->log,"END responses\n\n");
}

void initVars(struct LocalVars *this)
{
	int i;
	char path[128];
	initKeys(this->privateKey,&(this->pubKeys));
	randPoly(this->myPoly,g_np);

	fmpz_init(this->sharedSecret);

	init_fmpz_array(&(this->uVotes),g_np);
	
	init_fmpz_array(&(this->myPolyPoints),g_np);
	init_fmpz_array(&(this->myPolyPointsExp),g_np);
	init_fmpz_array(&(this->allPolyCoeffExp),g_np*g_np);
	init_fmpz_array(&(this->allPolyPointsEnc),g_np*g_np);

	init_fmpz_array(&(this->challenges),g_np);

	init_fmpz_array(&(this->wVec),g_np);
	for (i=0;i<g_np;++i) {
		randModQ((this->wVec)+i);
	}

	init_fmpz_array(&(this->allResponses),g_np*g_np);

	init_fmpz_array(&(this->shares),g_np*g_np);
	init_fmpz_array(&(this->shareProofs),g_np*g_np);
	init_fmpz_array(&(this->shareA1s),g_np*g_np);
	init_fmpz_array(&(this->shareA2s),g_np*g_np);

	snprintf(path,128,"logs/log_%d",g_rank);
	this->log=fopen(path,"w");
	fprintf(this->log,"BEGIN Log Rank: %d NP: %d\n",g_rank,g_np);
}

void freeVars(struct LocalVars *this)
{
	freeKeys(this->privateKey,this->pubKeys);
	fmpz_poly_clear(this->myPoly);

	fmpz_clear(this->sharedSecret);
	
	free_fmpz_array(this->uVotes,g_np);
	
	free_fmpz_array(this->myPolyPoints,g_np);
	free_fmpz_array(this->myPolyPointsExp,g_np);
	free_fmpz_array(this->allPolyCoeffExp,g_np*g_np);
	free_fmpz_array(this->allPolyPointsEnc,g_np*g_np);

	free_fmpz_array(this->challenges,g_np);

	free_fmpz_array(this->allResponses,g_np*g_np);

	free_fmpz_array(this->shares,g_np*g_np);
	free_fmpz_array(this->shareProofs,g_np*g_np);
	free_fmpz_array(this->shareA1s,g_np*g_np);
	free_fmpz_array(this->shareA2s,g_np*g_np);

	free_fmpz_array(this->wVec,g_np);
	fprintf(this->log,"END Log\n\n");
	fclose(this->log);
}

void computeCYU(struct LocalVars *this)
{
	fmpz_t x;
	int i;

	fmpz_init(x);
	for (i=0;i<g_np;++i) {
		fmpz_poly_get_coeff_fmpz(x,this->myPoly,i);
		fmpz_powm(x,g_generator1,x,g_groupMod);
		setFmpz((this->allPolyCoeffExp)+(g_np*g_rank)+i,x);

		evalPoly(x,this->myPoly,i+1);
		setFmpz((this->myPolyPoints)+i,x);
		fmpz_powm((this->myPolyPointsExp)+i,g_generator1,x,g_groupMod);
		fmpz_powm((this->allPolyPointsEnc)+(g_np*g_rank)+i,this->pubKeys+i,x,g_groupMod);
	}

	fmpz_poly_get_coeff_fmpz(x,this->myPoly,0);
	if (this->myVote == 1) {
		fmpz_add_ui(x,x,1);
	}
	fmpz_powm(this->uVotes+g_rank,g_generator2,x,g_groupMod);
	
	fmpz_clear(x);
}

void publishCYUR(struct LocalVars *this)
{
	bcast_all_n_fmpz(this->allPolyCoeffExp,g_np);
	bcast_all_n_fmpz(this->allPolyPointsEnc,g_np);
	bcast_fmpz_vec(this->challenges);
	bcast_fmpz_vec(this->uVotes);
	bcast_all_n_fmpz(this->allResponses,g_np);
}

void generateChallenge(struct LocalVars *this)
{
	fmpz *hashArr;
	fmpz_t x;
	int i;

	fmpz_init(x);
	init_fmpz_array(&hashArr,g_np*4);
	
	for (i=0;i<g_np;++i) {
		setFmpz(hashArr+(i*4),(this->myPolyPointsExp)+i);
		
		setFmpz(hashArr+(i*4)+1,(this->allPolyPointsEnc)+(g_rank*g_np)+i);
		
		fmpz_powm(x,g_generator1,(this->wVec)+i,g_groupMod);
		setFmpz(hashArr+(i*4)+2,x);

		fmpz_powm(x,(this->pubKeys)+i,(this->wVec)+i,g_groupMod);
		setFmpz(hashArr+(i*4)+3,x);
	}

	logHashArray(this,hashArr);
	hashElts((this->challenges)+g_rank,hashArr,g_np*4);

	free_fmpz_array(hashArr,g_np*4);
	fmpz_clear(x);
}

//Input: alpha=log_g1(h1)=log_g2(h2), challenge c, random w
//Output: proof r
void DLEQ_prove(fmpz_t g1, fmpz_t h1, fmpz_t g2, fmpz_t h2, fmpz_t alpha, fmpz_t c, fmpz_t w, fmpz_t r)
{
	fmpz_t a1,a2;
	fmpz_init(a1);
	fmpz_init(a2);
	
	fmpz_powm(a1,g1,w,g_groupMod);
	fmpz_powm(a2,g2,w,g_groupMod);
	fmpz_mul(r,alpha,c);
	fmpz_mod(r,r,g_groupOrder);
	fmpz_sub(r,g_groupOrder,r);
	fmpz_add(r,w,r);
	fmpz_mod(r,r,g_groupOrder);

	fmpz_clear(a1);
	fmpz_clear(a2);
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
		fmpz_sub_ui(y,g_groupMod,1);
		fmpz_divexact_ui(y,y,2);
		fmpz_powm(x,x,y,g_groupMod);
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
	fmpz_powm(gX,two,x,g_groupMod);

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
			fmpz_powm(gX,two,y,g_groupMod);
			if (fmpz_cmp(gX,commitments[i]) != 0) {
				fprintf(stderr,"Error, commitments do not agree\n");
				exit(-1);
			}
			fmpz_add(x,x,y);
			fmpz_mod(x,x,g_groupMod);
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

