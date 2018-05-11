#include "util.h"

#include <stdio.h>
#include <mpi.h>
#include <sys/random.h>

int g_rank, g_np;
fmpz_t g_groupOrder,g_generator1,g_generator2;

void init(int argc, char** argv)
{
	MPI_Init(&argc,&argv);

	ERR_load_crypto_strings();
	OpenSSL_add_all_algorithms();
	OPENSSL_config(NULL);
	seedRNG();

	MPI_Comm_rank(MPI_COMM_WORLD,&g_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&g_np);

	fmpz_init(g_groupOrder);
	fmpz_set_str(g_groupOrder,PRIME_GROUP_ORDER,16);

	fmpz_init(g_generator1);
	fmpz_set_str(g_generator1,GROUP_GENERATOR1,16);

	fmpz_init(g_generator2);
	fmpz_set_str(g_generator2,GROUP_GENERATOR2,16);
}

void read_private_key(fmpz_t key)
{
	char keyname[128];
	FILE *file;
	
	snprintf(keyname,128,"keys/private_key_%d",g_rank);
	file=fopen(keyname,"r");
	fmpz_inp_raw(key,file);
	fclose(file);
}

void read_public_key(fmpz_t key,int pid)
{
	char keyname[128];
	FILE *file;
	
	snprintf(keyname,128,"keys/public_key_%d",pid);
	file=fopen(keyname,"r");
	fmpz_inp_raw(key,file);
	fclose(file);
}

void cleanup()
{
	fmpz_clear(g_groupOrder);
	fmpz_clear(g_generator1);
	fmpz_clear(g_generator2);

	EVP_cleanup();
	CRYPTO_cleanup_all_ex_data();
	ERR_free_strings();
	
	MPI_Finalize();
}

void exportGroupElt(fmpz_t n,char* buf)
{
	FILE *file=tmpfile();
	fmpz_out_raw(file,n);
	rewind(file);
	fread(buf,1,GROUP_ELT_SIZE,file);
	fclose(file);
}

void importGroupElt(fmpz_t n,char* buf)
{
	FILE *file=tmpfile();
	fwrite(buf,1,GROUP_ELT_SIZE,file);
	rewind(file);
	fmpz_inp_raw(n,file);
	fclose(file);
}

void randGroupElt(fmpz_t n)
{
	randInt(n,g_groupOrder,1024);
}

void evalPoly(fmpz_t y,fmpz_poly_t p, ulong x)
{
	int i,d;
	fmpz *c;
	
	d=fmpz_poly_degree(p);

	fmpz_poly_get_coeff_fmpz(y,p,d);
	for (i=d-1;i>=0;--i) {
		fmpz_mul_ui(y,y,x);
		c=fmpz_poly_get_coeff_ptr(p,i);
		fmpz_add(y,y,c);
		fmpz_mod(y,y,g_groupOrder);
	}
}

void hashElts(fmpz_t result, fmpz *elts, int n)
{
	EVP_MD_CTX *mdctx;
	int i,md_len;
	char *buf,*digest;
	mdctx=EVP_MD_CTX_create();
	if (mdctx == NULL) {
		printf("Error creating context\n");
		exit(-1);
	}
	if (EVP_DigestInit_ex(mdctx,EVP_sha256(),NULL) != 1) {
		printf("Error in DigestInit\n");
		exit(-1);
	}
	buf=malloc(GROUP_ELT_SIZE);
	for (i=0;i<n;++i) {
		exportGroupElt(elts+i,buf);
		EVP_DigestUpdate(mdctx,buf,GROUP_ELT_SIZE);
	}
	digest=malloc(EVP_MD_size(EVP_sha256()));
	EVP_DigestFinal_ex(mdctx,digest,&md_len);
	EVP_MD_CTX_destroy(mdctx);
	
	hashToGroupElt(result,digest,md_len);
	free(buf);
	free(digest);
}



// take uninitialized p
void randPoly(fmpz_poly_t p, int degree)
{
	int i;
	fmpz_t x;
	
	fmpz_init(x);
	fmpz_poly_init2(p,degree);
	
	for (i=0;i<degree;++i) {
		randGroupElt(x);
		fmpz_poly_set_coeff_fmpz(p,i,x);
	}

	fmpz_clear(x);
}

void lagrangeCoeff(fmpz_t c,int n,int i)
{
	int j;
	fmpz_t d;
	fmpz_init(d);
	fmpz_set_ui(c,1);
	for (j=1;j<=n;++j) {
		if (i != j) {
			fmpz_set_si(d,j-i);
			if (j-i < 0) {
				fmpz_add(d,d,g_groupOrder);
			}
			fmpz_invmod(d,d,g_groupOrder);
			fmpz_mul_ui(d,d,j);
			fmpz_mul(c,d,c);
		}
	}
	fmpz_clear(d);
}

const char hexLookup[]="0123456789ABCDEF";

void hashToGroupElt(fmpz_t n,char* hash,int md_len)
{
	char str[EVP_MAX_MD_SIZE*2+1];
	int i;
	for (i=0;i<md_len;++i) {
		str[2*i]=hexLookup[hash[i]&0xF];
		str[2*i+1]=hexLookup[(hash[i]&0xF0)>>4];
	}
	str[2*md_len]='\0';
	fmpz_set_str(n,str,16);
}

//Return n in [0,maxVal)
//maxValLogCeil must be >= log_256(maxVal)
void randInt(fmpz_t n, fmpz_t maxVal, int maxValLogCeil)
{
	char *bytes,*str;
	int i;
	bytes=malloc(sizeof(char)*maxValLogCeil);
	str=malloc(sizeof(char)*maxValLogCeil*2+1);
	while (1) {
		RAND_bytes(bytes,maxValLogCeil);
		for (i=0;i<maxValLogCeil;++i) {
			str[2*i]=hexLookup[bytes[i]&0xF];
			str[2*i+1]=hexLookup[(bytes[i]&0xF0)>>4];
		}
		str[2*maxValLogCeil]='\0';
		fmpz_set_str(n,str,16);
		if (fmpz_cmp(n,maxVal) < 0) {
			break;
		}
	}
	free(str);
	free(bytes);
}

void seedRNG()
{
	int result;
	char rand[SEED_SIZE];
	result=getrandom(rand,SEED_SIZE,0);
	if (result < SEED_SIZE) {
		perror("getrandom");
		exit(-1);
	}
	RAND_seed(rand,SEED_SIZE);
	result=RAND_status();
	if (result != 1) {
		fprintf(stderr,"RAND_seed error\n");
		exit(-1);
	}
}

void print_hex_memory(void *mem)
{
	int i;
	unsigned char *p = (unsigned char *)mem;
	for (i=0;i<512;i++) {
		printf("0x%02x ", p[i]);
		if ((i%16==0) && i)
			printf("\n");
	}
	printf("\n");
}

void bcast_fmpz(fmpz_t n, int root)
{
	char buf[GROUP_ELT_SIZE];
	if (g_rank == root) {
		exportGroupElt(n,buf);
	}
	MPI_Bcast(buf,GROUP_ELT_SIZE,MPI_BYTE,root,MPI_COMM_WORLD);
	if (g_rank != root) {
		importGroupElt(n,buf);
	}
}

void bcast_n_fmpz(fmpz *arr, int n, int root)
{
	int i;
	for (i=0;i<n;++i) {
		bcast_fmpz(arr+i,root);
	}
}

void bcast_all_n_fmpz(fmpz *arr, int n)
{
	int i;
	for (i=0;i<g_np;++i) {
		bcast_n_fmpz(arr+i*n,n,i);
	}
}

void init_fmpz_array(fmpz **arr, int n)
{
	int i;
	(*arr)=malloc(sizeof(fmpz)*n);
	for (i=0;i<n;++i) {
		fmpz_init((*arr)+i);
	}
}

void free_fmpz_array(fmpz *arr, int n)
{
	int i;
	for (i=0;i<n;++i) {
		fmpz_clear(arr+i);
	}
	free(arr);
}

int debugHash(fmpz_t n)
{
	int retVal;
	fmpz_t temp,temp2;
	fmpz_init(temp);
	fmpz_init(temp2);
	fmpz_add_ui(temp2,n,0);
	
	hashElts(temp,temp2,1);
	
	fmpz_mod_ui(temp,temp,4294967291);
	retVal=fmpz_get_ui(temp);
	fmpz_clear(temp2);
	fmpz_clear(temp);
	return retVal;
}
