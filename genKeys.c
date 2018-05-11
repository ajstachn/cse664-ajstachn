
#include "util.h"
#include <stdio.h>

int main(int argc, char** argv)
{
	fmpz_t x,y;
	char keyname[128];
	FILE *file;
	init(argc,argv);

	fmpz_init(x);
	fmpz_init(y);

	randGroupElt(x);
	fmpz_powm(y,g_generator1,x,g_groupOrder);

	snprintf(keyname,128,"keys/private_key_%d",g_rank);
	file=fopen(keyname,"w");
	fmpz_out_raw(file,x);
	fclose(file);

	snprintf(keyname,128,"keys/public_key_%d",g_rank);
	file=fopen(keyname,"w");
	fmpz_out_raw(file,y);
	fclose(file);

	fmpz_clear(x);
	fmpz_clear(y);
	
	cleanup();
	return 0;
}
