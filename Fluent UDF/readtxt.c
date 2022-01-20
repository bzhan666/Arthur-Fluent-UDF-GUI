#include "udf.h"

DEFINE_ON_DEMAND(read_file)
{
	int i = 0;
	int npost;
	FILE *fpin;
	int a[10];
	int b[10];
	int c[10];
	float d[10];
	float e[10];
	float f[10];
	float g[10];
	float x[10];
	float y[10];

	if ((fpin = fopen("input.txt", "r")) == NULL)
		{
			Message0("Input file does not exist!");
		}
	fscanf(fpin, "%d", &npost);
	for (i = 0; i < npost; i++)
		{
			fscanf(fpin, "%d %d %d %f %f %f %f %f %f", &a[i], &b[i], &c[i], &d[i], &e[i], &f[i], &g[i], &x[i], &y[i]);
		}
	for (i = 0; i< npost; i++)
		{
			Message0("%d\t %d\t %d\t %f\t %f\t %f\t %f\t %f\t %f\n ", a[i], b[i], c[i], d[i], e[i], f[i], g[i], x[i], y[i]);
		}
}

