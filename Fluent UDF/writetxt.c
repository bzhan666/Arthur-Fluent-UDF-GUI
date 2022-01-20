#include "udf.h"

DEFINE_PROFILE(temp_wall_profile, t, i)
{
	real x[ND_ND];
	face_t f;
	begin_f_loop(f, t)
	{
		F_CENTROID(x, f, t);
		F_PROFILE(f, t, i) = -3200*x[0]*x[0]+1600*x[0]+300;
	}
	end_f_loop(f, t)
}

DEFINE_SOURCE(heat_source, c, t, dS, eqn)
{
	const real C1 = 1.0;
	const real C2 = 0.03;
	real source;
	real cst_aux = C1 / (C2-C_U(c, t));
	source = cst_aux * C_T(c, t);
	dS[eqn] = cst_aux;
	return source;
}

DEFINE_ON_DEMAND(write_temperature)

{
	FILE *fp;
	Thread *t;
	Domain *d;
	cell_t c;

	const char FILENAME[] = "temperature.txt";
	real position[ND_ND];

	
	if(N_UDM < 1)
	{
		Message0("\n Error: No UDM defined! Abort UDF execution.\n");
		return;
	}

	if((fp = fopen(FILENAME, "w")) == NULL)
	{
		Message("\n Error: No write access to file %s. Abort UDF execution.\n", FILENAME);
		return;
	}
	else
	{
		fprintf(fp, "x m\ty m\tTemperature degC\n");
	}

	d = Get_Domain(1);
	thread_loop_c(t, d)
	{
		begin_c_loop(c, t)
		{
			C_CENTROID(position, c, t);
			fprintf(fp, "%10.5e\t%10.5e\t%10.5e\n", position[0], position[1], C_T(c, t)-273.15);
		}
		end_c_loop(c, t)
	}

	fclose(fp);
	Message0(" \n Finished writing to file temperature.txt\n");
}	

