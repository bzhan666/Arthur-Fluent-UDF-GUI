#include "udf.h"
# define FLUID_ID 2

DEFINE_ON_DEMAND(pressures_to_file)
{
#if !RP_HOST
	Domain *domain=Get_Domain(1);
	Thread *thread;
	cell_t c;
#else
	int i;
#endif

#if !RP_NODE
	FILE *fp = NULL;
	char filename[]="press_out.txt";
#endif

	int size; /* data passing variables */
	real *array;
	int pe;

#if !RP_HOST
	thread=Lookup_Thread(domain,FLUID_ID);
#endif 
#if !RP_NODE 
	if ((fp = fopen(filename, "w"))==NULL)
		Message("\n Warning: Unable to open %s for writing\n",filename);
	else
		Message("\nWriting Pressure to %s...",filename);
#endif

#if RP_NODE
	size=THREAD_N_ELEMENTS_INT(thread);
	array = (real *)malloc(size * sizeof(real));
	begin_c_loop_int(c,thread)
		array[c]= C_P(c,thread);
	end_c_loop_int(c,thread)

	pe = (I_AM_NODE_ZERO_P) ? host : node_zero;
	PRF_CSEND_INT(pe, &size, 1, myid);
	PRF_CSEND_REAL(pe, array, size, myid);
	free(array);

if (I_AM_NODE_ZERO_P) 
	compute_node_loop_not_zero (pe)
	{
		PRF_CRECV_INT(pe, &size, 1, pe);
		array = (real *)malloc(size * sizeof(real));
		PRF_CRECV_REAL(pe, array, size, pe);
		PRF_CSEND_INT(host, &size, 1, myid);
		PRF_CSEND_REAL(host, array, size, myid);
		free((char *)array);
	}
#endif 

#if RP_HOST
	compute_node_loop (pe)
	{
		PRF_CRECV_INT(node_zero, &size, 1, node_zero);
		array = (real *)malloc(size * sizeof(real));
		PRF_CRECV_REAL(node_zero, array, size, node_zero);

	for (i=0; i<size; i++)
		fprintf(fp, "%g\n", array[i]);
	free(array);
	}
#endif 

#if !RP_NODE
	fclose(fp);
	Message("Done\n");
#endif
}
