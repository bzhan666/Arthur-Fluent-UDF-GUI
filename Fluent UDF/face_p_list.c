#include "udf.h"
 #define WALLID 3
   DEFINE_ON_DEMAND(face_p_list)
   {
   #if !RP_HOST /* Host will do nothing in this udf. */
		face_t f;
		Thread *tf;
		Domain *domain;
		real *p_array;
		real x[ND_ND], (*x_array)[ND_ND];
		int n_faces, i, j;
		domain=Get_Domain(1); /* Each Node will be able to access
			   its part of the domain */
		tf=Lookup_Thread(domain, WALLID); /* Get the thread from the domain */
		/* The number of faces of the thread on nodes 1,2... needs to be sent
		to compute node-0 so it knows the size of the arrays to receive
		from each */
		n_faces=THREAD_N_ELEMENTS_INT(tf);
		/* No need to check for Principal Faces as this UDF
		will be used for boundary zones only */
	  if(! I_AM_NODE_ZERO_P) /* Nodes 1,2... send the number of faces */
		{
		 PRF_CSEND_INT(node_zero, &n_faces, 1, myid);
		}
	 /* Allocating memory for arrays on each node */  
	 p_array=(real *)malloc(n_faces*sizeof(real));
	  x_array=(real (*)[ND_ND])malloc(ND_ND*n_faces*sizeof(real));
	  begin_f_loop(f, tf)
		 /* Loop over interior faces in the thread, filling p_array
		 with face pressure and x_array with centroid  */
		{
		   p_array[f] = F_P(f, tf);
		   F_CENTROID(x_array[f], f, tf);
		}
	  end_f_loop(f, tf)
   /* Send data from node 1,2, ... to node 0 */
   Message0("\nstart\n");
	  if(! I_AM_NODE_ZERO_P) /* Only SEND data from nodes 1,2... */
	  {
		 PRF_CSEND_REAL(node_zero, p_array, n_faces, myid);
		 PRF_CSEND_REAL(node_zero, x_array[0], ND_ND*n_faces, myid);
	  }
	else
	  {/* Node-0 has its own data,
		so list it out first */
	   Message0("\n\nList of Pressures...\n");
	   for(j=0; j<n_faces; j++)
		  /* n_faces is currently node-0 value */
		 {
 # if RP_3D
		   Message0("%12.4e %12.4e %12.4e %12.4e\n",
		   x_array[j][0], x_array[j][1], x_array[j][2], p_array[j]);
 # else /* 2D */
		   Message0("%12.4e %12.4e %12.4e\n",
		   x_array[j][0], x_array[j][1], p_array[j]);
 # endif
	   }
	}
   /* Node-0 must now RECV data from the other nodes and list that too */
	  if(I_AM_NODE_ZERO_P)
	  {
		compute_node_loop_not_zero(i)
		 /* See para.h for definition of this loop */
		{
		   PRF_CRECV_INT(i, &n_faces, 1, i);
			/* n_faces now value for node-i */
		   /* Reallocate memory for arrays for node-i */
				 p_array=(real *)realloc(p_array, n_faces*sizeof(real));
		   x_array=(real(*)[ND_ND])realloc(x_array,ND_ND*n_faces*sizeof(real));
		   /* Receive data */
		   PRF_CRECV_REAL(i, p_array, n_faces, i);
		   PRF_CRECV_REAL(i, x_array[0], ND_ND*n_faces, i);
		   for(j=0; j<n_faces; j++)
			 {
 # if RP_3D
				Message0("%12.4e %12.4e %12.4e %12.4e\n",
				x_array[j][0], x_array[j][1], x_array[j][2], p_array[j]);
 # else /* 2D */
				Message0("%12.4e %12.4e %12.4e\n",
				x_array[j][0], x_array[j][1], p_array[j]);
 # endif
			 }
		 }
	  }
	  free(p_array); /* Each array has to be freed before function exit */
	  free(x_array);
 #endif /* ! RP_HOST */
   }

