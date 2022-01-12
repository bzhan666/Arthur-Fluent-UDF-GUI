#include "udf.h"

DEFINE_ADJUST(face_av, domain)
{
int surface_thread_id=0.0;
real total_area=0.0;
real total_force=0.0;

#if !RP_HOST
	Thread* thread;
	face_t face;
	real area[ND_ND];
#endif

#if !RP_NODE
	surface_thread_id = RP_Get_Integer("pres_av/thread-id");
	Message("\nCalculating on Thread # %d\n",surface_thread_id);
#endif

host_to_node_int_1(surface_thread_id);

#if RP_NODE 
	Message("\nNode %d is calculating on thread # %d\n",myid,surface_thread_id);
#endif

#if !RP_HOST
thread = Lookup_Thread(domain,surface_thread_id);
begin_f_loop(face, thread)
	
	if (PRINCIPAL_FACE_P(face, thread))
		{
		F_AREA(area, face, thread);
		total_area += NV_MAG(area);
		total_force += NV_MAG(area)*F_P(face,thread);
		}
	
end_f_loop(face, thread)
Message("Total Area Before Summing %f\n",total_area);
Message("Total Normal Force Before Summing %\n",total_force);

#if RP_NODE
	total_area = PRF_GRSUM1(total_area);
	total_force = PRF_GRSUM1(total_force);
#endif

#endif

#if !RP_NODE
Message("Total Area After Summing: %f (m2)\n",total_area);
Message("Total Normal Force After Summing %f (N)\n",total_force);
Message("Average pressure on Surface %d is %f (Pa)\n",
	surface_thread_id,(total_force/total_area));
#endif
}












