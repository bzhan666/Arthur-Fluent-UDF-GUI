#include "udf.h"

#define ZoneID 10
#define TargetStaticPreesure 500
#define InitialTotalPressure 600

real staPA = 0;
real totPA = InitialTotalPressure;
real Target = TargetStaticPreesure;
DEFINE_EXECUTE_AT_END(inavgPa)

{
	real NV_VEC(A);
	real tot_p = 0.0, tot_a = 0.0;

#if !RP_HOST

	face_t f;
	Thread *t;
	Domain *d;
	d = Get_Domain(1);
	t = Lookup_Thread(d, ZoneID);
	{
	begin_f_loop(f, t)

	if(PRINCIPAL_FACE_P(f, t))
	{
		
		F_AREA(A, f, t);
		tot_a += NV_MAG(A);
		tot_p += NV_MAG(A)*F_P(f, t);
		
	}
	
	end_f_loop(f, t)
	
	}
	
tot_a = PRF_GRSUM1(tot_a);
tot_p = PRF_GRSUM1(tot_p);

#endif
	
node_to_host_real_2(tot_p, tot_a);

	staPA = tot_p / tot_a; /*进口壁面平均压力 */
	totPA -=0.5*(staPA-Target);

#if RP_HOST
	Message("\n static is %f Pa.\n",staPA);
	Message("\n total is %f Pa.\n",totPA);
#endif
	node_to_host_real_1(totPA);
	host_to_node_real_1(staPA);
}

DEFINE_PROFILE(tp, t, i)
{
	face_t f;
	
	begin_f_loop(f, t)
	{
		F_PROFILE(f, t, i) = totPA;
	}
	end_f_loop(f, t)
}
