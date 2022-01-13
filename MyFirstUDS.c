#include "udf.h"
DEFINE_UDS_UNSTEADY(MyUnsteady,c,t,i,apu,su)
{
real physical_dt, vol, phi_old;
physical_dt = RP_Get_Real("physical-time-step");
vol = C_VOLUME(c,t);
*apu = -vol / physical_dt; /*implicit part*/
phi_old = C_STORAGE_R(c,t,SV_UDSI_M1(i));
*su = vol*phi_old/physical_dt; /*explicit part*/
}

DEFINE_UDS_FLUX(MyFlux,f,t,i)
{
real NV_VEC(unit_vec), NV_VEC(A); //声明矢量变量
F_AREA(A, f, t);
NV_DS(unit_vec, =, 1, 1, 1, *, 1); //单位矢量赋值
return NV_DOT(unit_vec, A); //矢量点积
}

DEFINE_DIFFUSIVITY(MyDiff,c,t,i)
{
return 1.0;
}

DEFINE_SOURCE(MySource,c,t,dS,eqn)
{
 dS[eqn]=0;
 return 2.0;
}

DEFINE_PROFILE(MyProfile,thread,index)
{
real x[ND_ND]; /* this will hold the position vector */
real xx,yy;
face_t f;
begin_f_loop(f,thread) 
{
F_CENTROID(x,f,thread);
xx = x[0];
yy=x[1];
F_PROFILE(f,thread,index) = (xx+yy); 
}
end_f_loop(f,thread)
}

