#include "udf.h"

DEFINE_DIFFUSIVITY(conduction, c, t, i)
{
	return 2;
}

DEFINE_UDS_UNSTEADY(temp_unsteady, c, t, i, apu, su)
{
	real dt,pre_temp,vol;
	dt = CURRENT_TIMESTEP;  //获取时间不长
	vol = C_VOLUME(c, t);   //获取单元格体积
	*apu = -C_R(c, t) * vol / dt;  //隐式部分
	pre_temp = C_UDSI_M1(c, t, i);
	*su = C_R(c,t) * vol * pre_temp / dt;  //显示部分
}


