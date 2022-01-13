#include "udf.h"

DEFINE_PROPERTY(saturation_temperature,c,thread)
{
	real temp = C_T(c,thread);
	real sat_temp,sa_temp;
	real w;

	Thread *t_air = THREAD_SUB_THREAD(thread,0); /* 给定相域索引，THREAD_SUB_THREAD宏可用于检索相级线程（子线程）指针。 THREAD_SUB_THREAD具有两个参数：mixture_thread和phase_domain_index。 该函数返回给定phase_domain_index的阶段级线程指针。*/
	Thread *t_wv = THREAD_SUB_THREAD(thread,1);
	Thread *t_wl = THREAD_SUB_THREAD(thread,2);

	real p_wv = C_P(c,t_wv); /* 定义一个水蒸气相的压力并获取 */

	real operating_p = RP_Get_Real ("operating-pressure");
	real vapor_p = p_wv + operating_p; // absolute pressure = static pressure + operating pressure
	//printf("vapor_p %f",p_wv);

	if(p_wv>0)
	{
		sa_temp = (1723.6425/(8.05573-((double)log10((double)p_wv/133.322)))-233.08)+273.15; // Pressure in Pa to mmHg (Yaws)

		sat_temp=sa_temp;  /* 把sa_temp赋值给sat_temp */


	}
	else
		sat_temp=273.15;
	//printf("sat_temp %f",sat_temp);

	return sat_temp;

}