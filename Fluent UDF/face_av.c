#include "udf.h"

DEFINE_ADJUST(face_av, domain)
{
int surface_thread_id=0;
real total_area=0.0;
real total_force=0.0;  /*变量由主机和节点进程使用。计算节点将变量用于计算目的，主机将变量用于消息传递和显示目的*/

#if !RP_HOST
	Thread* thread;
	face_t face;
	real area[ND_ND];  /*预处理器仅针对节点版本（而不是主机）编译线程，面部和区域变量，因为面部和线程仅在FLUENT的节点版本中定义*/
#endif 

#if !RP_NODE
	surface_thread_id = RP_Get_Integer("surface_thread_id");
	Message("\nCalculating on Thread # %d\n",surface_thread_id);/*使用RP_GET_INTEGER实用程序的主机进程获得名为PRES_AV/thread_ID的用户定义方案变量,并分配给变量surface_thread_id*/
#endif 

host_to_node_int_1(surface_thread_id);

#if RP_NODE 
Message("\nNode %d is calculating on thread # %d\n",myid,surface_thread_id);
#endif

#if !RP_HOST
thread = Lookup_Thread(domain,surface_thread_id);
begin_f_loop(face,thread)
	if (PRINCIPAL_FACE_P(face,thread)) /*用于确保分区边界处的面部不计算两次*/
		{
		F_AREA(area,face,thread);
		total_area += NV_MAG(area);
		total_force += NV_MAG(area)*F_P(face,thread);
		}
	
end_f_loop(face, thread) /*使用#if！rp_host，节点进程循环到与surface_thread_id相关联的线程中的所有面部循环，并计算总区和总力*/
Message("Total Area Before Summing %f\n",total_area);
Message("Total Normal Force Before Summing %f\n",total_force);
	
# if RP_NODE 
	total_area = PRF_GRSUM1(total_area);   /*全局求和宏*/
    total_force = PRF_GRSUM1(total_force);
# endif

#endif 

node_to_host_real_2(total_area,total_force);

#if !RP_NODE
Message("Total Area After Summing: %f (m2)\n",total_area);
Message("Total Normal Force After Summing %f (N)\n",total_force);
Message("Average pressure on Surface %d is %f (Pa)\n",
     surface_thread_id,(total_force/total_area));
#endif 
}
