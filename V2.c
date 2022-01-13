#include "udf.h" 
 
DEFINE_PROFILE(unsteady_velocity, thread, position) 
{   
  face_t f; 
  begin_f_loop(f, thread)     
  {
       real t=CURRENT_TIME;      
       F_PROFILE(f, thread, position) = 20. + 5.0*sin(10.*t);
  }
   end_f_loop(f, thread) 
}