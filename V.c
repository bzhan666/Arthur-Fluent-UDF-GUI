 #include "udf.h"
 DEFINE_PROFILE(vel_profile,t,i) 
 {  
 real x[ND_ND];  
 real y;
 face_t f;    
 begin_f_loop(f,t)     
  {           
  	F_CENTROID(x,f,t);
  	y=x[1];
  	F_PROFILE(f,t,i) = 20-y*y/(.03*.03)*20;      
  }    
  	end_f_loop(f,t) 
 }
 
