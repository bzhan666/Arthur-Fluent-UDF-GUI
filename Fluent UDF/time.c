 #include "udf.h"
 DEFINE_DELTAT(mydeltat,d)
 {
 	real time_step;
 	real flow_time = CURRENT_TIME;
 	if (flow_time < 0.5)
 		time_step = 0.1;
 	else
 		time_step = 0.2;
 	return time_step;
 } 
