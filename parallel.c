/*****************************************************
  Simple UDF that uses compiler directives
 *****************************************************/
 #include "udf.h"
 DEFINE_ADJUST(where_am_i, domain)
 {
 #if RP_HOST
    Message("I am in the host process\n");
 #endif /* RP_HOST */
 #if RP_NODE
    Message("I am in the node process with ID %d\n",myid);
    /* myid is a global variable which is set to the multiport ID for
    each node */
 #endif /* RP_NODE */
 }

