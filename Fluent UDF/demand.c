#include "udf.h"

DEFINE_ON_DEMAND(shuge)
{
#if RP_HOST
	Message("hello a shuge!\n");
#endif
}
