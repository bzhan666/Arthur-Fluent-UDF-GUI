#include "udf.h"
DEFINE_SOURCE(cell_x_source,c,t,dS,eqn)
{
    real con=20;
    real source;
    if(C_T(c,t)<=288)
    {
        source=-con*C_U(c,t);
        dS[eqn]=-con;
    }
    else
    {
        source=dS[eqn]=0;
    }
    return source;
}

DEFINE_PROPERTY(cell_viscosity,c,t)
{
    real mu_lam;
    real temp=C_T(c,t);
    if (temp>288.)
    mu_lam=5.5e-3;
    else if (temp>286.)
    mu_lam=143.2135-0.49725*temp;
    else
    mu_lam=1.;
    return mu_lam;
}
