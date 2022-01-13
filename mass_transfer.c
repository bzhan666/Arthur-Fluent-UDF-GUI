#include "udf.h"
 DEFINE_MASS_TRANSFER(liq_gas_source, cell, thread, from_index,from_species_index, to_index, to_species_index)
 {
 	real m_lg;
 	real T_SAT = 373.15;
 	Thread *gas, *liq;
 	gas = THREAD_SUB_THREAD(thread, from_index);
 	liq = THREAD_SUB_THREAD(thread, to_index);
 	m_lg = 0.0;
 	if (C_T(cell, liq) > T_SAT)
 		{ /* Evaporating */
 			m_lg = -0.1*C_VOF(cell,liq)*C_R(cell,liq)*
 					(C_T(cell,liq)-T_SAT)/T_SAT;
 		}
 	else if (C_T(cell, gas) < T_SAT)
 		{ /* Condensing */
 		m_lg = 0.1*C_VOF(cell,gas)*C_R(cell,gas)*
					(T_SAT-C_T(cell,gas))/T_SAT;
 		}
 return (m_lg);
 }

