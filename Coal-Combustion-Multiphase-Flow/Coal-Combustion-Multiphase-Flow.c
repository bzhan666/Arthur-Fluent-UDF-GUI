/*
Fluent-UDF_煤炭燃烧多相流过程
---各种非均相反应的代码
按照不同的物理化学过程（如水分释放、脱挥发分、焦炭燃烧、气化等）组织成
多个独立的 DEFINE_HET_RXN_RATE 宏，每个函数负责一个特定的反应。
考虑了煤燃烧过程中的多个复杂现象，包括水分蒸发、脱挥发分、焦炭与多种气体
（O2, H2O, CO2, H2）的非均相反应、以及挥发物（焦油、CH4等）的均相燃烧。
*/
#include "udf.h"
#include "stdio.h"
#include "time.h"

#define SMALL_S 1.e-29 // 一个很小的数，用于避免除以零
#define eps_g_small 1. // 气体体积分数下限
#define eps_s_small 1.e-6 // 固体体积分数下限
#define spe_small 1.e-6 // 组分质量分数下限
#define spe_small_comb 1.e-6 // 燃烧过程中组分质量分数下限
#define TMAX 3000. // 最高温度
#define TMIN  280. // 最低温度
#define Rgas 1.987 /* cal/mol.K */ /* 通用气体常数 */

#define PR_NUMBER(cp,mu,k) ((cp)*(mu)/(k)) // 普朗特数计算
#define IP_HEAT_COEFF(vofP,vofS,k,nu,d) ((vofP)*(vofS)*6.*(k)*(Nu)/(d)/(d)) // 相间传热系数

#define Carbon_Init_MF 0.541 // 初始碳质量分数
#define Volatile_Init_MF 0.418 // 初始挥发分质量分数
#define Moisture_Init_MF 0.026 // 初始水分质量分数
#define Ash_Init_MF (1.-Carbon_Init_MF-Volatile_Init_MF-Moisture_Init_MF) // 初始灰分质量分数

#define rho_c  (1000.) // 碳密度
#define rho_ash (2931) // 灰分密度
#define rho_liq_water (998.2) // 液态水密度
#define rho_volatile (1000.) // 挥发分密度
#define solid_rho (1./(Carbon_Init_MF/rho_c+Volatile_Init_MF/rho_volatile+Moisture_Init_MF/rho_liq_water+Ash_Init_MF/rho_ash)) // 固体密度

/* 函数定义 */

// 固体燃料反应物
void SolidFuel_Reactant(cell_t c, Thread *t, Hetero_Reaction *hr, double* y_carbon, double* mol_weight);
// 挥发分质量分数
void volatile_mass_fractions();
// 设置组分索引
void SetSpeciesIndex();
// 读取c3m数据
void read_c3m_data();
// 饱和压力
double satPressure(double T);
// 获取相索引
double Get_Phase_Index(Hetero_Reaction *hr);
// 湍流反应速率
double Turbulent_rr(cell_t c, Thread *t, Hetero_Reaction *r, real yi[MAX_PHASES][MAX_SPE_EQNS]);
// 质量传递系数
double Mass_Transfer_Coeff(cell_t c, Thread *tp, Thread *ts);
// Gunn传热模型UDF
double heat_gunn_udf(cell_t c, Thread *ti, Thread *tj);
// 燃烧反应速率
double  rr_combustion(cell_t c, Thread *t, Thread *ts, Thread *tp, double yi_O2, double y_ash, double y_carbon);
// 水蒸气气化反应速率
double  rr_steam_gasif(cell_t c, Thread *t, Thread *ts, Thread *tp, double p_h2o, double p_co, double p_h2, double y_carbon, double mol_weight, double* direction);
// CO2气化反应速率
double  rr_co2_gasif(cell_t c, Thread *t, Thread *ts, Thread *tp, double y_co, double y_co2, double y_carbon, double mol_weight, double* direction);
// H2气化反应速率
double  rr_h2_gasif(cell_t c, Thread *t, Thread *ts, Thread *tp, double y_h2, double y_ch4, double y_carbon, double mol_weight, double* direction);

// 组分索引变量 (IP: 相索引, IS: 组分索引)
int IP_CH4 = 0, IS_CH4 = 0, IP_CO = 0, IS_CO = 0, 
    IP_CO2 = 0, IS_CO2 = 0, IP_H2 = 0, IS_H2 = 0,
    IP_H2O = 0, IS_H2O = 0, IP_O2 = 0, IS_O2 = 0, 
    IP_H2S = 0, IS_H2S = 0, IP_CL2 = 0, IS_CL2 = 0, 
    IP_NH3 = 0, IS_NH3 = 0, IP_N2 = 0, IS_N2 = 0, 
    IP_TAR = 0, IS_TAR = 0, IP_C = 0, IS_C = 0, 
    IP_VOL = 0, IS_VOL = 0, IP_MOISTURE = 0, IS_MOISTURE = 0, 
    IP_ASH = 0, IS_ASH = 0, IP_SOOT = 0, IS_SOOT = 0,
    IP_PAH = 0, IS_PAH = 0, IP_OIL = 0, IS_OIL = 0,
    IP_ASH_R = 0, IS_ASH_R = 0, IP_C_R = 0, IS_C_R = 0,
    IP_C2H2 = 0, IS_C2H2 = 0, IP_SLURRY = 0, IS_SLURRY = 0,    
    IP_C2H4 = 0, IS_C2H4 = 0, IP_C2H6 = 0, IS_C2H6 = 0,
    IP_C3H6 = 0, IS_C3H6 = 0, IP_C3H8 = 0, IS_C3H8 = 0,  
    IP_SLURRY_D = 0, IS_SLURRY_D = 0, IP_SAND = 0, IS_SAND = 0, 
    IP_ASH_D = 0, IS_ASH_D = 0;

int current_c3m_solid_phase = 0; // 当前c3m固体相

cxboolean init_flag                = TRUE; // 初始化标志
cxboolean PCCL_Devol               = FALSE; // PCCL脱挥发分模型开关
cxboolean MGAS_Devol               = FALSE; // MGAS脱挥发分模型开关
cxboolean CPD_Devol                = FALSE; // CPD脱挥发分模型开关
cxboolean FGDVC_Devol              = FALSE; // FGDVC脱挥发分模型开关
cxboolean HPTR_Devol               = FALSE; // HPTR脱挥发分模型开关
cxboolean MGAS_Moisture            = FALSE; // MGAS水分释放模型开关
cxboolean PCCL_Moisture            = FALSE; // PCCL水分释放模型开关
cxboolean MGAS_TarCracking         = FALSE; // MGAS焦油裂解模型开关
cxboolean PCCL_2nd_Pyro            = FALSE; // PCCL二次热解模型开关
cxboolean MGAS_Gasif               = FALSE; // MGAS气化模型开关
cxboolean PCCL_Gasif               = FALSE; // PCCL气化模型开关
cxboolean PCCL_TarCracking         = FALSE; // PCCL焦油裂解模型开关
cxboolean MGAS_WGS                 = FALSE; // MGAS水煤气变换反应模型开关
cxboolean PCCL_soot_gasif          = FALSE; // PCCL碳烟气化模型开关
cxboolean MGAS_char_combustion     = FALSE; // MGAS焦炭燃烧模型开关
cxboolean PCCL_char_combustion     = FALSE; // PCCL焦炭燃烧模型开关
cxboolean PCCL_soot_oxidation      = FALSE; // PCCL碳烟氧化模型开关
cxboolean TAR_oxidation            = FALSE; // 焦油氧化模型开关
cxboolean MGAS_gas_phase_oxidation = FALSE; // MGAS气相氧化模型开关
 
double avg_mf_h2o,avg_mf_co,avg_mf_h2,avg_mf_ch4,avg_mf_co2,avg_c,avg_volatile,avg_moisture,avg_ash; // 平均质量分数

/* 煤分析变量 */
double fc_ar=0.,vm_ar=0.,ash_ar=0.,moist_ar=0.; // 固定碳、挥发分、灰分、水分 (收到基)
double f_ep_a = 0.;
double mw[MAX_PHASES][MAX_SPE_EQNS]; // 分子量

double A1_devolatilization=0.0, E1_devolatilization=0.0;  /* pan : Oct 2012 */ // 脱挥发分反应动力学参数
double A2_devolatilization=0.0, E2_devolatilization=0.0;  /* pan : Oct 2012 */ // 脱挥发分反应动力学参数

double A_tar_cracking=0.0 , E_tar_cracking=0.0; // 焦油裂解反应动力学参数

double A_steam_gasification=0.0, E_steam_gasification=0.0; // 水蒸气气化反应动力学参数
double K_steam_gasification=0.0, N_steam_gasification=0.0;  /* pan : Oct 2012 */
double Annealing_steam_gasification=1.0;                /* pan : Oct 2012 */

double A_co2_gasification=0.0, E_co2_gasification=0.0; // CO2气化反应动力学参数
double K_co2_gasification=0.0, N_co2_gasification=0.0;    /* pan : Oct 2012 */
double Annealing_co2_gasification=1.0;                /* pan : Oct 2012 */

double A_h2_gasification=0.0, E_h2_gasification=0.0; // H2气化反应动力学参数
double N_h2_gasification=0.0;                        /* pan : Oct 2012 */
double Annealing_h2_gasification=1.0;                /* pan : Oct 2012 */

double A_soot_steam_gasification=0.0, E_soot_steam_gasification=0.0;  /* pan : oct 2012 */ // 碳烟水蒸气气化反应动力学参数
double K_soot_steam_gasification=0.0, N_soot_steam_gasification=0.0;  /* pan : Oct 2012 */
double Annealing_soot_steam_gasification=0.0;                /* pan : Oct 2012 */

double A_soot_co2_gasification=0.0, E_soot_co2_gasification=0.0;    /* pan : Oct 2012 */ // 碳烟CO2气化反应动力学参数
double K_soot_co2_gasification=0.0, N_soot_co2_gasification=0.0;    /* pan : Oct 2012 */
double Annealing_soot_co2_gasification=0.0;                /* pan : Oct 2012 */

double A_soot_h2_gasification=0.0, E_soot_h2_gasification=0.0;     /* pan : Oct 2012 */ // 碳烟H2气化反应动力学参数
double N_soot_h2_gasification=0.0;                        /* pan : Oct 2012 */
double Annealing_soot_h2_gasification=0.0;                /* pan : Oct 2012 */

double A_Soot_Combustion = 0.0, E_Soot_Combustion = 0.0; // 碳烟燃烧反应动力学参数

double A_c_combustion = 8710. /* g/(atm.cm^2.s) */, E_c_combustion = 27000. /* cal/mole */; // 碳燃烧反应动力学参数
double Annealing_c_combustion=0.0 , N_c_combustion=0.0;    /* pan : Oct 2012 */

double A_moisture_release = 0.0, E_moisture_release = 0.0; // 水分释放动力学参数
double wg3 = 0.014;
double Moisture_Flux; // 水分通量
   
   
DEFINE_ADJUST(gasification,domain)
{ 
 if(init_flag)
  {
#if !RP_HOST
    if(0) SetSpeciesIndex(); 
#endif

#if !RP_NODE       
    volatile_mass_fractions();
#endif
    // 将主机数据传输到计算节点
    host_to_node_real_5(fc_ar,vm_ar,ash_ar,moist_ar,f_ep_a);
    host_to_node_int_4(PCCL_Devol,MGAS_Devol,MGAS_Moisture,PCCL_2nd_Pyro);
    host_to_node_int_3(MGAS_Gasif,PCCL_Gasif,PCCL_Moisture); 
    host_to_node_int_3(CPD_Devol,FGDVC_Devol,HPTR_Devol);    /* pan : Oct 2012 */
    host_to_node_int_3(PCCL_TarCracking,MGAS_WGS,PCCL_soot_gasif);    /* pan : Oct 2012 */
    host_to_node_int_1(PCCL_char_combustion);    /* pan : Oct 2012 */
    host_to_node_int_3(MGAS_char_combustion,TAR_oxidation,MGAS_gas_phase_oxidation);    /* pan : Oct 2012 */
    host_to_node_real_2(A_moisture_release, E_moisture_release);
    host_to_node_real_2(A1_devolatilization, E1_devolatilization);  /* pan : Oct 2012 ... added the "1" */
    host_to_node_real_2(A2_devolatilization, E2_devolatilization);  /* pan : Oct 2012 */
    host_to_node_real_2(A_tar_cracking, E_tar_cracking);
    host_to_node_real_2(A_steam_gasification, E_steam_gasification);
    host_to_node_real_2(K_steam_gasification, N_steam_gasification);               /* pan : Oct 2012 */
    host_to_node_real_1(Annealing_steam_gasification);                               /* pan : Oct 2012 */
    host_to_node_real_2(A_co2_gasification, E_co2_gasification);
    host_to_node_real_2(K_co2_gasification, N_co2_gasification);                   /* pan : Oct 2012 */
    host_to_node_real_1(Annealing_co2_gasification);                                 /* pan : Oct 2012 */
    host_to_node_real_2(A_h2_gasification, E_h2_gasification);
    host_to_node_real_2(Annealing_h2_gasification, N_h2_gasification);             /* pan : Oct 2012 */
    host_to_node_real_2(A_soot_steam_gasification, E_soot_steam_gasification);     /* pan : Oct 2012 */
    host_to_node_real_2(K_soot_steam_gasification, N_soot_steam_gasification);     /* pan : Oct 2012 */
    host_to_node_real_1(Annealing_soot_steam_gasification);                          /* pan : Oct 2012 */
    host_to_node_real_2(A_soot_co2_gasification, E_soot_co2_gasification);         /* pan : Oct 2012 */
    host_to_node_real_2(K_soot_co2_gasification, N_soot_co2_gasification);         /* pan : Oct 2012 */
    host_to_node_real_1(Annealing_soot_co2_gasification);                            /* pan : Oct 2012 */
    host_to_node_real_2(A_soot_h2_gasification, E_soot_h2_gasification);           /* pan : Oct 2012 */  
    host_to_node_real_2(Annealing_soot_h2_gasification, N_soot_h2_gasification);   /* pan : Oct 2012 */
    host_to_node_real_2(A_Soot_Combustion, E_Soot_Combustion);
    host_to_node_real_2(A_c_combustion, E_c_combustion);                           /* pan : Oct 2012 */
    host_to_node_real_2(N_c_combustion, Annealing_c_combustion);                   /* pan : Oct 2012 */
    host_to_node_real_2(wg3,Moisture_Flux);
    init_flag = FALSE; // 初始化完成
  }  /* end init_flag */     
}

DEFINE_ON_DEMAND(Devol_and_Tar_Cracking)
{
#if !RP_NODE
       
       volatile_mass_fractions();
       
#endif
}

// 设置组分索引
void SetSpeciesIndex()
{
  Domain *domain = Get_Domain(1);
   int n, ns, MAX_SPE_EQNS_PRIM = 0;
   Domain *subdomain;
   
   /*int n_phases = DOMAIN_N_DOMAINS(domain);*/

        /* 搜索所有组分并保存分子量 */
        sub_domain_loop(subdomain, domain, n)
           {
               Material *m_mat, *s_mat;
               if (DOMAIN_NSPE(subdomain) > 0)
                  {
                     m_mat = Pick_Material(DOMAIN_MATERIAL_NAME(subdomain),NULL);
                     mixture_species_loop(m_mat,s_mat,ns)
                        {
                            if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"ch4"))
                                {
                                   IP_CH4 = n; 
                                   IS_CH4 = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"co"))
                                {
                                   IP_CO = n; 
                                   IS_CO = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"c2h2"))
                                {
                                   IP_C2H2 = n; 
                                   IS_C2H2 = ns;
 				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                               }
                             else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"c2h4"))
                                {
                                   IP_C2H4 = n; 
                                   IS_C2H4 = ns;
 				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                               }
                             else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"c2h6"))
                                {
                                   IP_C2H6 = n; 
                                   IS_C2H6 = ns;
 				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                               }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"c3h6"))
                                {
                                   IP_C3H6 = n; 
                                   IS_C3H6 = ns;
 				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                               }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"c3h8"))
                                {
                                   IP_C3H8 = n; 
                                   IS_C3H8 = ns;
 				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                               }
                             else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"co2"))
                                {
                                   IP_CO2 = n; 
                                   IS_CO2 = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"h2"))
                                {
                                   IP_H2 = n; 
                                   IS_H2 = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"h2o"))
                                {
                                   IP_H2O = n; 
                                   IS_H2O = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"o2"))
                                {
                                   IP_O2 = n; 
                                   IS_O2 = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"h2s"))
                                {
                                   IP_H2S = n; 
                                   IS_H2S = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"cl2"))
                                {
                                   IP_CL2 = n; 
                                   IS_CL2 = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"nh3"))
                                {
                                   IP_NH3 = n; 
                                   IS_NH3 = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"n2"))
                                {
                                   IP_N2 = n; 
                                   IS_N2 = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"oil"))
                                {
                                   IP_OIL = n; 
                                   IS_OIL = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"pah"))
                                {
                                   IP_PAH = n; 
                                   IS_PAH = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"tar"))
                                {
                                   IP_TAR = n; 
                                   IS_TAR = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"c"))
                                {
                                   IP_C = n; 
                                   IS_C = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"c_recycle"))
                                {
                                   IP_C_R = n; 
                                   IS_C_R = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"soot"))
                                {
                                   IP_SOOT = n; 
                                   IS_SOOT = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"volatile"))
                                {
                                   IP_VOL = n; 
                                   IS_VOL = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"h2o<l>"))
                                {
                                   IP_MOISTURE = n; 
                                   IS_MOISTURE = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"h2o<l>-slurry"))
                                {
                                   IP_SLURRY = n; 
                                   IS_SLURRY = ns;
                            mw[n][ns] = MATERIAL_PROP(s_mat,PROP_mwi);			    
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"h2o<l>-dummy"))
                                {
                                   IP_SLURRY_D = n; 
                                   IS_SLURRY_D = ns;
                            mw[n][ns] = MATERIAL_PROP(s_mat,PROP_mwi);			    
                                }
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"ash-coal"))
                                {
                                   IP_ASH = n; 
                                   IS_ASH = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
				}
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"ash-dummy"))
                                {
                                   IP_ASH_D = n; 
                                   IS_ASH_D = ns;
				}                   
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"ash-recycle"))
                                {
                                   IP_ASH_R = n; 
                                   IS_ASH_R = ns;
				   if(n == 0)
				     MAX_SPE_EQNS_PRIM +=1;
                                }                             
                              else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"si<s>"))
                                {
                                   IP_SAND = n; 
                                   IS_SAND = ns;


// 定义水分释放的非均相反应速率
DEFINE_HET_RXN_RATE(moisture_release,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t); // 获取子线程
     int index_phase = Get_Phase_Index(hr); // 获取相索引
     Thread *ts = pt[index_phase]; /* 固体相 */
     double prod = 0.0, Ts = C_T(c,ts); // 产物，固体温度
     double Pt = MAX(0.1,(op_pres+C_P(c,t))); // 总压
     double Tsat = 1./(0.0727/log(Pt/611.) - 0.0042) + 273.; // 饱和温度
   
     *rr = 0;
  
 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
  */    
     if(IS_ASH == 0)
      SetSpeciesIndex(); 
	   
     if(C_VOF(c, ts) > 0.0 && Ts >= Tsat) // 如果固体体积分数大于0且固体温度大于等于饱和温度
      {      	   
       if(MGAS_Moisture) // 如果MGAS水分模型开启
        {
          prod  =  yi[IP_MOISTURE][IS_MOISTURE]*C_R(c,ts)/mw[IP_MOISTURE][IS_MOISTURE]; /* kg-mol/m^3*/
          A_moisture_release = A1_devolatilization;
          E_moisture_release = E1_devolatilization;
          *rr = A_moisture_release*exp(-E_moisture_release/(Rgas*Ts)) * prod*C_VOF(c, ts); /* kmol/(m3.s) */
        }
   
       if(PCCL_Moisture) // 如果PCCL水分模型开启
        {
          *rr = 6. * C_VOF(c, ts) / C_PHASE_DIAMETER(c,ts) * Moisture_Flux /mw[IP_MOISTURE][IS_MOISTURE] ;  /* kmol/(m3.s) */
        }
      } 
}


// 定义脱挥发分的非均相反应速率
DEFINE_HET_RXN_RATE(devolatilization,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);   
     int index_phase = Get_Phase_Index(hr);
     Thread *ts = pt[index_phase]; /* 固体相 */
     
     double x0_star = 0., x_star =0.;
     double Ts = C_T(c,ts);

     *rr = 0.0;
  
 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
  */
     
  if(IS_ASH == 0)
    SetSpeciesIndex(); 	   

  /* 
     挥发分 --> c1 焦油 + c2 CO + c3 CO2 + c4 CH4 + c5 H2 + c6 H2O + c7 H2S + c8 NH3
  */
  
     if(C_YI(c,ts,IS_MOISTURE) < spe_small) // 如果水分含量很小
      {
       if(MGAS_Devol) // 如果MGAS脱挥发分模型开启
        {
         if(Ts<1223)
          x0_star = pow((867.2/MAX((Ts-273.),1.)),3.914) / 100.;
          x_star = solid_rho * (fc_ar/100. + vm_ar/100.)/C_R(c,ts) * x0_star;
	  A2_devolatilization = 0.0;
	  E2_devolatilization = 0.0;
        }
	
       if((PCCL_Devol || CPD_Devol || FGDVC_Devol || HPTR_Devol) && yi[IP_VOL][IS_VOL] > x_star) // 如果其他脱挥发分模型开启且挥发分含量大于阈值
        {
         *rr = (A1_devolatilization *exp(-E1_devolatilization/(1.987*Ts)) + A2_devolatilization *exp(-E2_devolatilization/(1.987*Ts)))
               *(yi[IP_VOL][IS_VOL]-x_star)*C_VOF(c,ts)*C_R(c,ts)/mw[IP_VOL][IS_VOL];
        }
      }
}

// 定义焦油燃烧的非均相反应速率
DEFINE_HET_RXN_RATE(tar_comb,c,t,r,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0]; /* 气相 */

     double T_g = MAX(TMIN, C_T(c,tp)); // 气体温度
     double rho_g = C_R(c,tp)*1.e-3; /* g/cm^3 */ // 气体密度
     double rr_turb = 1e+20; // 湍流反应速率

     *rr = 0.0;
  
 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
  */
     
    if(IS_ASH == 0)
     SetSpeciesIndex(); 
	   
    if (rp_ke) // 如果使用k-epsilon模型
       rr_turb = Turbulent_rr(c, t, r, yi);	   
	 
    if(yi[IP_O2][IS_O2] > spe_small_comb) // 如果氧气含量足够
      {
       double tmp_exp, tmp_o2, tmp_tar;
       T_g = MIN(T_g, TMAX);
       tmp_exp = 3.8e11 * exp(-30000./(Rgas*T_g)); /* 30000 ---> 60000*/
       tmp_o2 = pow(rho_g*yi[IP_O2][IS_O2]/mw[IP_O2][IS_O2],1.5);
       tmp_tar = pow(rho_g*yi[IP_TAR][IS_TAR]/mw[IP_TAR][IS_TAR],0.25);

       *rr = tmp_exp * tmp_o2 * tmp_tar * C_VOF(c,tp);  /* mol/cm^3.s */
       *rr *= 1000.; /* kmol/(m^3 .s) */

       *rr = MIN(*rr, rr_turb); // 取动力学速率和湍流速率的较小值
      }
}

// 定义焦油裂解的非均相反应速率
DEFINE_HET_RXN_RATE(tar_cracking,c,t,r,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0]; /* 气相 */  
     double prod = 0.0;
     double T_g = MAX(TMIN, C_T(c,tp));
     double rr_turb = 1e+20;

     *rr = 0.0;
  
 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
  */
     
   if(IS_ASH == 0)
     SetSpeciesIndex();    
 
   if (rp_ke)
    rr_turb = Turbulent_rr(c, t, r, yi);	    
    
   if(yi[IP_TAR][IS_TAR] > spe_small) // 如果焦油含量足够
    {
      prod  =  yi[IP_TAR][IS_TAR]*C_R(c,tp)*C_VOF(c,tp)/mw[IP_TAR][IS_TAR];
      *rr = A_tar_cracking*exp(-E_tar_cracking/(1.987*T_g))
            * prod*C_VOF(c, tp); /* kmol/(m3.s) */
    }
}

// 定义CO燃烧的非均相反应速率
DEFINE_HET_RXN_RATE(co_comb,c,t,r,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0]; /* 气相 */
     
     double T_g = MAX(TMIN, C_T(c,tp));
     double rho_g = C_R(c,tp)*1.e-3; /* g/cm^3 */
     double rr_turb = 1e+20;
  
     *rr = 0;
  
 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
  */
     
     if(IS_ASH == 0)
      SetSpeciesIndex(); 
	   
     if(yi[IP_O2][IS_O2] > spe_small_comb) // 如果氧气含量足够
      {
       double tmp_exp, p_o2, p_co, p_h2o;
       if (rp_ke)
        rr_turb = Turbulent_rr(c, t, r, yi);	   
	  	  
       T_g = MIN(T_g, TMAX);
       tmp_exp = 3.98e+14 * exp(-40000./(Rgas*T_g)); /*40000 ---> 80000 */
       p_o2 = pow(rho_g*yi[IP_O2][IS_O2]/mw[IP_O2][IS_O2],0.25);
       p_co = rho_g*yi[IP_CO][IS_CO]/mw[IP_CO][IS_CO];
       p_h2o = pow(rho_g*yi[IP_H2O][IS_H2O]/mw[IP_H2O][IS_H2O], 0.5);
            
       *rr = tmp_exp * p_o2 * p_co * p_h2o * C_VOF(c,tp);  /* mol/cm^3.s */
       *rr *= 1000.; /* kmol/(m^3 .s) */
       *rr = MIN(*rr, rr_turb);	   
      }
}

// 定义CH4燃烧的非均相反应速率
DEFINE_HET_RXN_RATE(ch4_comb,c,t,r,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0]; /* 气相 */

     double T_g = MAX(TMIN, C_T(c,tp));
     double rho_g = C_R(c,tp)*1.e-3; /* g/cm^3 */
     double rr_turb = 1e+20;

     *rr = 0;
     
 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
  */
     
    if(IS_ASH == 0)
      SetSpeciesIndex(); 
	   
    if(yi[IP_O2][IS_O2] > spe_small_comb) // 如果氧气含量足够
     {
       double p_o2, p_ch4;
       if (rp_ke)
         rr_turb = Turbulent_rr(c, t, r, yi);	
	    
       T_g = MIN(T_g, TMAX);
       p_o2 = pow(rho_g*yi[IP_O2][IS_O2]/mw[IP_O2][IS_O2],1.3);
       p_ch4 = pow(rho_g*yi[IP_CH4][IS_CH4]/mw[IP_CH4][IS_CH4],0.2);

       *rr = 6.7e12 * exp(-48400./(Rgas*T_g)) * p_o2 * p_ch4 * C_VOF(c,tp);  /* mol/cm^3.s */
       *rr *= 1000.; /* kmol/(m^3 .s) */
       *rr = MIN(*rr, rr_turb);
     }
}    

// 定义H2燃烧的非均相反应速率
DEFINE_HET_RXN_RATE(h2_comb,c,t,r,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0]; /* 气相 */

     double T_g = MAX(TMIN, C_T(c,tp));
     double rho_g = C_R(c,tp)*1.e-3; /* g/cm^3 */
     double rr_turb = 1e+20;

     *rr = 0;
  
 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
  */     
     if(IS_ASH == 0)
       SetSpeciesIndex(); 
 
     if(yi[IP_O2][IS_O2] > spe_small_comb) // 如果氧气含量足够
      {
       double tmp_exp, tmp_o2, tmp_h2;
       if (rp_ke)
        rr_turb = Turbulent_rr(c, t, r, yi);   
     
       T_g = MIN(T_g, TMAX);
       tmp_exp = 1.08e16 * exp(-30000./(Rgas*T_g)); /* 30000 ----> 60000   cal/mole */    /* cm3/gmole-s */
       tmp_o2 = rho_g*yi[IP_O2][IS_O2]/mw[IP_O2][IS_O2];  /* gmole/cm3 */
       tmp_h2 = rho_g*yi[IP_H2][IS_H2]/mw[IP_H2][IS_H2];  /* gmole/cm3 */

       *rr = tmp_exp * tmp_o2 * tmp_h2 * C_VOF(c,tp);  /* mol/cm^3.s */
       *rr *= 1000.; /* kmol/(m^3 .s) */
       *rr = MIN(*rr, rr_turb);
      }
}

// 定义水煤气变换反应速率
DEFINE_HET_RXN_RATE(WGS_char,c,t,hr,mw,yi,rr,rr_t)
{
     Domain *domain = Get_Domain(1);
     Domain *subdomain;	
     Thread **pt = THREAD_SUB_THREADS(t);   
     Thread *tp = pt[0]; /* 气相 */
     double T_g = C_T(c,tp), f3 = 0., k3=1.;
     double Pt = MAX(0.1,(op_pres+C_P(c,t))/101325); // 总压 (atm)
     double p_co = 0., p_h2o = 0., p_co2 = 0., p_h2=0., sum = 0.0;
     int phase_domain_index, ns;
     double rr_turb = 1e+20;

     *rr = 0;
  
 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
  */     
     if(IS_ASH == 0)
       SetSpeciesIndex(); 
    
     sub_domain_loop(subdomain, domain, phase_domain_index)
      {
       if(phase_domain_index > 0)
        {
          Thread *ts = pt[phase_domain_index]; /* 固体相 */ 
          Material *m_mat, *s_mat;
          if (DOMAIN_NSPE(subdomain) > 0)
           {
             m_mat = Pick_Material(DOMAIN_MATERIAL_NAME(subdomain),NULL);
             mixture_species_loop(m_mat,s_mat,ns)
              {
               if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"ash-coal"))
                {
                  sum +=C_VOF(c,ts)*C_YI(c,ts,IS_ASH)*C_R(c,ts)*1.e-3;       
                }
	      }
           }
        }
      }
     
     if (rp_ke)
       rr_turb = Turbulent_rr(c, t, hr, yi);    
            
     if(C_VOF(c, tp) < 1.0)
      {
       T_g = MIN(T_g, TMAX);
       f3 = sum *exp(-8.91+5553/T_g); 
       k3 = exp( -3.63061 + 3955.71/T_g );

       p_co = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp)
              * yi[IP_CO][IS_CO]/mw[IP_CO][IS_CO] / 101325.;
       p_co2 = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp)
              * yi[IP_CO2][IS_CO2]/mw[IP_CO2][IS_CO2] / 101325.;
       p_h2o = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp)
              * yi[IP_H2O][IS_H2O]/mw[IP_H2O][IS_H2O] / 101325.;
       p_h2 = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp)
              * yi[IP_H2][IS_H2]/mw[IP_H2][IS_H2] / 101325.;
            
       *rr = 2.877e+5 * wg3 * f3 * pow (Pt, 0.5-Pt/250.) * exp(-27760/Rgas/T_g)*
             (p_co * p_h2o /Pt/Pt - p_co2 * p_h2 /Pt/Pt / k3) * C_VOF(c,tp);  /* mol/cm^3.s */
		  
       *rr *= 1000.; /* kmol/(m^3 .s) */
       *rr = MIN(*rr, rr_turb);
      }
}

// 定义碳烟燃烧的非均相反应速率
DEFINE_HET_RXN_RATE(soot_comb,c,t,r,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0]; /* 气相 */
     double T_g = MIN(TMAX, MAX(TMIN, C_T(c,tp)));
     double rr_turb = 1e+20;
     double soot_dia = 8.e-6; // 碳烟直径

     *rr = 0;
  
 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
  */     
     if(IS_ASH == 0)
      SetSpeciesIndex(); 
    
     double y_o2 = yi[IP_O2][IS_O2];   
     if(y_o2 > spe_small_comb) // 如果氧气含量足够
      {
       double tmp1_exp, tmp2_exp;
 
       if (rp_ke)
         rr_turb = Turbulent_rr(c, t, r, yi);  
     
       double P_o2 = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * T_g*y_o2/mw[IP_O2][IS_O2]/ 101325.; // 氧气分压
       tmp1_exp = 240.*exp(-15100./T_g)*P_o2/(1.+21.3*exp(2060./T_g)); 
       tmp2_exp = 5.35e-2*exp(-7640./T_g)*P_o2;
       double tmp3_exp = 1.51e5*exp(-48800./T_g);
       double tmp4_exp = 4.46e-3*exp(-7640./T_g)*P_o2;
       double chi = tmp4_exp/(tmp3_exp+tmp4_exp);
       double prod = tmp1_exp*chi+tmp2_exp*(1.-chi);   /* g/cm2-s */
       prod *= 10.;  /* kg/m2-s */
       *rr = (6.*C_VOF(c, tp)/soot_dia)*prod/mw[IP_SOOT][IS_SOOT];  /* kmol/m^3.s */
       *rr = MIN(*rr, rr_turb);
      }
}

// 定义水蒸气气化的非均相反应速率
DEFINE_HET_RXN_RATE(SteamGasif,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0]; /* 气相 */
     int index_phase = Get_Phase_Index(hr);
     Thread *ts = pt[index_phase]; /* 固体相 */

     *rr = 0;
     double direction = 0.0, mol_weight, y_carbon;

/*
                              C(s) + H2O ---> CO + H2 
 
设置相（phase）和物种（species）的索引。灰分（Ash）物种的索引初始值为 0，其余所有索引也初始为 0。
此处将灰分索引（IS_ASH）用作一个标志位：仅当其值仍为 0 时，才执行一次 SetSpeciesIndex() 函数，以完成所有物种索引的自动识别与赋值。
该初始化动作由用户在 Fluent 图形界面（GUI）中定义的第一个异相反应（heterogeneous reaction）触发执行，从而确保整个计算过程中仅初始化一次。
  */
     
     if(IS_ASH == 0)
      SetSpeciesIndex(); 
	   
     double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp);
     double p_h2o = RoRT * yi[IP_H2O][IS_H2O]/mw[IP_H2O][IS_H2O]/ 101325.; // 水蒸气分压
     double p_co = RoRT * yi[IP_CO][IS_CO]/mw[IP_CO][IS_CO] / 101325.; // CO分压
     double p_h2 = RoRT * yi[IP_H2][IS_H2]/mw[IP_H2][IS_H2] / 101325.;		       // H2分压
 

     SolidFuel_Reactant(c, t, hr, &y_carbon, &mol_weight); // 获取固体燃料反应物（碳）的质量分数和摩尔质量
 
     if(C_VOF(c, ts) >= eps_s_small) // 如果固体相的体积分数大于一个很小的值
      {
        *rr = rr_steam_gasif(c, t, ts, tp, p_h2o, p_co, p_h2, y_carbon, mol_weight, &direction); /* 计算水蒸气气化反应速率，单位 mol/(cm^3 .s) */

        if( direction < 0.0)   /* 负值表示逆向水蒸气气化 */
         *rr = 0.0; // 正向反应速率设为0
      }  
}

// 定义一个多相非均相反应速率的UDF，用于计算烟灰的水蒸气气化
DEFINE_HET_RXN_RATE(Soot_H2O_Gasif,c,t,r,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t); // 获取子线程
     Thread *tp = pt[0]; /* 气相 */

     *rr = 0; // 初始化反应速率
     double rr_turb = 1e+20; // 初始化湍流反应速率为一个大值
     double T_g = MIN((MAX(TMIN,C_T(c,tp))),TMAX); // 获取气体温度，并限制在最小和最大温度之间
 
/*
                              1/25 Soot + H2O ---> CO + H2 
 
    设置相和组分的索引。灰分(Ash)组分的索引被初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只被执行一次。这是由FLUENT GUI中
    非均相反应面板中定义的第一个反应完成的。
  */
     
     if(IS_ASH == 0) // 如果灰分组分索引未设置
      SetSpeciesIndex(); // 设置所有组分的索引
	   
     double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp); // 计算 Ro * R * T
     double p_h2o = RoRT * yi[IP_H2O][IS_H2O]/mw[IP_H2O][IS_H2O]/ 101325.; // 计算H2O的分压（atm）
     double p_h2 = RoRT * yi[IP_H2][IS_H2]/mw[IP_H2][IS_H2] / 101325.; // 计算H2的分压（atm）
 
     if (rp_ke) // 如果启用了k-epsilon湍流模型相关的反应模型
       rr_turb = Turbulent_rr(c, t, r, yi); // 计算湍流混合限制的反应速率
       
     double prod  =  yi[IP_SOOT][IS_SOOT]*(C_R(c,tp)*1e-03)/mw[IP_SOOT][IS_SOOT]*C_VOF(c,tp);  /* 计算反应物浓度，1e-3用于将密度从kg/m^3转换为g/cm^3 */
     *rr = A_soot_steam_gasification*exp(-E_soot_steam_gasification/Rgas/T_g)* Annealing_soot_steam_gasification * prod * 
            pow(p_h2o, N_soot_steam_gasification)/(1.+K_soot_steam_gasification*p_h2);  /* 计算反应速率，单位 mol/cm^3.s */	 
     *rr *= 1000.; /* 转换为 kmol/(m^3 .s) */     
     *rr = MIN(*rr, rr_turb); // 反应速率取化学反应速率和湍流混合速率的较小值
  }


  
// 定义一个多相非均相反应速率的UDF，用于MGAS模型中的逆水蒸气气化反应
DEFINE_HET_RXN_RATE(SteamGasif_Rev_MGAS,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t); // 获取子线程
     Thread *tp = pt[0]; /* 气相 */
     int index_phase = Get_Phase_Index(hr); // 获取反应的相索引
     Thread *ts = pt[index_phase]; /* 固相 */

     *rr = 0; // 初始化反应速率
     double direction = 0.0, mol_weight, y_carbon, rr_turb = 1e+20; // 初始化方向、摩尔质量、碳质量分数和湍流速率

/*
                              CO + H2 ---> H2O + 1/25 Soot 
			      
    逆水蒸气气化反应，即 CO + H2 ---> 1/25 Soot + H2O，被写成上面的反应形式。
    因此，负的速率意味着CO和H2被消耗，而H2O和Soot被生成。注意，没有C(s)生成，
    上述反应中C(s)的化学计量系数为零。

    设置相和组分的索引。灰分(Ash)组分的索引被初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只被执行一次。这是由FLUENT GUI中
    非均相反应面板中定义的第一个反应完成的。
  */
     
     if(IS_ASH == 0) // 如果灰分组分索引未设置
      SetSpeciesIndex(); // 设置所有组分的索引
     if(MGAS_Gasif) // 如果启用了MGAS气化模型
       {	   
        double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp); // 计算 Ro * R * T
        double p_h2o = RoRT * yi[IP_H2O][IS_H2O]/mw[IP_H2O][IS_H2O]/ 101325.; // 计算H2O的分压（atm）
        double p_co = RoRT * yi[IP_CO][IS_CO]/mw[IP_CO][IS_CO] / 101325.; // 计算CO的分压（atm）
        double p_h2 = RoRT * yi[IP_H2][IS_H2]/mw[IP_H2][IS_H2] / 101325.; // 计算H2的分压（atm）
 
	y_carbon = yi[IP_SOOT][IS_SOOT]; // 反应物是烟灰
	mol_weight = mw[IP_SOOT][IS_SOOT]; // 烟灰的摩尔质量
	
     if (rp_ke) // 如果启用了k-epsilon湍流模型相关的反应模型
        rr_turb = Turbulent_rr(c, t, hr, yi); // 计算湍流混合限制的反应速率
 
        if(C_VOF(c, ts) >= eps_s_small) // 如果固体相的体积分数大于一个很小的值
         {
           *rr = rr_steam_gasif(c, t, ts, tp, p_h2o, p_co, p_h2, y_carbon, mol_weight, &direction); /* 计算水蒸气气化反应速率，单位 mol/(cm^3 .s) */

           if( direction > 0.0)  /* 正值表示 C(s) + H2O ---> CO + H2 */
            *rr = 0.0; // 逆反应速率设为0
	   else                  /* 负值表示 CO + H2 ---> H2O + 1/25 Soot */
	    {
	     *rr = abs(*rr); // 取速率的绝对值
	     *rr = MIN(*rr, rr_turb); // 反应速率取化学反应速率和湍流混合速率的较小值
	    }
         }  
       }
}

// 计算水蒸气气化反应速率的辅助函数
double rr_steam_gasif(cell_t c, Thread *t, Thread *ts, Thread *tp, double p_h2o, double p_co, double p_h2, double y_carbon, double mol_weight, double* direction)
{
     double rate, prod, T_g = MIN((MAX(TMIN,C_T(c,tp))),TMAX); // 获取气体温度，并限制在最小和最大温度之间
     double p_h2o_star = p_h2 * p_co / ( exp(17.29 - 16326/T_g) ); // 计算平衡时的H2O分压

     if(MGAS_Gasif) *direction = p_h2o - p_h2o_star; // MGAS模型中，反应方向由当前H2O分压与平衡分压的差值决定
     if(PCCL_Gasif) *direction = pow(p_h2o, N_steam_gasification)/(1.+K_steam_gasification*p_h2); // PCCL模型中，反应方向（驱动力）的表达式
     
     prod  =  y_carbon*(C_R(c,ts)*1e-03)/mol_weight*C_VOF(c,ts);  /* 计算反应物浓度，1e-3用于将密度从kg/m^3转换为g/cm^3 */
     if(MGAS_Gasif && *direction < 0.0)  /* 这意味着逆向H2O气化 */
       prod = y_carbon*(C_R(c,tp)*1e-03)/mol_weight*C_VOF(c,tp);  /* 逆反应的反应物在气相中，1e-3用于将密度从kg/m^3转换为g/cm^3 */
       
     rate = A_steam_gasification*exp(-E_steam_gasification/Rgas/T_g)* Annealing_steam_gasification * prod * *direction;  /* 计算反应速率，单位 mol/cm^3.s */	 
     rate *= 1000.; /* 转换为 kmol/(m^3 .s) */
     return rate;
}

// 定义一个多相非均相反应速率的UDF，用于计算焦炭的CO2气化
DEFINE_HET_RXN_RATE(Co2Gasif,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t); // 获取子线程
     Thread *tp = pt[0]; /* 气相 */
     int index_phase = Get_Phase_Index(hr); // 获取反应的相索引
     Thread *ts = pt[index_phase]; /* 固相 */

     *rr = 0; // 初始化反应速率
     double direction = 0.0, mol_weight, y_carbon; // 初始化方向、摩尔质量和碳质量分数
  
/*
                              C(s) + CO2 ---> 2CO 
 
    设置相和组分的索引。灰分(Ash)组分的索引被初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只被执行一次。这是由FLUENT GUI中
    非均相反应面板中定义的第一个反应完成的。
  */
     
     if(IS_ASH == 0) // 如果灰分组分索引未设置
       SetSpeciesIndex(); // 设置所有组分的索引
	   
     double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp); // 计算 Ro * R * T
     double p_co = RoRT * yi[IP_CO][IS_CO]/mw[IP_CO][IS_CO] / 101325.; // 计算CO的分压（atm）
     double p_co2 = RoRT * yi[IP_CO2][IS_CO2]/mw[IP_CO2][IS_CO2] / 101325.; // 计算CO2的分压（atm）
       
     SolidFuel_Reactant(c, t, hr, &y_carbon, &mol_weight); // 获取固体燃料反应物（碳）的质量分数和摩尔质量

     if(C_VOF(c, ts) >= eps_s_small) // 如果固体相的体积分数大于一个很小的值
       {
        *rr = rr_co2_gasif(c, t, ts, tp, p_co, p_co2, y_carbon, mol_weight, &direction); /* 计算CO2气化反应速率，单位 mol/(cm^3 .s) */ 
	
        if( direction < 0.0)   /* 负值表示逆向CO2气化 */
         *rr = 0.0; // 正向反应速率设为0
       }	   
}     

// 定义一个多相非均相反应速率的UDF，用于计算烟灰的CO2气化
DEFINE_HET_RXN_RATE(Soot_CO2_Gasif,c,t,r,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t); // 获取子线程
     Thread *tp = pt[0]; /* 气相 */

     *rr = 0; // 初始化反应速率
     double rr_turb = 1e+20; // 初始化湍流反应速率为一个大值
     double T_g = MIN((MAX(TMIN,C_T(c,tp))),TMAX); // 获取气体温度，并限制在最小和最大温度之间
 
/*
                              1/25 Soot + CO2 ---> 2CO  
 
    设置相和组分的索引。灰分(Ash)组分的索引被初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只被执行一次。这是由FLUENT GUI中
    非均相反应面板中定义的第一个反应完成的。
  */
     
     if(IS_ASH == 0) // 如果灰分组分索引未设置
      SetSpeciesIndex(); // 设置所有组分的索引
	   
     double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp); // 计算 Ro * R * T
     double p_co = RoRT * yi[IP_CO][IS_CO]/mw[IP_CO][IS_CO]/ 101325.; // 计算CO的分压（atm）
     double p_co2 = RoRT * yi[IP_CO2][IS_CO2]/mw[IP_CO2][IS_CO2] / 101325.; // 计算CO2的分压（atm）
 
     if (rp_ke) // 如果启用了k-epsilon湍流模型相关的反应模型
       rr_turb = Turbulent_rr(c, t, r, yi); // 计算湍流混合限制的反应速率
       
     double prod  =  yi[IP_SOOT][IS_SOOT]*(C_R(c,tp)*1e-03)/mw[IP_SOOT][IS_SOOT]*C_VOF(c,tp);  /* 计算反应物浓度，1e-3用于将密度从kg/m^3转换为g/cm^3 */
     *rr = A_soot_co2_gasification*exp(-E_soot_co2_gasification/Rgas/T_g)* Annealing_soot_co2_gasification * prod * 
            pow(p_co2, N_soot_co2_gasification)/(1.+K_soot_co2_gasification*p_co);  /* 计算反应速率，单位 mol/cm^3.s */	 
     *rr *= 1000.; /* 转换为 kmol/(m^3 .s) */     
     *rr = MIN(*rr, rr_turb); // 反应速率取化学反应速率和湍流混合速率的较小值
  }

// 定义一个多相非均相反应速率的UDF，用于MGAS模型中的逆CO2气化反应
DEFINE_HET_RXN_RATE(Co2Gasif_Rev_MGAS,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t); // 获取子线程
     Thread *tp = pt[0]; /* 气相 */
     int index_phase = Get_Phase_Index(hr); // 获取反应的相索引
     Thread *ts = pt[index_phase]; /* 固相 */

     *rr = 0; // 初始化反应速率
     double direction = 0.0, mol_weight, y_carbon, rr_turb = 1e+20; // 初始化方向、摩尔质量、碳质量分数和湍流速率
  
/*
                              2CO ---> CO2 + 1/25 Soot  

    逆CO2气化反应，即 2CO ---> 1/25 Soot + CO2，被写成上面的反应形式。
    因此，负的速率意味着CO被消耗，而CO2和Soot被生成。注意，没有C(s)生成，
    上述反应中C(s)的化学计量系数为零。
 
    设置相和组分的索引。灰分(Ash)组分的索引被初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只被执行一次。这是由FLUENT GUI中
    非均相反应面板中定义的第一个反应完成的。
  */
     
     if(IS_ASH == 0) // 如果灰分组分索引未设置
       SetSpeciesIndex(); // 设置所有组分的索引
       
     if(MGAS_Gasif) // 如果启用了MGAS气化模型
       {   
         double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp); // 计算 Ro * R * T
         double p_co = RoRT * yi[IP_CO][IS_CO]/mw[IP_CO][IS_CO] / 101325.; // 计算CO的分压（atm）
         double p_co2 = RoRT * yi[IP_CO2][IS_CO2]/mw[IP_CO2][IS_CO2] / 101325.; // 计算CO2的分压（atm）
     
         y_carbon = yi[IP_SOOT][IS_SOOT]; // 反应物是烟灰
         mol_weight = mw[IP_SOOT][IS_SOOT]; // 烟灰的摩尔质量
       
     if (rp_ke) // 如果启用了k-epsilon湍流模型相关的反应模型
         rr_turb = Turbulent_rr(c, t, hr, yi); // 计算湍流混合限制的反应速率

         if(C_VOF(c, ts) >= eps_s_small) // 如果固体相的体积分数大于一个很小的值
           {
            *rr = rr_co2_gasif(c, t, ts, tp, p_co, p_co2, y_carbon, mol_weight, &direction); /* 计算CO2气化反应速率，单位 mol/(cm^3 .s) */ 
	
            if( direction > 0.0)  /* 正值表示 C(s) + CO2 ---> 2CO */
             *rr = 0.0; // 逆反应速率设为0
	    else                  /* 负值表示 2CO ---> CO2 + 1/25 Soot */
	    {
	     *rr = abs(*rr); // 取速率的绝对值
	     *rr = MIN(*rr, rr_turb); // 反应速率取化学反应速率和湍流混合速率的较小值
	    }
           }
       }
}

// 计算CO2气化反应速率的辅助函数
double rr_co2_gasif(cell_t c, Thread *t, Thread *ts, Thread *tp, double p_co, double p_co2, double y_carbon, double mol_weight, double* direction)
{
     double T_g = MIN(MAX(TMIN,C_T(c,tp)), TMAX), prod; // 获取气体温度并限制范围，声明产物浓度变量
     double p_co2_star = p_co * p_co/(exp(20.92 - 20282/T_g)); // 计算平衡时的CO2分压
 
     if(MGAS_Gasif)  *direction = p_co2-p_co2_star; // MGAS模型中，反应方向由当前CO2分压与平衡分压的差值决定
     if(PCCL_Gasif)  *direction = pow(p_co2, N_co2_gasification)/(1. + K_co2_gasification * p_co); // PCCL模型中，反应方向（驱动力）的表达式
	
     prod  =  y_carbon*C_R(c,ts)*1.e-3/mol_weight* C_VOF(c,ts);  /* 计算反应物浓度，1e-3用于将密度从kg/m^3转换为g/cm^3 */
     if(MGAS_Gasif && *direction < 0.0)  /* 这意味着逆向CO2气化 */
       prod = y_carbon*(C_R(c,tp)*1e-03)/mol_weight*C_VOF(c,tp);  /* 逆反应的反应物在气相中，1e-3用于将密度从kg/m^3转换为g/cm^3 */

     double rate = A_co2_gasification*exp(-E_co2_gasification/Rgas/T_g)*Annealing_co2_gasification * prod * (*direction);  /* 计算反应速率，单位 mol/cm^3.s */
     rate *= 1000.; /* 转换为 kmol/(m^3 .s) */
     return rate;
}     
            
// 定义一个多相非均相反应速率的UDF，用于计算焦炭的H2气化（加氢气化）
DEFINE_HET_RXN_RATE(H2Gasif,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t); // 获取子线程
     Thread *tp = pt[0]; /* 气相 */
     int index_phase = Get_Phase_Index(hr); // 获取反应的相索引
     Thread *ts = pt[index_phase]; /* 固相 */

     *rr = 0; // 初始化反应速率
     double direction = 0.0, mol_weight, y_carbon; // 初始化方向、摩尔质量和碳质量分数
  
/*
                              1/2 C(s) + H2 ---> 1/2 CH4 
 
    设置相和组分的索引。灰分(Ash)组分的索引被初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只被执行一次。这是由FLUENT GUI中
    非均相反应面板中定义的第一个反应完成的。
  */
     
     if(IS_ASH == 0) // 如果灰分组分索引未设置
       SetSpeciesIndex(); // 设置所有组分的索引
	   
     double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp); // 计算 Ro * R * T
     double p_h2 = RoRT * yi[IP_H2][IS_H2]/mw[IP_H2][IS_H2] / 101325.; // 计算H2的分压（atm）
     double p_ch4 = RoRT * yi[IP_CH4][IS_CH4]/mw[IP_CH4][IS_CH4] / 101325.; // 计算CH4的分压（atm）
     
     SolidFuel_Reactant(c, t, hr, &y_carbon, &mol_weight); // 获取固体燃料反应物（碳）的质量分数和摩尔质量
     	   
     if(C_VOF(c, ts) >= eps_s_small) // 如果固体相的体积分数大于一个很小的值
       {
        *rr = rr_h2_gasif(c, t, ts, tp, p_h2, p_ch4, y_carbon, mol_weight, &direction); /* 计算H2气化反应速率，单位 mol/(cm^3 .s) */
	
        if( direction < 0.0)   /* 负值表示逆向H2气化 */
         *rr = 0.0; // 正向反应速率设为0
       }  
} 

// 定义一个多相非均相反应速率的UDF，用于计算烟灰的H2气化（加氢气化）
DEFINE_HET_RXN_RATE(Soot_H2_Gasif,c,t,r,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t); // 获取子线程
     Thread *tp = pt[0]; /* 气相 */

     *rr = 0; // 初始化反应速率
     double rr_turb = 1e+20; // 初始化湍流反应速率为一个大值
     double T_g = MIN((MAX(TMIN,C_T(c,tp))),TMAX); // 获取气体温度，并限制在最小和最大温度之间
 
/*
                              1/25 Soot + 2H2 ---> CH4 
 
    设置相和组分的索引。灰分(Ash)组分的索引被初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只被执行一次。这是由FLUENT GUI中
    非均相反应面板中定义的第一个反应完成的。
  */
     
     if(IS_ASH == 0) // 如果灰分组分索引未设置
      SetSpeciesIndex(); // 设置所有组分的索引
	   
     double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp); // 计算 Ro * R * T
     double p_h2 = RoRT * yi[IP_H2][IS_H2]/mw[IP_H2][IS_H2] / 101325.; // 计算H2的分压（atm）
 
     if (rp_ke) // 如果启用了k-epsilon湍流模型相关的反应模型
       rr_turb = Turbulent_rr(c, t, r, yi); // 计算湍流混合限制的反应速率
       
     double prod  =  yi[IP_SOOT][IS_SOOT]*(C_R(c,tp)*1e-03)/mw[IP_SOOT][IS_SOOT]*C_VOF(c,tp);  /* 计算反应物浓度，1e-3用于将密度从kg/m^3转换为g/cm^3 */
     *rr = A_soot_h2_gasification*exp(-E_soot_h2_gasification/Rgas/T_g)* Annealing_soot_h2_gasification * prod * 
            pow(p_h2, N_soot_h2_gasification);  /* 计算反应速率，单位 mol/cm^3.s */	 
     *rr *= 1000.; /* 转换为 kmol/(m^3 .s) */     
     *rr = MIN(*rr, rr_turb); // 反应速率取化学反应速率和湍流混合速率的较小值
  }
  
// 定义一个多相非均相反应速率的UDF，用于MGAS模型中的逆H2气化反应
DEFINE_HET_RXN_RATE(H2Gasif_Rev_MGAS,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t); // 获取子线程
     Thread *tp = pt[0]; /* 气相 */
     int index_phase = Get_Phase_Index(hr); // 获取反应的相索引
     Thread *ts = pt[index_phase]; /* 固相 */

     *rr = 0; // 初始化反应速率
     double direction = 0.0, mol_weight, y_carbon, rr_turb = 1e+20; // 初始化方向、摩尔质量、碳质量分数和湍流速率
  
/*
                              1/2 CH4 ---> H2 + (0.5)*1/25 Soot 
			      
    逆H2气化反应，即 1/2 CH4 ---> 1/25 Soot + H2，被写成上面的反应形式。
    因此，负的速率意味着CH4被消耗，而H2和Soot被生成。注意，没有C(s)生成，
    上述反应中C(s)的化学计量系数为零。
 
    设置相和组分的索引。灰分(Ash)组分的索引被初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只被执行一次。这是由FLUENT GUI中
    非均相反应面板中定义的第一个反应完成的。
  */
     
     if(IS_ASH == 0) // 如果灰分组分索引未设置
       SetSpeciesIndex(); // 设置所有组分的索引
     if(MGAS_Gasif) // 如果启用了MGAS气化模型
       {    
        double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp); // 计算 Ro * R * T
        double p_h2 = RoRT * yi[IP_H2][IS_H2]/mw[IP_H2][IS_H2] / 101325.; // 计算H2的分压（atm）
        double p_ch4 = RoRT * yi[IP_CH4][IS_CH4]/mw[IP_CH4][IS_CH4] / 101325.; // 计算CH4的分压（atm）

        y_carbon = yi[IP_SOOT][IS_SOOT]; // 反应物是烟灰
        mol_weight = mw[IP_SOOT][IS_SOOT]; // 烟灰的摩尔质量
     	   
        if(C_VOF(c, ts) >= eps_s_small) // 如果固体相的体积分数大于一个很小的值
          {
            if (rp_ke) // 如果启用了k-epsilon湍流模型相关的反应模型
            rr_turb = Turbulent_rr(c, t, hr, yi); // 计算湍流混合限制的反应速率
	    
           *rr = rr_h2_gasif(c, t, ts, tp, p_h2, p_ch4, y_carbon, mol_weight, &direction); /* 计算H2气化反应速率，单位 mol/(cm^3 .s) */
	
           if( direction > 0.0)  /* 正值表示 1/2 C(s) + H2 ---> 1/2 CH4 */
            *rr = 0.0; // 逆反应速率设为0
	   else                  /* 负值表示 1/2 CH4 ---> H2 + (0.5)*1/25 Soot */
	    {
	     *rr = abs(*rr); // 取速率的绝对值
	     *rr = MIN(*rr, rr_turb); // 反应速率取化学反应速率和湍流混合速率的较小值
	    } 
          }  
       }
}

// 计算H2气化反应速率的辅助函数
double rr_h2_gasif(cell_t c, Thread *t, Thread *ts, Thread *tp, double p_h2, double p_ch4, double y_carbon, double mol_weight, double* direction)
{
     double rate = 0.0, prod; // 初始化速率和产物浓度变量
     double T_g = MIN((MAX(TMIN,C_T(c,tp))), TMAX); // 获取气体温度，并限制在最小和最大温度之间
     double p_h2_star = pow ((p_ch4/(exp(-13.43 + 10999/T_g))), 0.5); // 计算平衡时的H2分压

     prod  =  y_carbon*C_R(c,ts)*1.e-3/mol_weight * C_VOF(c,ts);  /* 计算反应物浓度，1e-3用于将密度从kg/m^3转换为g/cm^3 */
 
     if(MGAS_Gasif) // 如果启用了MGAS气化模型
       { 
         *direction = p_h2-p_h2_star; // 反应方向由当前H2分压与平衡分压的差值决定
         if(*direction < 0.0)  /* 这意味着逆向H2气化 */
           prod = y_carbon*(C_R(c,tp)*1e-03)/mol_weight*C_VOF(c,tp);  /* 逆反应的反应物在气相中，1e-3用于将密度从kg/m^3转换为g/cm^3 */
         rate = exp( -7.087 - 8078/T_g )* prod * *direction ;  /* 计算反应速率，单位 mol/cm^3.s */
       }
     if(PCCL_Gasif) // 如果启用了PCCL气化模型
       {
        *direction = p_h2; // 反应驱动力为H2分压
        rate = A_h2_gasification*exp(-E_h2_gasification/Rgas/T_g)*Annealing_h2_gasification * prod * *direction;  /* 计算反应速率，单位 mol/cm^3.s */
       }

     rate *= 1000.; /* 转换为 kmol/(m^3 .s) */
     return rate;
}

// 定义一个多相非均相反应速率的UDF，用于计算煤的燃烧
DEFINE_HET_RXN_RATE(coal_combustion,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t); // 获取子线程
     Thread *tp = pt[0]; /* 气相 */
     int index_phase = Get_Phase_Index(hr); // 获取反应的相索引
     Thread *ts = pt[index_phase]; /* 固相 */
     double mol_weight, y_carbon, y_ash; // 摩尔质量，碳质量分数，灰分质量分数
     *rr = 0.0; // 初始化反应速率
         
 /* 设置相和组分的索引。灰分(Ash)组分的索引被初始化为零，以及所有其他索引。
    灰分组分索引被用作一个标志，以确保SetSpeciesIndex只被执行一次。这是由FLUENT GUI中
    非均相反应面板中定义的第一个反应完成的。
  */
     if(IS_ASH == 0) // 如果灰分组分索引未设置
       SetSpeciesIndex(); // 设置所有组分的索引
       
     if( C_YI(c,tp,IS_O2) >= spe_small) // 如果氧气的质量分数大于一个很小的值
      {
        SolidFuel_Reactant(c, t, hr, &y_carbon, &mol_weight); // 获取固体燃料反应物（碳）的质量分数和摩尔质量
	y_ash = yi[index_phase][IS_ASH]; // 获取灰分的质量分数

        *rr = rr_combustion(c, t, ts, tp, yi[IP_O2][IS_O2], y_ash, y_carbon); /* 计算燃烧反应速率，单位 mol/(cm^3 .s) */
        *rr *=  1000.; /* 转换为 kmol/(m^3 .s) */ 
      } 
}

// 计算燃烧反应速率的辅助函数
double rr_combustion(cell_t c, Thread *t, Thread *ts, Thread *tp, double yi_O2, double y_ash, double y_carbon)
{		   
		   
     double rd, k_f, k_r, factor, k_a, rate = 0.0, vrel;		   
     double Pt = MAX(0.1, (op_pres+C_P(c,t))/101325); // 计算总压（atm）
     double gas_constant = 82.06; /* 气体常数，单位 atm.cm^3/mol.K */
     double T = C_T(c,tp), T_s = C_T(c,ts), D_p = C_PHASE_DIAMETER(c,ts)*100.; // 获取气相温度、固相温度和颗粒直径（cm）
     double p_o2 = C_R(c,tp)*UNIVERSAL_GAS_CONSTANT* T *yi_O2/mw[IP_O2][IS_O2] / 101325.; /* 计算氧气分压（atm） */

     if(fc_ar > 0.) // 如果固定碳含量大于0
       {
        if (y_ash > 0.) // 如果灰分含量大于0
          {
            rd = pow( (y_carbon * ash_ar/100.)/(y_ash * fc_ar/100.), (1./3.) ); // 计算收缩核模型的无量纲半径
            rd = MIN(1., rd); // 半径不能大于1
          }
        else 
            rd = 1.;  
       }
     else 
        rd = 0.;
       
     double diff = MAX((4.26 * pow((T/1800.),1.75)/Pt), 1.e-10); /* 计算扩散系数 (cm^2/s) */
     double Sc1o3 = pow(C_MU_L(c,tp)/(C_R(c,tp) * diff * 1.e-4), 1./3.); // 计算施密特数的1/3次方
#if RP_2D
     vrel = pow(( (C_U(c,tp)-C_U(c,ts))*(C_U(c,tp)-C_U(c,ts)) +
            (C_V(c,tp)-C_V(c,ts))*(C_V(c,tp)-C_V(c,ts))), 0.5);  // 计算二维相对速度
#endif
#if RP_3D
     vrel = pow(( (C_U(c,tp)-C_U(c,ts))*(C_U(c,tp)-C_U(c,ts)) +
            (C_V(c,tp)-C_V(c,ts))*(C_V(c,tp)-C_V(c,ts)) +  
            (C_W(c,tp)-C_W(c,ts))*(C_W(c,tp)-C_W(c,ts)) ), 0.5);  // 计算三维相对速度
#endif
     double Re = C_VOF(c,tp) * D_p/100. * vrel * C_R(c,tp)/(C_MU_L(c,tp)+SMALL_S);  // 计算雷诺数
     double N_sherwood = (7. - 10. * C_VOF(c,tp) + 5. * C_VOF(c,tp) * C_VOF(c,tp) )*
                         (1. + 0.7 * pow(Re, 0.2) * Sc1o3) + 
                         (1.33 - 2.4 * C_VOF(c,tp) + 1.2 * C_VOF(c,tp) * C_VOF(c,tp)) *
                         pow(Re, 0.7) * Sc1o3;   // 计算舍伍德数

     if ( rd <= 0. || C_VOF(c, ts) <= 0. ) // 如果无量纲半径或固相体积分数小于等于0
       {
          rate = 0.; // 反应速率为0
       }
     else
       {   
          k_f = diff * N_sherwood / (D_p * gas_constant/mw[IP_O2][IS_O2] * T ); /* 外部传质系数 g/(atm.cm^2.s) */
          k_r = A_c_combustion * exp( -E_c_combustion/Rgas/T_s ) * rd * rd; // 表面化学反应速率系数
          if ( rd >= 1.) // 如果是初始燃烧阶段
           {
            rate = 1. / (1./k_f + 1./k_r); // 速率由外部传质和化学反应控制
           }
          else
           {
            k_a = 2. * rd * diff * f_ep_a / (D_p * (1.-rd) * gas_constant/mw[IP_O2][IS_O2] * T_s );  // 灰层扩散系数
            rate = 1. / (1./k_f + 1./k_r + 1./k_a); // 速率由外部传质、化学反应和灰层扩散共同控制
           }

          factor = y_carbon / (y_carbon + 1.e-6); // 碳质量分数因子
          rate *= p_o2 * 6. * C_VOF(c,ts) * factor / (D_p * 32.); /* 计算最终反应速率 mol/(cm^3 .s) */ 
       }  	   
      return rate;
}

#if !RP_NODE || !PARALLEL // 仅在串行或非计算节点上执行

// 读取挥发分质量分数和其他参数
void volatile_mass_fractions()
{            
       read_c3m_data(); // 读取c3m数据
       
       /* pan2 : Oct 2012 ... 为调试添加了CX_Messages */
       // 打印各种反应模型的布尔标志值
       CX_Message("PCCL_Devol = %d\n",PCCL_Devol);
       CX_Message("MGAS_Devol = %d\n",MGAS_Devol);
       // ... (省略类似的打印语句)
       
       // 打印煤的工业分析数据
       CX_Message("fc_ar = %f\n",fc_ar);
       CX_Message("vm_ar = %f\n",vm_ar);
       CX_Message("ash_ar = %f\n",ash_ar);
       CX_Message("moist_ar = %f\n",moist_ar);

       // 打印各种反应的动力学参数
       CX_Message("a1_devolatilization = %f\n",A1_devolatilization);
       CX_Message("e1_devolatilization = %f\n",E1_devolatilization);
       // ... (省略类似的打印语句)
       
       
  /* f_ep_a 用于收缩核模型，在煤燃烧模型中 */
       double ep_a = 0.25 + 0.75*(1-ash_ar/100.); // 计算灰分孔隙率
       f_ep_a = pow(ep_a,2.5); // 计算孔隙率因子
}       


// 设置布尔值的静态函数
static void SetBooleanValue(char * var , char * svalue)  /* pan : Oct 2012 ... new function */
{
	cxboolean value;

	if ( strcmp(svalue,"true") == 0)  /* pan2 : Oct 2012 : correction */
		value = TRUE;
	else
		value = FALSE;
	
	// 根据变量名设置相应的全局布尔变量
    if (strcmp(var,"pccl_devol")               == 0)  PCCL_Devol               = value;
    if (strcmp(var,"mgas_devol")               == 0)  MGAS_Devol               = value;
    // ... (省略类似的比较和赋值语句)
}

// 设置浮点数值的静态函数
static void SetValue(char * var , char * svalue) /* pan : oct 2012 ... replace entire function */
{
	char * pEnd;

	double value = strtod(svalue,&pEnd); // 将字符串转换为double
	

	// 根据变量名设置相应的全局浮点变量
	if (strcmp(var,"fc_ar")     == 0) fc_ar    = value;
	if (strcmp(var,"vm_ar")     == 0) vm_ar    = value;
	// ... (省略类似的比较和赋值语句)
}

// 从文件读取c3m数据的函数
void read_c3m_data()
{
    CX_Message("start of read_c3m_data \n"); // 打印开始信息

	FILE * pFile; // 文件指针
	char line[80]; // 行缓冲区
	char field1[80]; // 字段1缓冲区
	char field2[80]; // 字段2缓冲区
	char *pch; // 字符串指针
	int  i , field_index; // 循环变量和字段索引

	pFile = fopen("fluent_c3m_udf.inp","r"); // 打开输入文件

	if (pFile != NULL) // 如果文件成功打开
	{
	    CX_Message("fopen OK \n"); // 打印成功信息
		char * pEnd = line;
		
		while (pEnd != NULL) // 循环读取每一行
		{
                        /* CX_Message("read a line - start \n"); */
			pEnd = fgets(line,80,pFile); // 读取一行

			if (pEnd != NULL)
			{
				for (i=0; i<80; ++i) line[i] = tolower(line[i]); // 将行内容转为小写

				pch = strtok(line," ,\t\n="); // 使用分隔符分割字符串

				field_index = 0;
				while (pch != NULL)
				{
					if (field_index == 0)
					{
						strcpy(field1,pch); // 复制第一个字段（变量名）
						field_index = 1;
                                             /*   CX_Message("token %s \n",field1); */
					}
					else
					{
						strcpy(field2,pch); // 复制第二个字段（值）
						field_index = 0;
						SetValue(field1,field2); // 设置浮点数值
						SetBooleanValue(field1,field2);   /* pan : Oct 2012 */ // 设置布尔值
						
						/*
						 CX_Message("\n\n");
                                                 CX_Message("token %s \n",field1); 
                                                 CX_Message("token %s \n",field2); 
						 */
					}
					pch = strtok(NULL," ,\t\n="); // 获取下一个token
				}
			}
		}

		fclose(pFile); // 关闭文件
	}
}
/* pan c3m end */
#endif


// 定义一个交换属性UDF，用于计算煤颗粒与气相之间的传热
DEFINE_EXCHANGE_PROPERTY(Heat_Trans_Coal, c, t, i, j)
{
  Thread *ti = THREAD_SUB_THREAD(t,i); // 获取相i的线程
  Thread *tj = THREAD_SUB_THREAD(t,j); // 获取相j的线程
  double val;

  val = heat_gunn_udf(c,ti, tj); // 调用Gunn模型的传热计算函数
  return val;
}

// 定义一个交换属性UDF，用于计算再循环颗粒与气相之间的传热
DEFINE_EXCHANGE_PROPERTY(Heat_Trans_Recy, c, t, i, j)
{
  Thread *ti = THREAD_SUB_THREAD(t,i); // 获取相i的线程
  Thread *tj = THREAD_SUB_THREAD(t,j); // 获取相j的线程
  double val;

  val = heat_gunn_udf(c,ti, tj); // 调用Gunn模型的传热计算函数
  return val;
}

// Gunn模型的传热计算函数
double heat_gunn_udf(cell_t c, Thread *ti, Thread *tj) 
{                       
  double h;               
  double d = C_PHASE_DIAMETER(c,tj); // 颗粒直径
  double k = C_K_L(c,ti); // 气相导热系数
  double vf = C_VOF(c,ti); // 气相体积分数
  double vf2 = vf*vf;
  double vel, Re, Pr, Nu; // 相对速度，雷诺数，普朗特数，努塞尔数
      
#if RP_2D 
    vel = pow(((C_U(c,tj)-C_U(c,ti))*(C_U(c,tj)-C_U(c,ti)) + (C_V(c,tj)-C_V(c,ti))*(C_V(c,tj)-C_V(c,ti))),0.5); // 二维相对速度
#endif 

#if RP_3D 
    vel = pow(((C_U(c,tj)-C_U(c,ti))*(C_U(c,tj)-C_U(c,ti)) + (C_V(c,tj)-C_V(c,ti))*(C_V(c,tj)-C_V(c,ti)) + 
          (C_W(c,tj)-C_W(c,ti))*(C_W(c,tj)-C_W(c,ti))),0.5); // 三维相对速度
#endif 
    
  Re = RE_NUMBER(C_R(c,ti),vel,d,C_MU_L(c,ti)); // 计算雷诺数
  Pr = PR_NUMBER (C_CP(c,ti),C_MU_L(c,ti),k); // 计算普朗特数
  Pr = pow (Pr,1./3.);
  Nu = (7. - 10*vf + 5.*vf2)*(1. + 0.7*pow(Re,0.2)*Pr) +
    (1.33 - 2.4*vf + 1.2*vf2)*pow(Re,0.7)*Pr; // 计算努塞尔数
                            
  h =IP_HEAT_COEFF(C_VOF(c,ti),C_VOF(c,tj),k,Nu,d); // 计算相间传热系数
  return h;
}
 
// 计算饱和蒸汽压的函数
double satPressure(double T)
{
  const double Tstd = 273.15; // 标准温度（K）
  double TT = T - Tstd; // 摄氏温度
  double p_sat;

    p_sat = (100.0 * 6.1121 * exp((18.678-TT/234.5)*(TT/(257.14+TT)))); // 计算饱和蒸汽压（Pa）
  return p_sat;
}

// 获取非均相反应中次要相的索引
double Get_Phase_Index(Hetero_Reaction *hr)
{
  /* 这个例程返回以下反应的次要相的相索引
    
         A(s) + B(g) ---> ....
	 
     如果反应物中有多个次要相，下面获得的相索引将始终是
     反应物列表中最后一个次要组分的相索引 */	 

     Domain **dr = hr->domain_reactant; // 获取反应物域
     int i, index_phase = 0, iphase;
     for (i = 0; i < hr->n_reactants; i++) // 遍历所有反应物
     {
       iphase = DOMAIN_INDEX(dr[i]); // 获取当前反应物的相索引
       if(iphase > 0) // 如果不是主相（气相）
         index_phase = iphase; // 更新相索引
     } 
     return index_phase;    
}

// 计算湍流混合限制的反应速率
double Turbulent_rr(cell_t c, Thread *t, Hetero_Reaction *r, real yi[MAX_PHASES][MAX_SPE_EQNS])
{
  Thread **pt = THREAD_SUB_THREADS(t);
  Thread *tp = pt[0]; /* 气相 */
  double minR = 1.e+20, numerator = 0.;
  double denom =0.,ci;

  int i;
  double Amix = 4., Bmix = 0.5; // Eddy-Dissipation模型常数
	   
  for(i=0; i< r->n_reactants; ++i) // 遍历反应物
    {
      if(r->stoich_reactant[i] != 0.0)
        {
          ci = C_R(c,tp)*yi[0][r->reactant[i]]/mw[0][r->reactant[i]]; // 计算摩尔浓度
          SETMIN(minR, ci/r->stoich_reactant[i]);	   // 找到限制反应的组分
	}
     }
    
	     
  for(i=0; i< r->n_products; ++i) // 遍历产物
    {
      if(r->stoich_product[i] != 0.0)
        {
          ci = C_R(c,tp)*yi[0][r->product[i]]/mw[0][r->product[i]]; // 计算摩尔浓度
          numerator +=ci*mw[0][r->product[i]];
          denom +=r->stoich_product[i]*mw[0][r->product[i]];
        }
    }
    	     
  double rr_turb = Amix * C_D(c,tp)/C_K(c,tp) * MIN(minR, Bmix*numerator/denom);  // 计算湍流反应速率
  return rr_turb;
}

// 获取固体燃料反应物的信息
void SolidFuel_Reactant(cell_t c, Thread *t, Hetero_Reaction *hr, double* y_carbon, double* mol_weight)
{
  Thread **pt = THREAD_SUB_THREADS(t);
  int i;
	   
  for(i=0; i< hr->n_reactants; ++i) // 遍历反应物
    {
      if(hr->stoich_reactant[i] != 0.0)
        {
	  int phase = DOMAIN_INDEX(hr->domain_reactant[i]); // 获取相索引
          Thread *thread = pt[phase]; /* 固相 */ 
	  if(phase != 0) // 如果是固相
	   {
	     *y_carbon = C_YI(c,thread,hr->reactant[i]); // 获取碳的质量分数
	     *mol_weight = mw[phase][hr->reactant[i]]; // 获取摩尔质量
	   }
        }
    }
}