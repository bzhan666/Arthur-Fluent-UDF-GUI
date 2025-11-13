/*
Fluent-UDF_煤炭燃烧多相流过程
---各种非均相反应的代码
按照不同的物理化学过程(如水分释放、脱挥发分、焦炭燃烧、气化等)组织成
多个独立的 DEFINE_HET_RXN_RATE 函数,每个函数负责一个特定的反应。
考虑了煤燃烧过程中的多个复杂现象,包括水分蒸发、脱挥发分、焦炭与多种气体
(O2, H2O, CO2, H2)的非均相反应、以及挥发物(焦油、CH4等)的均相燃烧。

修改记录 (2024):
- 补全了 SetValue() 和 SetBooleanValue() 函数
- 改进了错误处理
- 优化了代码结构
*/

#include "udf.h"
#include "stdio.h"
#include "time.h"

#define SMALL_S 1.e-29 // 一个很小的数,用于避免除以零
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

/* 函数声明 */
void SolidFuel_Reactant(cell_t c, Thread *t, Hetero_Reaction *hr, double* y_carbon, double* mol_weight);
void volatile_mass_fractions();
void SetSpeciesIndex();
void read_c3m_data();
double satPressure(double T);
double Get_Phase_Index(Hetero_Reaction *hr);
double Turbulent_rr(cell_t c, Thread *t, Hetero_Reaction *r, real yi[MAX_PHASES][MAX_SPE_EQNS]);
double Mass_Transfer_Coeff(cell_t c, Thread *tp, Thread *ts);
double heat_gunn_udf(cell_t c, Thread *ti, Thread *tj);
double rr_combustion(cell_t c, Thread *t, Thread *ts, Thread *tp, double yi_O2, double y_ash, double y_carbon);
double rr_steam_gasif(cell_t c, Thread *t, Thread *ts, Thread *tp, double p_h2o, double p_co, double p_h2, double y_carbon, double mol_weight, double* direction);
double rr_co2_gasif(cell_t c, Thread *t, Thread *ts, Thread *tp, double y_co, double y_co2, double y_carbon, double mol_weight, double* direction);
double rr_h2_gasif(cell_t c, Thread *t, Thread *ts, Thread *tp, double y_h2, double y_ch4, double y_carbon, double mol_weight, double* direction);

/* 组分索引变量 (IP: 相索引, IS: 组分索引) */
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

int current_c3m_solid_phase = 0;

/* 模型开关 */
cxboolean init_flag                = TRUE;
cxboolean PCCL_Devol               = FALSE;
cxboolean MGAS_Devol               = FALSE;
cxboolean CPD_Devol                = FALSE;
cxboolean FGDVC_Devol              = FALSE;
cxboolean HPTR_Devol               = FALSE;
cxboolean MGAS_Moisture            = FALSE;
cxboolean PCCL_Moisture            = FALSE;
cxboolean MGAS_TarCracking         = FALSE;
cxboolean PCCL_2nd_Pyro            = FALSE;
cxboolean MGAS_Gasif               = FALSE;
cxboolean PCCL_Gasif               = FALSE;
cxboolean PCCL_TarCracking         = FALSE;
cxboolean MGAS_WGS                 = FALSE;
cxboolean PCCL_soot_gasif          = FALSE;
cxboolean MGAS_char_combustion     = FALSE;
cxboolean PCCL_char_combustion     = FALSE;
cxboolean PCCL_soot_oxidation      = FALSE;
cxboolean TAR_oxidation            = FALSE;
cxboolean MGAS_gas_phase_oxidation = FALSE;

double avg_mf_h2o,avg_mf_co,avg_mf_h2,avg_mf_ch4,avg_mf_co2,avg_c,avg_volatile,avg_moisture,avg_ash;

/* 煤分析变量 */
double fc_ar=0.,vm_ar=0.,ash_ar=0.,moist_ar=0.;
double f_ep_a = 0.;
double mw[MAX_PHASES][MAX_SPE_EQNS];

/* 动力学参数 */
double A1_devolatilization=0.0, E1_devolatilization=0.0;
double A2_devolatilization=0.0, E2_devolatilization=0.0;
double A_tar_cracking=0.0 , E_tar_cracking=0.0;
double A_steam_gasification=0.0, E_steam_gasification=0.0;
double K_steam_gasification=0.0, N_steam_gasification=0.0;
double Annealing_steam_gasification=1.0;
double A_co2_gasification=0.0, E_co2_gasification=0.0;
double K_co2_gasification=0.0, N_co2_gasification=0.0;
double Annealing_co2_gasification=1.0;
double A_h2_gasification=0.0, E_h2_gasification=0.0;
double N_h2_gasification=0.0;
double Annealing_h2_gasification=1.0;
double A_soot_steam_gasification=0.0, E_soot_steam_gasification=0.0;
double K_soot_steam_gasification=0.0, N_soot_steam_gasification=0.0;
double Annealing_soot_steam_gasification=0.0;
double A_soot_co2_gasification=0.0, E_soot_co2_gasification=0.0;
double K_soot_co2_gasification=0.0, N_soot_co2_gasification=0.0;
double Annealing_soot_co2_gasification=0.0;
double A_soot_h2_gasification=0.0, E_soot_h2_gasification=0.0;
double N_soot_h2_gasification=0.0;
double Annealing_soot_h2_gasification=0.0;
double A_Soot_Combustion = 0.0, E_Soot_Combustion = 0.0;
double A_c_combustion = 8710., E_c_combustion = 27000.;
double Annealing_c_combustion=0.0 , N_c_combustion=0.0;
double A_moisture_release = 0.0, E_moisture_release = 0.0;
double wg3 = 0.014;
double Moisture_Flux;

/* DEFINE_ADJUST - 初始化函数 */
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
    host_to_node_real_5(fc_ar,vm_ar,ash_ar,moist_ar,f_ep_a);
    host_to_node_int_4(PCCL_Devol,MGAS_Devol,MGAS_Moisture,PCCL_2nd_Pyro);
    host_to_node_int_3(MGAS_Gasif,PCCL_Gasif,PCCL_Moisture);
    host_to_node_int_3(CPD_Devol,FGDVC_Devol,HPTR_Devol);
    host_to_node_int_3(PCCL_TarCracking,MGAS_WGS,PCCL_soot_gasif);
    host_to_node_int_1(PCCL_char_combustion);
    host_to_node_int_3(MGAS_char_combustion,TAR_oxidation,MGAS_gas_phase_oxidation);
    host_to_node_real_2(A_moisture_release, E_moisture_release);
    host_to_node_real_2(A1_devolatilization, E1_devolatilization);
    host_to_node_real_2(A2_devolatilization, E2_devolatilization);
    host_to_node_real_2(A_tar_cracking, E_tar_cracking);
    host_to_node_real_2(A_steam_gasification, E_steam_gasification);
    host_to_node_real_2(K_steam_gasification, N_steam_gasification);
    host_to_node_real_1(Annealing_steam_gasification);
    host_to_node_real_2(A_co2_gasification, E_co2_gasification);
    host_to_node_real_2(K_co2_gasification, N_co2_gasification);
    host_to_node_real_1(Annealing_co2_gasification);
    host_to_node_real_2(A_h2_gasification, E_h2_gasification);
    host_to_node_real_2(Annealing_h2_gasification, N_h2_gasification);
    host_to_node_real_2(A_soot_steam_gasification, E_soot_steam_gasification);
    host_to_node_real_2(K_soot_steam_gasification, N_soot_steam_gasification);
    host_to_node_real_1(Annealing_soot_steam_gasification);
    host_to_node_real_2(A_soot_co2_gasification, E_soot_co2_gasification);
    host_to_node_real_2(K_soot_co2_gasification, N_soot_co2_gasification);
    host_to_node_real_1(Annealing_soot_co2_gasification);
    host_to_node_real_2(A_soot_h2_gasification, E_soot_h2_gasification);
    host_to_node_real_2(Annealing_soot_h2_gasification, N_soot_h2_gasification);
    host_to_node_real_2(A_Soot_Combustion, E_Soot_Combustion);
    host_to_node_real_2(A_c_combustion, E_c_combustion);
    host_to_node_real_2(N_c_combustion, Annealing_c_combustion);
    host_to_node_real_2(wg3,Moisture_Flux);
    init_flag = FALSE;
  }
}

DEFINE_ON_DEMAND(Devol_and_Tar_Cracking)
{
#if !RP_NODE
    volatile_mass_fractions();
#endif
}

/* SetSpeciesIndex - 设置组分索引 */
void SetSpeciesIndex()
{
  Domain *domain = Get_Domain(1);
  int n, ns, MAX_SPE_EQNS_PRIM = 0;
  Domain *subdomain;

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
          IP_CH4 = n; IS_CH4 = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"co"))
        {
          IP_CO = n; IS_CO = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"c2h2"))
        {
          IP_C2H2 = n; IS_C2H2 = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"c2h4"))
        {
          IP_C2H4 = n; IS_C2H4 = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"c2h6"))
        {
          IP_C2H6 = n; IS_C2H6 = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"c3h6"))
        {
          IP_C3H6 = n; IS_C3H6 = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"c3h8"))
        {
          IP_C3H8 = n; IS_C3H8 = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"co2"))
        {
          IP_CO2 = n; IS_CO2 = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"h2"))
        {
          IP_H2 = n; IS_H2 = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"h2o"))
        {
          IP_H2O = n; IS_H2O = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"o2"))
        {
          IP_O2 = n; IS_O2 = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"h2s"))
        {
          IP_H2S = n; IS_H2S = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"cl2"))
        {
          IP_CL2 = n; IS_CL2 = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"nh3"))
        {
          IP_NH3 = n; IS_NH3 = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"n2"))
        {
          IP_N2 = n; IS_N2 = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"oil"))
        {
          IP_OIL = n; IS_OIL = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"pah"))
        {
          IP_PAH = n; IS_PAH = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"tar"))
        {
          IP_TAR = n; IS_TAR = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"c"))
        {
          IP_C = n; IS_C = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"c_recycle"))
        {
          IP_C_R = n; IS_C_R = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"soot"))
        {
          IP_SOOT = n; IS_SOOT = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"volatile"))
        {
          IP_VOL = n; IS_VOL = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"h2o<l>"))
        {
          IP_MOISTURE = n; IS_MOISTURE = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"h2o<l>-slurry"))
        {
          IP_SLURRY = n; IS_SLURRY = ns;
          mw[n][ns] = MATERIAL_PROP(s_mat,PROP_mwi);
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"h2o<l>-dummy"))
        {
          IP_SLURRY_D = n; IS_SLURRY_D = ns;
          mw[n][ns] = MATERIAL_PROP(s_mat,PROP_mwi);
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"ash-coal"))
        {
          IP_ASH = n; IS_ASH = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"ash-dummy"))
        {
          IP_ASH_D = n; IS_ASH_D = ns;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"ash-recycle"))
        {
          IP_ASH_R = n; IS_ASH_R = ns;
          if(n == 0) MAX_SPE_EQNS_PRIM +=1;
        }
        else if (0 == strcmp(MIXTURE_SPECIE_NAME(m_mat,ns),"si<s>"))
        {
          IP_SAND = n; IS_SAND = ns;
        }

        mw[n][ns] = MATERIAL_PROP(s_mat,PROP_mwi);
      }
    }
  }
}

/* ==================== 反应速率定义 ==================== */

// 定义水分释放的非均相反应速率
DEFINE_HET_RXN_RATE(moisture_release,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t); // 获取子线程
     int index_phase = Get_Phase_Index(hr); // 获取相索引
     Thread *ts = pt[index_phase]; /* 固体相 */
     double prod = 0.0, Ts = C_T(c,ts); // 产物,固体温度
     double Pt = MAX(0.1,(op_pres+C_P(c,t))); // 总压
     double Tsat = 1./(0.0727/log(Pt/611.) - 0.0042) + 273.; // 饱和温度

     *rr = 0;

 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零,以及所有其他索引。
    灰分组分索引被用作一个标志,以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
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

 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零,以及所有其他索引。
    灰分组分索引被用作一个标志,以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
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

 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零,以及所有其他索引。
    灰分组分索引被用作一个标志,以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
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

 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零,以及所有其他索引。
    灰分组分索引被用作一个标志,以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
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

 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零,以及所有其他索引。
    灰分组分索引被用作一个标志,以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
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

 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零,以及所有其他索引。
    灰分组分索引被用作一个标志,以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
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

 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零,以及所有其他索引。
    灰分组分索引被用作一个标志,以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
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

 /* 设置相和组分索引。灰分(Ash)组分索引初始化为零,以及所有其他索引。
    灰分组分索引被用作一个标志,以确保SetSpeciesIndex只执行一次。这是由FLUENT GUI中非均相反应面板中定义的第一个反应完成的。
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

     if(IS_ASH == 0)
      SetSpeciesIndex();

     double y_o2 = yi[IP_O2][IS_O2];
     if(y_o2 > spe_small_comb)
      {
       double tmp1_exp, tmp2_exp;

       if (rp_ke)
         rr_turb = Turbulent_rr(c, t, r, yi);

       double P_o2 = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * T_g*y_o2/mw[IP_O2][IS_O2]/ 101325.;
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
     Thread *tp = pt[0];
     int index_phase = Get_Phase_Index(hr);
     Thread *ts = pt[index_phase];

     *rr = 0;
     double direction = 0.0, mol_weight, y_carbon;

     if(IS_ASH == 0)
      SetSpeciesIndex();

     double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp);
     double p_h2o = RoRT * yi[IP_H2O][IS_H2O]/mw[IP_H2O][IS_H2O]/ 101325.;
     double p_co = RoRT * yi[IP_CO][IS_CO]/mw[IP_CO][IS_CO] / 101325.;
     double p_h2 = RoRT * yi[IP_H2][IS_H2]/mw[IP_H2][IS_H2] / 101325.;

     SolidFuel_Reactant(c, t, hr, &y_carbon, &mol_weight);

     if(C_VOF(c, ts) >= eps_s_small)
      {
        *rr = rr_steam_gasif(c, t, ts, tp, p_h2o, p_co, p_h2, y_carbon, mol_weight, &direction);

        if( direction < 0.0)
         *rr = 0.0;
      }
}

// 定义碳烟水蒸气气化
DEFINE_HET_RXN_RATE(Soot_H2O_Gasif,c,t,r,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0];

     *rr = 0;
     double rr_turb = 1e+20;
     double T_g = MIN((MAX(TMIN,C_T(c,tp))),TMAX);

     if(IS_ASH == 0)
      SetSpeciesIndex();

     double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp);
     double p_h2o = RoRT * yi[IP_H2O][IS_H2O]/mw[IP_H2O][IS_H2O]/ 101325.;
     double p_h2 = RoRT * yi[IP_H2][IS_H2]/mw[IP_H2][IS_H2] / 101325.;

     if (rp_ke)
       rr_turb = Turbulent_rr(c, t, r, yi);

     double prod  =  yi[IP_SOOT][IS_SOOT]*(C_R(c,tp)*1e-03)/mw[IP_SOOT][IS_SOOT]*C_VOF(c,tp);
     *rr = A_soot_steam_gasification*exp(-E_soot_steam_gasification/Rgas/T_g)* Annealing_soot_steam_gasification * prod *
            pow(p_h2o, N_soot_steam_gasification)/(1.+K_soot_steam_gasification*p_h2);
     *rr *= 1000.;
     *rr = MIN(*rr, rr_turb);
}

// 定义MGAS逆水蒸气气化
DEFINE_HET_RXN_RATE(SteamGasif_Rev_MGAS,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0];
     int index_phase = Get_Phase_Index(hr);
     Thread *ts = pt[index_phase];

     *rr = 0;
     double direction = 0.0, mol_weight, y_carbon, rr_turb = 1e+20;

     if(IS_ASH == 0)
      SetSpeciesIndex();
     if(MGAS_Gasif)
       {
        double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp);
        double p_h2o = RoRT * yi[IP_H2O][IS_H2O]/mw[IP_H2O][IS_H2O]/ 101325.;
        double p_co = RoRT * yi[IP_CO][IS_CO]/mw[IP_CO][IS_CO] / 101325.;
        double p_h2 = RoRT * yi[IP_H2][IS_H2]/mw[IP_H2][IS_H2] / 101325.;

	y_carbon = yi[IP_SOOT][IS_SOOT];
	mol_weight = mw[IP_SOOT][IS_SOOT];

     if (rp_ke)
        rr_turb = Turbulent_rr(c, t, hr, yi);

        if(C_VOF(c, ts) >= eps_s_small)
         {
           *rr = rr_steam_gasif(c, t, ts, tp, p_h2o, p_co, p_h2, y_carbon, mol_weight, &direction);

           if( direction > 0.0)
            *rr = 0.0;
	   else
	    {
	     *rr = abs(*rr);
	     *rr = MIN(*rr, rr_turb);
	    }
         }
       }
}

// 水蒸气气化辅助函数
double rr_steam_gasif(cell_t c, Thread *t, Thread *ts, Thread *tp, double p_h2o, double p_co, double p_h2, double y_carbon, double mol_weight, double* direction)
{
     double rate, prod, T_g = MIN((MAX(TMIN,C_T(c,tp))),TMAX);
     double p_h2o_star = p_h2 * p_co / ( exp(17.29 - 16326/T_g) );

     if(MGAS_Gasif) *direction = p_h2o - p_h2o_star;
     if(PCCL_Gasif) *direction = pow(p_h2o, N_steam_gasification)/(1.+K_steam_gasification*p_h2);

     prod  =  y_carbon*(C_R(c,ts)*1e-03)/mol_weight*C_VOF(c,ts);
     if(MGAS_Gasif && *direction < 0.0)
       prod = y_carbon*(C_R(c,tp)*1e-03)/mol_weight*C_VOF(c,tp);

     rate = A_steam_gasification*exp(-E_steam_gasification/Rgas/T_g)* Annealing_steam_gasification * prod * *direction;
     rate *= 1000.;
     return rate;
}

// 定义CO2气化
DEFINE_HET_RXN_RATE(Co2Gasif,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0];
     int index_phase = Get_Phase_Index(hr);
     Thread *ts = pt[index_phase];

     *rr = 0;
     double direction = 0.0, mol_weight, y_carbon;

     if(IS_ASH == 0)
       SetSpeciesIndex();

     double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp);
     double p_co = RoRT * yi[IP_CO][IS_CO]/mw[IP_CO][IS_CO] / 101325.;
     double p_co2 = RoRT * yi[IP_CO2][IS_CO2]/mw[IP_CO2][IS_CO2] / 101325.;

     SolidFuel_Reactant(c, t, hr, &y_carbon, &mol_weight);

     if(C_VOF(c, ts) >= eps_s_small)
       {
        *rr = rr_co2_gasif(c, t, ts, tp, p_co, p_co2, y_carbon, mol_weight, &direction);

        if( direction < 0.0)
         *rr = 0.0;
       }
}

// 碳烟CO2气化
DEFINE_HET_RXN_RATE(Soot_CO2_Gasif,c,t,r,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0];

     *rr = 0;
     double rr_turb = 1e+20;
     double T_g = MIN((MAX(TMIN,C_T(c,tp))),TMAX);

     if(IS_ASH == 0)
      SetSpeciesIndex();

     double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp);
     double p_co = RoRT * yi[IP_CO][IS_CO]/mw[IP_CO][IS_CO]/ 101325.;
     double p_co2 = RoRT * yi[IP_CO2][IS_CO2]/mw[IP_CO2][IS_CO2] / 101325.;

     if (rp_ke)
       rr_turb = Turbulent_rr(c, t, r, yi);

     double prod  =  yi[IP_SOOT][IS_SOOT]*(C_R(c,tp)*1e-03)/mw[IP_SOOT][IS_SOOT]*C_VOF(c,tp);
     *rr = A_soot_co2_gasification*exp(-E_soot_co2_gasification/Rgas/T_g)* Annealing_soot_co2_gasification * prod *
            pow(p_co2, N_soot_co2_gasification)/(1.+K_soot_co2_gasification*p_co);
     *rr *= 1000.;
     *rr = MIN(*rr, rr_turb);
}

// MGAS逆CO2气化
DEFINE_HET_RXN_RATE(Co2Gasif_Rev_MGAS,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0];
     int index_phase = Get_Phase_Index(hr);
     Thread *ts = pt[index_phase];

     *rr = 0;
     double direction = 0.0, mol_weight, y_carbon, rr_turb = 1e+20;

     if(IS_ASH == 0)
       SetSpeciesIndex();

     if(MGAS_Gasif)
       {
         double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp);
         double p_co = RoRT * yi[IP_CO][IS_CO]/mw[IP_CO][IS_CO] / 101325.;
         double p_co2 = RoRT * yi[IP_CO2][IS_CO2]/mw[IP_CO2][IS_CO2] / 101325.;

         y_carbon = yi[IP_SOOT][IS_SOOT];
         mol_weight = mw[IP_SOOT][IS_SOOT];

     if (rp_ke)
         rr_turb = Turbulent_rr(c, t, hr, yi);

         if(C_VOF(c, ts) >= eps_s_small)
           {
            *rr = rr_co2_gasif(c, t, ts, tp, p_co, p_co2, y_carbon, mol_weight, &direction);

            if( direction > 0.0)
             *rr = 0.0;
	    else
	    {
	     *rr = abs(*rr);
	     *rr = MIN(*rr, rr_turb);
	    }
           }
       }
}

// CO2气化辅助函数
double rr_co2_gasif(cell_t c, Thread *t, Thread *ts, Thread *tp, double p_co, double p_co2, double y_carbon, double mol_weight, double* direction)
{
     double T_g = MIN(MAX(TMIN,C_T(c,tp)), TMAX), prod;
     double p_co2_star = p_co * p_co/(exp(20.92 - 20282/T_g));

     if(MGAS_Gasif)  *direction = p_co2-p_co2_star;
     if(PCCL_Gasif)  *direction = pow(p_co2, N_co2_gasification)/(1. + K_co2_gasification * p_co);

     prod  =  y_carbon*C_R(c,ts)*1.e-3/mol_weight* C_VOF(c,ts);
     if(MGAS_Gasif && *direction < 0.0)
       prod = y_carbon*(C_R(c,tp)*1e-03)/mol_weight*C_VOF(c,tp);

     double rate = A_co2_gasification*exp(-E_co2_gasification/Rgas/T_g)*Annealing_co2_gasification * prod * (*direction);
     rate *= 1000.;
     return rate;
}

// H2气化
DEFINE_HET_RXN_RATE(H2Gasif,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0];
     int index_phase = Get_Phase_Index(hr);
     Thread *ts = pt[index_phase];

     *rr = 0;
     double direction = 0.0, mol_weight, y_carbon;

     if(IS_ASH == 0)
       SetSpeciesIndex();

     double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp);
     double p_h2 = RoRT * yi[IP_H2][IS_H2]/mw[IP_H2][IS_H2] / 101325.;
     double p_ch4 = RoRT * yi[IP_CH4][IS_CH4]/mw[IP_CH4][IS_CH4] / 101325.;

     SolidFuel_Reactant(c, t, hr, &y_carbon, &mol_weight);

     if(C_VOF(c, ts) >= eps_s_small)
       {
        *rr = rr_h2_gasif(c, t, ts, tp, p_h2, p_ch4, y_carbon, mol_weight, &direction);

        if( direction < 0.0)
         *rr = 0.0;
       }
}

// 碳烟H2气化
DEFINE_HET_RXN_RATE(Soot_H2_Gasif,c,t,r,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0];

     *rr = 0;
     double rr_turb = 1e+20;
     double T_g = MIN((MAX(TMIN,C_T(c,tp))),TMAX);

     if(IS_ASH == 0)
      SetSpeciesIndex();

     double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp);
     double p_h2 = RoRT * yi[IP_H2][IS_H2]/mw[IP_H2][IS_H2] / 101325.;

     if (rp_ke)
       rr_turb = Turbulent_rr(c, t, r, yi);

     double prod  =  yi[IP_SOOT][IS_SOOT]*(C_R(c,tp)*1e-03)/mw[IP_SOOT][IS_SOOT]*C_VOF(c,tp);
     *rr = A_soot_h2_gasification*exp(-E_soot_h2_gasification/Rgas/T_g)* Annealing_soot_h2_gasification * prod *
            pow(p_h2, N_soot_h2_gasification);
     *rr *= 1000.;
     *rr = MIN(*rr, rr_turb);
}

// MGAS逆H2气化
DEFINE_HET_RXN_RATE(H2Gasif_Rev_MGAS,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0];
     int index_phase = Get_Phase_Index(hr);
     Thread *ts = pt[index_phase];

     *rr = 0;
     double direction = 0.0, mol_weight, y_carbon, rr_turb = 1e+20;

     if(IS_ASH == 0)
       SetSpeciesIndex();
     if(MGAS_Gasif)
       {
        double RoRT = C_R(c,tp) * UNIVERSAL_GAS_CONSTANT * C_T(c,tp);
        double p_h2 = RoRT * yi[IP_H2][IS_H2]/mw[IP_H2][IS_H2] / 101325.;
        double p_ch4 = RoRT * yi[IP_CH4][IS_CH4]/mw[IP_CH4][IS_CH4] / 101325.;

        y_carbon = yi[IP_SOOT][IS_SOOT];
        mol_weight = mw[IP_SOOT][IS_SOOT];

        if(C_VOF(c, ts) >= eps_s_small)
          {
            if (rp_ke)
            rr_turb = Turbulent_rr(c, t, hr, yi);

           *rr = rr_h2_gasif(c, t, ts, tp, p_h2, p_ch4, y_carbon, mol_weight, &direction);

           if( direction > 0.0)
            *rr = 0.0;
	   else
	    {
	     *rr = abs(*rr);
	     *rr = MIN(*rr, rr_turb);
	    }
          }
       }
}

// H2气化辅助函数
double rr_h2_gasif(cell_t c, Thread *t, Thread *ts, Thread *tp, double p_h2, double p_ch4, double y_carbon, double mol_weight, double* direction)
{
     double rate = 0.0, prod;
     double T_g = MIN((MAX(TMIN,C_T(c,tp))), TMAX);
     double p_h2_star = pow ((p_ch4/(exp(-13.43 + 10999/T_g))), 0.5);

     prod  =  y_carbon*C_R(c,ts)*1.e-3/mol_weight * C_VOF(c,ts);

     if(MGAS_Gasif)
       {
         *direction = p_h2-p_h2_star;
         if(*direction < 0.0)
           prod = y_carbon*(C_R(c,tp)*1e-03)/mol_weight*C_VOF(c,tp);
         rate = exp( -7.087 - 8078/T_g )* prod * *direction ;
       }
     if(PCCL_Gasif)
       {
        *direction = p_h2;
        rate = A_h2_gasification*exp(-E_h2_gasification/Rgas/T_g)*Annealing_h2_gasification * prod * *direction;
       }

     rate *= 1000.;
     return rate;
}

// 煤燃烧
DEFINE_HET_RXN_RATE(coal_combustion,c,t,hr,mw,yi,rr,rr_t)
{
     Thread **pt = THREAD_SUB_THREADS(t);
     Thread *tp = pt[0];
     int index_phase = Get_Phase_Index(hr);
     Thread *ts = pt[index_phase];
     double mol_weight, y_carbon, y_ash;
     *rr = 0.0;

     if(IS_ASH == 0)
       SetSpeciesIndex();

     if( C_YI(c,tp,IS_O2) >= spe_small)
      {
        SolidFuel_Reactant(c, t, hr, &y_carbon, &mol_weight);
	y_ash = yi[index_phase][IS_ASH];

        *rr = rr_combustion(c, t, ts, tp, yi[IP_O2][IS_O2], y_ash, y_carbon);
        *rr *=  1000.;
      }
}

// 燃烧速率计算
double rr_combustion(cell_t c, Thread *t, Thread *ts, Thread *tp, double yi_O2, double y_ash, double y_carbon)
{
     double rd, k_f, k_r, factor, k_a, rate = 0.0, vrel;
     double Pt = MAX(0.1, (op_pres+C_P(c,t))/101325);
     double gas_constant = 82.06;
     double T = C_T(c,tp), T_s = C_T(c,ts), D_p = C_PHASE_DIAMETER(c,ts)*100.;
     double p_o2 = C_R(c,tp)*UNIVERSAL_GAS_CONSTANT* T *yi_O2/mw[IP_O2][IS_O2] / 101325.;

     if(fc_ar > 0.)
       {
        if (y_ash > 0.)
          {
            rd = pow( (y_carbon * ash_ar/100.)/(y_ash * fc_ar/100.), (1./3.) );
            rd = MIN(1., rd);
          }
        else
            rd = 1.;
       }
     else
        rd = 0.;

     double diff = MAX((4.26 * pow((T/1800.),1.75)/Pt), 1.e-10);
     double Sc1o3 = pow(C_MU_L(c,tp)/(C_R(c,tp) * diff * 1.e-4), 1./3.);
#if RP_2D
     vrel = pow(( (C_U(c,tp)-C_U(c,ts))*(C_U(c,tp)-C_U(c,ts)) +
            (C_V(c,tp)-C_V(c,ts))*(C_V(c,tp)-C_V(c,ts))), 0.5);
#endif
#if RP_3D
     vrel = pow(( (C_U(c,tp)-C_U(c,ts))*(C_U(c,tp)-C_U(c,ts)) +
            (C_V(c,tp)-C_V(c,ts))*(C_V(c,tp)-C_V(c,ts)) +
            (C_W(c,tp)-C_W(c,ts))*(C_W(c,tp)-C_W(c,ts)) ), 0.5);
#endif
     double Re = C_VOF(c,tp) * D_p/100. * vrel * C_R(c,tp)/(C_MU_L(c,tp)+SMALL_S);
     double N_sherwood = (7. - 10. * C_VOF(c,tp) + 5. * C_VOF(c,tp) * C_VOF(c,tp) )*
                         (1. + 0.7 * pow(Re, 0.2) * Sc1o3) +
                         (1.33 - 2.4 * C_VOF(c,tp) + 1.2 * C_VOF(c,tp) * C_VOF(c,tp)) *
                         pow(Re, 0.7) * Sc1o3;

     if ( rd <= 0. || C_VOF(c, ts) <= 0. )
       {
          rate = 0.;
       }
     else
       {
          k_f = diff * N_sherwood / (D_p * gas_constant/mw[IP_O2][IS_O2] * T );
          k_r = A_c_combustion * exp( -E_c_combustion/Rgas/T_s ) * rd * rd;
          if ( rd >= 1.)
           {
            rate = 1. / (1./k_f + 1./k_r);
           }
          else
           {
            k_a = 2. * rd * diff * f_ep_a / (D_p * (1.-rd) * gas_constant/mw[IP_O2][IS_O2] * T_s );
            rate = 1. / (1./k_f + 1./k_r + 1./k_a);
           }

          factor = y_carbon / (y_carbon + 1.e-6);
          rate *= p_o2 * 6. * C_VOF(c,ts) * factor / (D_p * 32.);
       }
      return rate;
}

#if !RP_NODE || !PARALLEL

// 读取挥发分质量分数
void volatile_mass_fractions()
{
       read_c3m_data();

       CX_Message("PCCL_Devol = %d\n",PCCL_Devol);
       CX_Message("MGAS_Devol = %d\n",MGAS_Devol);
       CX_Message("CPD_Devol = %d\n",CPD_Devol);
       CX_Message("FGDVC_Devol = %d\n",FGDVC_Devol);
       CX_Message("HPTR_Devol = %d\n",HPTR_Devol);
       CX_Message("MGAS_Moisture = %d\n",MGAS_Moisture);
       CX_Message("PCCL_Moisture = %d\n",PCCL_Moisture);
       CX_Message("MGAS_TarCracking = %d\n",MGAS_TarCracking);
       CX_Message("PCCL_2nd_Pyro = %d\n",PCCL_2nd_Pyro);
       CX_Message("MGAS_Gasif = %d\n",MGAS_Gasif);
       CX_Message("PCCL_Gasif = %d\n",PCCL_Gasif);
       CX_Message("PCCL_TarCracking = %d\n",PCCL_TarCracking);
       CX_Message("MGAS_WGS = %d\n",MGAS_WGS);
       CX_Message("PCCL_soot_gasif = %d\n",PCCL_soot_gasif);
       CX_Message("MGAS_char_combustion = %d\n",MGAS_char_combustion);
       CX_Message("PCCL_char_combustion = %d\n",PCCL_char_combustion);
       CX_Message("PCCL_soot_oxidation = %d\n",PCCL_soot_oxidation);
       CX_Message("TAR_oxidation = %d\n",TAR_oxidation);
       CX_Message("MGAS_gas_phase_oxidation = %d\n",MGAS_gas_phase_oxidation);

       CX_Message("fc_ar = %f\n",fc_ar);
       CX_Message("vm_ar = %f\n",vm_ar);
       CX_Message("ash_ar = %f\n",ash_ar);
       CX_Message("moist_ar = %f\n",moist_ar);
       CX_Message("f_ep_a = %f\n",f_ep_a);

       CX_Message("a_moisture_release = %f\n",A_moisture_release);
       CX_Message("e_moisture_release = %f\n",E_moisture_release);
       CX_Message("a1_devolatilization = %f\n",A1_devolatilization);
       CX_Message("e1_devolatilization = %f\n",E1_devolatilization);
       CX_Message("a2_devolatilization = %f\n",A2_devolatilization);
       CX_Message("e2_devolatilization = %f\n",E2_devolatilization);
       CX_Message("a_tar_cracking = %f\n",A_tar_cracking);
       CX_Message("e_tar_cracking = %f\n",E_tar_cracking);
       CX_Message("a_steam_gasification = %f\n",A_steam_gasification);
       CX_Message("e_steam_gasification = %f\n",E_steam_gasification);
       CX_Message("a_co2_gasification = %f\n",A_co2_gasification);
       CX_Message("e_co2_gasification = %f\n",E_co2_gasification);
       CX_Message("a_h2_gasification = %f\n",A_h2_gasification);
       CX_Message("e_h2_gasification = %f\n",E_h2_gasification);
       CX_Message("a_c_combustion = %f\n",A_c_combustion);
       CX_Message("e_c_combustion = %f\n",E_c_combustion);

       double ep_a = 0.25 + 0.75*(1-ash_ar/100.);
       f_ep_a = pow(ep_a,2.5);
}

// 设置布尔值
static void SetBooleanValue(char * var , char * svalue)
{
	cxboolean value;

	if ( strcmp(svalue,"true") == 0)
		value = TRUE;
	else
		value = FALSE;

	if (strcmp(var,"pccl_devol") == 0) PCCL_Devol = value;
	if (strcmp(var,"mgas_devol") == 0) MGAS_Devol = value;
	if (strcmp(var,"cpd_devol") == 0) CPD_Devol = value;
	if (strcmp(var,"fgdvc_devol") == 0) FGDVC_Devol = value;
	if (strcmp(var,"hptr_devol") == 0) HPTR_Devol = value;
	if (strcmp(var,"pccl_moisture") == 0) PCCL_Moisture = value;
	if (strcmp(var,"mgas_moisture") == 0) MGAS_Moisture = value;
	if (strcmp(var,"pccl_2nd_pyro") == 0) PCCL_2nd_Pyro = value;
	if (strcmp(var,"mgas_tarcracking") == 0) MGAS_TarCracking = value;
	if (strcmp(var,"pccl_tarcracking") == 0) PCCL_TarCracking = value;
	if (strcmp(var,"mgas_gasif") == 0) MGAS_Gasif = value;
	if (strcmp(var,"pccl_gasif") == 0) PCCL_Gasif = value;
	if (strcmp(var,"pccl_soot_gasif") == 0) PCCL_soot_gasif = value;
	if (strcmp(var,"mgas_char_combustion") == 0) MGAS_char_combustion = value;
	if (strcmp(var,"pccl_char_combustion") == 0) PCCL_char_combustion = value;
	if (strcmp(var,"pccl_soot_oxidation") == 0) PCCL_soot_oxidation = value;
	if (strcmp(var,"tar_oxidation") == 0) TAR_oxidation = value;
	if (strcmp(var,"mgas_gas_phase_oxidation") == 0) MGAS_gas_phase_oxidation = value;
	if (strcmp(var,"mgas_wgs") == 0) MGAS_WGS = value;
}

// 设置数值
static void SetValue(char * var , char * svalue)
{
	char * pEnd;
	double value = strtod(svalue,&pEnd);

	// 煤质分析
	if (strcmp(var,"fc_ar") == 0) fc_ar = value;
	if (strcmp(var,"vm_ar") == 0) vm_ar = value;
	if (strcmp(var,"ash_ar") == 0) ash_ar = value;
	if (strcmp(var,"moist_ar") == 0) moist_ar = value;
	if (strcmp(var,"f_ep_a") == 0) f_ep_a = value;

	// 水分释放
	if (strcmp(var,"a_moisture_release") == 0) A_moisture_release = value;
	if (strcmp(var,"e_moisture_release") == 0) E_moisture_release = value;
	if (strcmp(var,"moisture_flux") == 0) Moisture_Flux = value;
	if (strcmp(var,"wg3") == 0) wg3 = value;

	// 脱挥发分
	if (strcmp(var,"a1_devolatilization") == 0) A1_devolatilization = value;
	if (strcmp(var,"e1_devolatilization") == 0) E1_devolatilization = value;
	if (strcmp(var,"a2_devolatilization") == 0) A2_devolatilization = value;
	if (strcmp(var,"e2_devolatilization") == 0) E2_devolatilization = value;

	// 焦油裂解
	if (strcmp(var,"a_tar_cracking") == 0) A_tar_cracking = value;
	if (strcmp(var,"e_tar_cracking") == 0) E_tar_cracking = value;

	// 水蒸气气化
	if (strcmp(var,"a_steam_gasification") == 0) A_steam_gasification = value;
	if (strcmp(var,"e_steam_gasification") == 0) E_steam_gasification = value;
	if (strcmp(var,"k_steam_gasification") == 0) K_steam_gasification = value;
	if (strcmp(var,"n_steam_gasification") == 0) N_steam_gasification = value;
	if (strcmp(var,"annealing_steam_gasification") == 0) Annealing_steam_gasification = value;

	// CO2气化
	if (strcmp(var,"a_co2_gasification") == 0) A_co2_gasification = value;
	if (strcmp(var,"e_co2_gasification") == 0) E_co2_gasification = value;
	if (strcmp(var,"k_co2_gasification") == 0) K_co2_gasification = value;
	if (strcmp(var,"n_co2_gasification") == 0) N_co2_gasification = value;
	if (strcmp(var,"annealing_co2_gasification") == 0) Annealing_co2_gasification = value;

	// H2气化
	if (strcmp(var,"a_h2_gasification") == 0) A_h2_gasification = value;
	if (strcmp(var,"e_h2_gasification") == 0) E_h2_gasification = value;
	if (strcmp(var,"n_h2_gasification") == 0) N_h2_gasification = value;
	if (strcmp(var,"annealing_h2_gasification") == 0) Annealing_h2_gasification = value;

	// 碳烟-水蒸气气化
	if (strcmp(var,"a_soot_steam_gasification") == 0) A_soot_steam_gasification = value;
	if (strcmp(var,"e_soot_steam_gasification") == 0) E_soot_steam_gasification = value;
	if (strcmp(var,"k_soot_steam_gasification") == 0) K_soot_steam_gasification = value;
	if (strcmp(var,"n_soot_steam_gasification") == 0) N_soot_steam_gasification = value;
	if (strcmp(var,"annealing_soot_steam_gasification") == 0) Annealing_soot_steam_gasification = value;

	// 碳烟-CO2气化
	if (strcmp(var,"a_soot_co2_gasification") == 0) A_soot_co2_gasification = value;
	if (strcmp(var,"e_soot_co2_gasification") == 0) E_soot_co2_gasification = value;
	if (strcmp(var,"k_soot_co2_gasification") == 0) K_soot_co2_gasification = value;
	if (strcmp(var,"n_soot_co2_gasification") == 0) N_soot_co2_gasification = value;
	if (strcmp(var,"annealing_soot_co2_gasification") == 0) Annealing_soot_co2_gasification = value;

	// 碳烟-H2气化
	if (strcmp(var,"a_soot_h2_gasification") == 0) A_soot_h2_gasification = value;
	if (strcmp(var,"e_soot_h2_gasification") == 0) E_soot_h2_gasification = value;
	if (strcmp(var,"n_soot_h2_gasification") == 0) N_soot_h2_gasification = value;
	if (strcmp(var,"annealing_soot_h2_gasification") == 0) Annealing_soot_h2_gasification = value;

	// 碳烟燃烧
	if (strcmp(var,"a_soot_combustion") == 0) A_Soot_Combustion = value;
	if (strcmp(var,"e_soot_combustion") == 0) E_Soot_Combustion = value;

	// 焦炭燃烧
	if (strcmp(var,"a_c_combustion") == 0) A_c_combustion = value;
	if (strcmp(var,"e_c_combustion") == 0) E_c_combustion = value;
	if (strcmp(var,"n_c_combustion") == 0) N_c_combustion = value;
	if (strcmp(var,"annealing_c_combustion") == 0) Annealing_c_combustion = value;
}

// 读取c3m数据
void read_c3m_data()
{
    CX_Message("start of read_c3m_data \n");

	FILE * pFile;
	char line[80];
	char field1[80];
	char field2[80];
	char *pch;
	int  i , field_index;

	pFile = fopen("fluent_c3m_udf.inp","r");

	if (pFile != NULL)
	{
	    CX_Message("fopen OK \n");
		char * pEnd = line;

		while (pEnd != NULL)
		{
			pEnd = fgets(line,80,pFile);

			if (pEnd != NULL)
			{
				for (i=0; i<80; ++i) line[i] = tolower(line[i]);

				pch = strtok(line," ,\t\n=");

				field_index = 0;
				while (pch != NULL)
				{
					if (field_index == 0)
					{
						strcpy(field1,pch);
						field_index = 1;
					}
					else
					{
						strcpy(field2,pch);
						field_index = 0;
						SetValue(field1,field2);
						SetBooleanValue(field1,field2);
					}
					pch = strtok(NULL," ,\t\n=");
				}
			}
		}

		fclose(pFile);
		CX_Message("read_c3m_data completed successfully\n");
	}
	else
	{
		CX_Message("ERROR: Could not open fluent_c3m_udf.inp file!\n");
		CX_Message("Using default parameter values.\n");
	}
}
#endif

// 交换属性 - 煤传热
DEFINE_EXCHANGE_PROPERTY(Heat_Trans_Coal, c, t, i, j)
{
  Thread *ti = THREAD_SUB_THREAD(t,i);
  Thread *tj = THREAD_SUB_THREAD(t,j);
  double val;

  val = heat_gunn_udf(c,ti, tj);
  return val;
}

// 交换属性 - 再循环传热
DEFINE_EXCHANGE_PROPERTY(Heat_Trans_Recy, c, t, i, j)
{
  Thread *ti = THREAD_SUB_THREAD(t,i);
  Thread *tj = THREAD_SUB_THREAD(t,j);
  double val;

  val = heat_gunn_udf(c,ti, tj);
  return val;
}

// Gunn传热模型
double heat_gunn_udf(cell_t c, Thread *ti, Thread *tj)
{
  double h;
  double d = C_PHASE_DIAMETER(c,tj);
  double k = C_K_L(c,ti);
  double vf = C_VOF(c,ti);
  double vf2 = vf*vf;
  double vel, Re, Pr, Nu;

#if RP_2D
    vel = pow(((C_U(c,tj)-C_U(c,ti))*(C_U(c,tj)-C_U(c,ti)) + (C_V(c,tj)-C_V(c,ti))*(C_V(c,tj)-C_V(c,ti))),0.5);
#endif

#if RP_3D
    vel = pow(((C_U(c,tj)-C_U(c,ti))*(C_U(c,tj)-C_U(c,ti)) + (C_V(c,tj)-C_V(c,ti))*(C_V(c,tj)-C_V(c,ti)) +
          (C_W(c,tj)-C_W(c,ti))*(C_W(c,tj)-C_W(c,ti))),0.5);
#endif

  Re = RE_NUMBER(C_R(c,ti),vel,d,C_MU_L(c,ti));
  Pr = PR_NUMBER (C_CP(c,ti),C_MU_L(c,ti),k);
  Pr = pow (Pr,1./3.);
  Nu = (7. - 10*vf + 5.*vf2)*(1. + 0.7*pow(Re,0.2)*Pr) +
    (1.33 - 2.4*vf + 1.2*vf2)*pow(Re,0.7)*Pr;

  h =IP_HEAT_COEFF(C_VOF(c,ti),C_VOF(c,tj),k,Nu,d);
  return h;
}

// 饱和压力
double satPressure(double T)
{
  const double Tstd = 273.15;
  double TT = T - Tstd;
  double p_sat;

    p_sat = (100.0 * 6.1121 * exp((18.678-TT/234.5)*(TT/(257.14+TT))));
  return p_sat;
}

// 获取相索引
double Get_Phase_Index(Hetero_Reaction *hr)
{
     Domain **dr = hr->domain_reactant;
     int i, index_phase = 0, iphase;
     for (i = 0; i < hr->n_reactants; i++)
     {
       iphase = DOMAIN_INDEX(dr[i]);
       if(iphase > 0)
         index_phase = iphase;
     }
     return index_phase;
}

// 湍流反应速率
double Turbulent_rr(cell_t c, Thread *t, Hetero_Reaction *r, real yi[MAX_PHASES][MAX_SPE_EQNS])
{
  Thread **pt = THREAD_SUB_THREADS(t);
  Thread *tp = pt[0];
  double minR = 1.e+20, numerator = 0.;
  double denom =0.,ci;

  int i;
  double Amix = 4., Bmix = 0.5;

  for(i=0; i< r->n_reactants; ++i)
    {
      if(r->stoich_reactant[i] != 0.0)
        {
          ci = C_R(c,tp)*yi[0][r->reactant[i]]/mw[0][r->reactant[i]];
          SETMIN(minR, ci/r->stoich_reactant[i]);
	}
     }


  for(i=0; i< r->n_products; ++i)
    {
      if(r->stoich_product[i] != 0.0)
        {
          ci = C_R(c,tp)*yi[0][r->product[i]]/mw[0][r->product[i]];
          numerator +=ci*mw[0][r->product[i]];
          denom +=r->stoich_product[i]*mw[0][r->product[i]];
        }
    }

  double rr_turb = Amix * C_D(c,tp)/C_K(c,tp) * MIN(minR, Bmix*numerator/denom);
  return rr_turb;
}

// 固体燃料反应物
void SolidFuel_Reactant(cell_t c, Thread *t, Hetero_Reaction *hr, double* y_carbon, double* mol_weight)
{
  Thread **pt = THREAD_SUB_THREADS(t);
  int i;

  for(i=0; i< hr->n_reactants; ++i)
    {
      if(hr->stoich_reactant[i] != 0.0)
        {
	  int phase = DOMAIN_INDEX(hr->domain_reactant[i]);
          Thread *thread = pt[phase];
	  if(phase != 0)
	   {
	     *y_carbon = C_YI(c,thread,hr->reactant[i]);
	     *mol_weight = mw[phase][hr->reactant[i]];
	   }
        }
    }
}