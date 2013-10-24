// header_main.h

// main プログラムのみで
// include する(multiple definition を防ぐため)

#ifndef _HEADER_MAIN_H_DEFINED_

#include <stdio.h>
#include "atom_struct.h" 



//@@@@@@@@@   Clラジカル/イオン の flux比   @@@@@@@@@@@
// (argument によって変更されうるので、const にしない)

double NEUTRAL_ION_RATIO       =  100.0;//50.0;//10.0;//10;//1;//20;//5;// 

double INCIDENT_ENERGY_ION     =  100.0;//200.0;//50.0;//55.0;//150.0;//75.0;//(eV)

const double INCIDENT_ENERGY_NEUTRAL =   2.0 ; 


// yield の係数β = C(√E - √Eth): C = 0.90537, Eth = 21.719987
// （入射エネルギーの関数）
// modified on 14 Aug 2004   : C = 0.1 , Eth = 10.0
// 
// Cの値をどう決定するかについて、2004/Sep に再検討。
// Chang et al.(1998) のFig. 3 によると、incident energy: 35eV
// Γn/Γi = 200 , normal incidence で yield はほぼ 1.0
// Eth = 10.0 として C を求めると、C = 0.363 となる。


double  YIELD_BETA    
//=         0.36 * (sqrt( INCIDENT_ENERGY_ION ) -  sqrt(10.0 )) ;
//=         0.90537 * (sqrt( INCIDENT_ENERGY_ION ) -  sqrt(10.0 )) ;
=         0.77 * (sqrt( INCIDENT_ENERGY_ION ) -  sqrt(20.0 )) ;

//==================================================
//== 前方散乱の導入（隣接セルでの衝突判定）の flag 
const bool   FLAG_FORWARD_SCATTER = true;//false;//

//==================================================
//== マスクの浸食の flag 
const bool   FLAG_MASK_EROSION = true;//false;//

//==================================================
//== 片側open space条件の flag 
const bool   FLAG_OPEN_SPACE   = false;//true;//

//==================================================
//== Chemical etching の Yield の計算    ===========
//   Sano's master's thesis に基づく。
//   （反応 SiCl + Cl(g) -> SiCl2 (g) の反応確率）
//   係数は Ogryzlo, JAP(1990)より
//   実験では基板温度〜ガス温度が等しいとみなしてよい？
//   ここでは、flux 計算も基板温度に基づいて行う。

//==== Chemical etching を導入するかどうかのflag ===
const bool   FLAG_CHEMICAL_ETCH  = false;//true;// 

const double N_DOPANT    = 1.0e+20; //== ドーパント密度(cm-3)
const double T_SUBSTRATE = 0.1;//300.0;//350.0;//== 基板温度 (K)
const double T_GAS       = 300.0;//== 気相の温度 (K)

const double YIELD_CHEMICAL  //= 0.01;//0.005 ; 
= (4.04e-18 * pow(N_DOPANT, 0.39) * sqrt(T_SUBSTRATE) 
   *  exp( - 4.70 / (K_BOLTZMANN_kcal * T_SUBSTRATE) )) //<- Ogryzloの式
  *   (1.0e-8 / 60.0) * Si.density_atom * 1.0e+22   //<- 原子数に換算
/  ( 0.25 * 1.0e+2 *                                //<- フラックス 
     sqrt( 8.0 * K_BOLTZMANN * T_GAS / (PI * CL_MASS ))) ; 

// ==== flux のカウントを行うか否かの flag（TRUE: カウントする）====
const bool FLAG_FLUX_COUNT = false;//true;//

//************************************
//**  standard I/O
//************************************
// --- 標準出力する頻度
const int N_COUT_INTERVAL   = 100;// 1;//

//************************************
//**  file I/O
//************************************

// === 形状ファイルを読み込むか（TRUE: 読み込む）=====
const bool FLAG_INPUT_PROFILE = false;//true;//
char CL_INPUT_FILE[50] ; 
int input_sprintf1 = sprintf(CL_INPUT_FILE , "Cl_bond543950.dat") ;
char SI_INPUT_FILE[50] ; 
int input_sprintf2 = sprintf(SI_INPUT_FILE , "Si_matrix543950.dat") ;

char CL_BOND_FILE[50] ; 
int tmp_sprintf1 = sprintf(CL_BOND_FILE   , "Cl_bond") ;
char SI_MATRIX_FILE[50] ; 
int tmp_sprintf2 = sprintf(SI_MATRIX_FILE , "Si_matrix") ;
char CONDITION_FILE[50] ; 
int tmp_sprintf3 = sprintf(CONDITION_FILE , "condition.txt") ;
char YIELD_FILE[50] ; 
int tmp_sprintf4 = sprintf(YIELD_FILE ,     "yield") ;

char ION_COUNTER_FILE[50] ;
//int tmp_sprintf5 = sprintf(ION_COUNTER_FILE ,"ion_count") ;
int tmp_sprintf5 = sprintf(ION_COUNTER_FILE ,"ion_count20R") ;
char NEUTRAL_COUNTER_FILE[50] ;
//int tmp_sprintf6 = sprintf(NEUTRAL_COUNTER_FILE ,"neutral_count") ;
int tmp_sprintf6 = sprintf(NEUTRAL_COUNTER_FILE ,"neutral_Sn1_") ;

/*
  #define CL_BOND_FILE    "Cl_bond"
  #define SI_MATRIX_FILE  "Si_matrix"
  
  #define CONDITION_FILE  "condition.dat"
  #define YIELD_FILE      "yield.dat"
*/

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//#define IADF_FILE       "iadfR50.dat"
//#define IADF_FILE       "iadf50.dat"
#define IADF_FILE       "iadf100.dat"
//#define IADF_FILE       "iadf150.dat"
//#define IADF_FILE       "iadf200.dat"
#define CL_DENSITY_DEPTH_FILE  "cl_density_depth" 


#define _HEADER_MAIN_H_DEFINED_
#endif

