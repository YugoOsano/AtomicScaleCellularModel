// header_main.h

// main �ץ����Τߤ�
// include ����(multiple definition ���ɤ�����)

#ifndef _HEADER_MAIN_H_DEFINED_

#include <stdio.h>
#include "atom_struct.h" 



//@@@@@@@@@   Cl�饸����/������ �� flux��   @@@@@@@@@@@
// (argument �ˤ�ä��ѹ����줦��Τǡ�const �ˤ��ʤ�)

double NEUTRAL_ION_RATIO       =  100.0;//50.0;//10.0;//10;//1;//20;//5;// 

double INCIDENT_ENERGY_ION     =  100.0;//200.0;//50.0;//55.0;//150.0;//75.0;//(eV)

const double INCIDENT_ENERGY_NEUTRAL =   2.0 ; 


// yield �η����� = C(��E - ��Eth): C = 0.90537, Eth = 21.719987
// �����ͥ��ͥ륮���δؿ���
// modified on 14 Aug 2004   : C = 0.1 , Eth = 10.0
// 
// C���ͤ�ɤ����ꤹ�뤫�ˤĤ��ơ�2004/Sep �˺Ƹ�Ƥ��
// Chang et al.(1998) ��Fig. 3 �ˤ��ȡ�incident energy: 35eV
// ��n/��i = 200 , normal incidence �� yield �Ϥۤ� 1.0
// Eth = 10.0 �Ȥ��� C �����ȡ�C = 0.363 �Ȥʤ롣


double  YIELD_BETA    
//=         0.36 * (sqrt( INCIDENT_ENERGY_ION ) -  sqrt(10.0 )) ;
//=         0.90537 * (sqrt( INCIDENT_ENERGY_ION ) -  sqrt(10.0 )) ;
=         0.77 * (sqrt( INCIDENT_ENERGY_ION ) -  sqrt(20.0 )) ;

//==================================================
//== ���������Ƴ�������ܥ���Ǥξ���Ƚ��ˤ� flag 
const bool   FLAG_FORWARD_SCATTER = true;//false;//

//==================================================
//== �ޥ����ο����� flag 
const bool   FLAG_MASK_EROSION = true;//false;//

//==================================================
//== ��¦open space���� flag 
const bool   FLAG_OPEN_SPACE   = false;//true;//

//==================================================
//== Chemical etching �� Yield �η׻�    ===========
//   Sano's master's thesis �˴�Ť���
//   ��ȿ�� SiCl + Cl(g) -> SiCl2 (g) ��ȿ����Ψ��
//   ������ Ogryzlo, JAP(1990)���
//   �¸��Ǥϴ��Ĳ��١��������٤��������Ȥߤʤ��Ƥ褤��
//   �����Ǥϡ�flux �׻�����Ĳ��٤˴�Ť��ƹԤ���

//==== Chemical etching ��Ƴ�����뤫�ɤ�����flag ===
const bool   FLAG_CHEMICAL_ETCH  = false;//true;// 

const double N_DOPANT    = 1.0e+20; //== �ɡ��ѥ��̩��(cm-3)
const double T_SUBSTRATE = 0.1;//300.0;//350.0;//== ���Ĳ��� (K)
const double T_GAS       = 300.0;//== ����β��� (K)

const double YIELD_CHEMICAL  //= 0.01;//0.005 ; 
= (4.04e-18 * pow(N_DOPANT, 0.39) * sqrt(T_SUBSTRATE) 
   *  exp( - 4.70 / (K_BOLTZMANN_kcal * T_SUBSTRATE) )) //<- Ogryzlo�μ�
  *   (1.0e-8 / 60.0) * Si.density_atom * 1.0e+22   //<- ���ҿ��˴���
/  ( 0.25 * 1.0e+2 *                                //<- �ե�å��� 
     sqrt( 8.0 * K_BOLTZMANN * T_GAS / (PI * CL_MASS ))) ; 

// ==== flux �Υ�����Ȥ�Ԥ����ݤ��� flag��TRUE: ������Ȥ����====
const bool FLAG_FLUX_COUNT = false;//true;//

//************************************
//**  standard I/O
//************************************
// --- ɸ����Ϥ�������
const int N_COUT_INTERVAL   = 100;// 1;//

//************************************
//**  file I/O
//************************************

// === �����ե�������ɤ߹��फ��TRUE: �ɤ߹����=====
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

