 
#ifndef _ATOM_STRUCT_H_DEFINED_

#include <math.h>

//====================================
//==  utility constants       ========
//====================================
const int RANDOM_NUMBER_SEED = 529;//3 ;//4; //
//const int TRUE  = 1 ; 
//const int FALSE = 0 ;
 
//==================================
//==  ��拾��ѥ��������
//==================================
// Cl �饸��������ͤ���  <-@@@@ ���� @@@@ 
#define _INJECT_CL_RADICAL_

//====== ����������� ========
// -- ���֤ϥ����ࡢ��ľ���͡�         ��
// -- �ΰ����������Τߡʿ�ľ�ˡ����� ����
// -- ���֤ϥ����ࡢIADF���ѡ�         ��
const int ION_INJECT_FLAG  = 3 ;

const int N_IADF_ARRAY     = 901 ; // �夬���ξ��Ρ�IADF�ζ��ڤ��

// ==== �ޥ���̵���ξ���yield �γ��ٰ�¸������פ����
//#define _NO_MASK_
 
// ==== yield�Υ��������ͳѰ�¸�����θ������ ====
const bool FLAG_INCIDENT_ANGLE = true;//false;//

// ==== (3rd paper)���Ǥ����夹�뤳�Ȥˤ������԰�����ɤ����ᡢ
//      �����αƶ��ʼ�Ψ���㲼�ˤ�ʬ��������
const bool FLAG_DISPERSE_OXIDATION = true;

#include <math.h>

// ===================================
// ==  ʪ�����               ========
// ===================================

const  double  PI    = 3.1415926535897932384626433 ; // �߼�Ψ
const  double  Q_ELEMENTAL    =  1.602e-19  ;    // ���Ų� ( c )
const  double  MASS_ELEC      =  9.1095e-31 ;   // �ŻҼ��� ( kg )
const  double  E_PERMITTIVITY =  8.8542e-12  ;  // ������Ͷ��Ψ( C2 / J m )
const  double  K_BOLTZMANN    =  1.3806503e-23; // Boltzmann���(m2 kg s-2 K-1)
const  double  K_BOLTZMANN_kcal = 1.987e-3  ;  // Boltzmann���(kcal mol-1 K-1)

const  double  BOHR_RADIUS    = 5.29177e-11 ;  // BohrȾ�� (m)



// --SCOFF ��¤���������
struct Atom_struct
{
  int atomic_n ;   // �����ֹ�
  
  // Atomic mass of most abundant isotope
  int atomic_mass ;

  // Atomic weight of most abundant isotope 
  double atomic_weight ;
  // Atomic weight of solid with normal isotopic abundance
  double atomic_weight_solid ;

  // Density of solid in grams/cm3
  double density_mass ;
  // Density of solid in units of 10^22 atoms/cm3
  double density_atom ;

  // Fermi velocity of solid
  double fermi_v ;

  // Ion screening length factor
  double factor_ion_screening;
  
};

const struct Atom_struct Si = {14, 28,  27.977,  28.086,
			       2.3212,  4.977, 0.97411,  0.88 }  ;

// �� L = N^(-1/3)  (p-118)�ˤ��Ǥ�ñ��ʼ�ͳ����������(m)
//                   ( 2.718592533e-10 m )
const double FREE_FLIGHT_PATH_Si = pow( (Si.density_atom * 1.0e+28),
                                        - 1.0 / 3.0 ) ;

// ���Ҵֵ�Υ����ͳ������Ʊ���Ȥ���
const double L_INTER_ATOMIC      = FREE_FLIGHT_PATH_Si ;

//const struct Atom_struct Cl = {0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}  ;
const struct Atom_struct Cl = {17, 35,  34.969,  35.453, 
			       1.8956,  3.22,  0.70827,  0.97 }  ;

const double CL_MASS = Cl.atomic_weight * 1.6726231e-27  ;

//====================================
//==  modeling geometry     ==========
//====================================
// ion stopping �κ�ɸ�˹�碌�ơ�x-z ʿ�̤򰷤����Ȥˤ���(z�ϲ�����)
// 2.718592533e-10 m / 1 cell
// 30nm  �����Ȥǡ�110.3512191 cells
// 50nm  �����Ȥǡ�183.9186983 cells
// 70nm  �����Ȥǡ�257.4861777 cells
// 100nm �����Ȥǡ�367.8373967 cells
// 200nm : 735 cells

const double CELL_SIZE = L_INTER_ATOMIC ;
const int    N_CELL_X  = 368 + 368 ;
//120;//184;//257;//110 111;//735;//368 �嵭�Υ������70­�� @@@ ���� @@@ ��
//�ޥ����ʤ��Ǽ¹Ԥ������������ 184 �Ȥ���
const int    N_CELL_Z  = 1200;//2570;//1000;// 

const int    N_TOP_MASK   = 10 + 184 ;// �ޥ��������ΰ���
const int    N_TOP_Si  = 184 + N_TOP_MASK;//735;// 368
//Si���ľ����ΰ��֡ʥ�������ޥ��������ǤϤʤ�������Ū��Z��ɸ��
const int    N_BOTTOM_Si = N_TOP_Si + 735;//N_CELL_Z ;//N_TOP_Si+368;//
//Si ��������θ������

const int    N_LEFT_MASK  = 184;//60;//35;//70;// �ޥ��������̤ΰ���
const int    N_RIGHT_MASK = N_CELL_X - 184;//60;//35;//70;// �ޥ��������̤ΰ���

const double SYSTEM_WIDTH_X  = L_INTER_ATOMIC * N_CELL_X ;
const double SYSTEM_HEIGHT_Z = L_INTER_ATOMIC * N_CELL_Z ;

//== �ޥ���¦�ɤη��г� ==
const double SLOPE_ANGLE_SIDEWALL = 0.0 * PI / 180.0 ;//2.5 


//----------------------------------
//-- flux �ˤĤ��Ƥιͻ�
// Sano's thesis �Ǥϡ������󡧦�i  = 1��10^16 /cm2 s (1��10^20 /m2 s)
//                   �饸���롧��Cl = 2��10^18 /cm2 s (2��10^22 /m2 s)
// �Ȥ��Ƥ��롣�����Ǥϡ�N_ION_INJECT���ɤ�����λ��֤�
// �������뤫�׻����Ƥ�����

// ���Υ��ߥ�졼�����Ǥ�(x-yʿ�̤�)ɽ���Ѥϡ�
// SYSTEM_WIDTH_X * L_INTER_ATOMIC

// �������äơ�����������֤�
// N_ION_INJECT / (��i �� SYSTEM_WIDTH_X �� L_INTER_ATOMIC )
// 
// �Ȥʤ롣

const double ION_FLUX = 1.0e+20 ;//(/m2 s)
//const double REAL_ELAPSED_TIME = 135.305 ;
//33.82625 ; //67.6525; //13.5305 ; //N_ION_INJECT 
// 135.305 �Ȥ��Ƥ����ȡ�50nm, 1.0e+20 ion flux �Ǥ��礦�� �ˤʤ�




// @@@@@@  Cl -> Si �ε����Ψ  @@@@@@@
// @@@@@@  ���ߥ�졼������Ǥε����Ψ�ϡ������
// @@@@@@  1 / Ymax �򤫤�����ΤȤ��롣
const  double  ADSORPTION_PROBABILITY = 1.0;//0.01;//0.1;//0.55 ;


// == ���ĤΥ����Ȥ˵������Cl adatom �ο�
const  int     NMAX_ADATOM  = 4 ;

// @@@@@@  Etch Yield �κ����͡�
// @@@@@@  ���ͥ륮����¸����Ƴ�����뤿�ᡢγ�Ҥε��塢
// @@@@@@  æΥ�ʤɤ�ȿ���ϡ�1 / Ymax �򤫤�����Ψ�ǵ���Ȥ���
const  double  YIELD_MAX     = 5.0;//8.0;//3.0;//7.0;//4.0;
const  int     YIELD_MAX_INT = int (YIELD_MAX) ;

const  double  ENERGY_ION_STOPPING    = 4.0  ; // ( eV )

//=================================================
//==== Si -> SiO, SiO2 �Ȼ�������ȡ�etch yield ��������
// === �������ʤ� [Tuda, 1996]�ˤ��� 50��
// �����Ǥϡ�ñ��� SiO/SiO2��yield �� Si �� x�ܤˤʤ�Ȳ��ꤹ�롣
// �����ܿ������ꤷ�ơ�desorption �δؿ���Ƴ�����롧���̤��������
// �տ���ɽ�����롣

// �����餯50�ܤȤ����Τ�bulk��SiO2��¬�ꤷ����������ȹͤ�����
// �������̤Ϥ������Ѥ��� �ʵտ��� 0.02 �ˡ�
// SiO �Ǥ�25�ܤȤ���ʵտ��� 0.04�ˡ�

const double SELECTIVITY_SiO  =  1.0/7.0;//0.05;//0.2;//0.1;//0.02 ; //0.04 ;
const double SELECTIVITY_SiO2 =  1.0/7.0;//0.05;//0.2;//0.1;//0.02 ;

//== �ϡ��ɥޥ��������������嵭��SiO2����˷���
const double SELECTIVITY_HARDMASK =  0.05 /  SELECTIVITY_SiO2 ;// 0.05 / SELECTIVITY_SiO2 ;

//============================================================
// --- �ե�����˵�Ͽ�������١ʥ��ߥ�졼������λ��֤Ƿפ�褦���ѹ���
// int( 5.0  <- ���ο������ÿ���ɽ����
const int INTERVAL_FILE_OUTPUT 
=         int( 5.0 *   ION_FLUX * SYSTEM_WIDTH_X * L_INTER_ATOMIC ) ;
//=         int( 30.0 *   ION_FLUX * SYSTEM_WIDTH_X * L_INTER_ATOMIC ) ;
//const int INTERVAL_FILE_OUTPUT = 40000 ;//int(N_ION_INJECT / 1) ; //10) ; //


const int N_ION_INJECT = INTERVAL_FILE_OUTPUT * 10;//5;//15;//20;//2;//1;//3;//
// 10000000 ;//-- �����������
//=int(ION_FLUX * SYSTEM_WIDTH_X * 
//     L_INTER_ATOMIC * REAL_ELAPSED_TIME) ; 
const double REAL_ELAPSED_TIME 
=     N_ION_INJECT /  (ION_FLUX * SYSTEM_WIDTH_X * L_INTER_ATOMIC ) ;





/* ---------------------------------- 
  ��ɸ�μ������®�٤ˤϵ��̺�ɸ���Ѥ����
     
   \  |  /
    \ | /
     \|/) ��
      \--------------------->z
      |\
      | \  z �����:��
*/
struct Particle_location_velocity_struct
{
  double x ;
  double y ;
  double z ; 
  
  double v_x ;
  double v_y ;
  double v_z ;

  double v_r ;
  double v_theta ;
  double v_psi ;

  /*double v_x ;
  double v_y ;
  double v_z ;
  */
};


#define _ATOM_STRUCT_H_DEFINED_
#endif

//===========================================
// update record
//===========================================

/*
  
  17Feb2004
  Ʊ��directory ��shell�ˤ��Ϣ³�¹Ԥ��Ǥ���褦��
  ����(argc, argv)��Ƴ��

  archive: 03Feb2004
  Yield �����ͥ��ͥ륮����¸����Ƴ�����뤿��ˡ�
  Si���Ҥ�æΥ���ΨŪ�ˤ��� (Y / Ymax) ����
  ����������
*/
