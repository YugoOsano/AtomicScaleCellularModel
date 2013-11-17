 
#ifndef _ATOM_STRUCT_H_DEFINED_

#include <math.h>

//====================================
//==  utility constants       ========
//====================================
const int RANDOM_NUMBER_SEED = 529;//3 ;//4; //
//const int TRUE  = 1 ; 
//const int FALSE = 0 ;
 
//==================================
//==  条件コンパイルの設定
//==================================
// Cl ラジカルを入射する  <-@@@@ 重要 @@@@ 
#define _INJECT_CL_RADICAL_

//====== イオンの入射 ========
// -- 位置はランダム、垂直入射：         １
// -- 領域中央部からのみ（垂直に）入射 ：２
// -- 位置はランダム、IADF使用：         ３
const int ION_INJECT_FLAG  = 3 ;

const int N_IADF_ARRAY     = 901 ; // 上が３の場合の、IADFの区切り数

// ==== マスク無しの場合（yield の角度依存性等を計る場合）
//#define _NO_MASK_
 
// ==== yieldのイオン入射角依存性を考慮する場合 ====
const bool FLAG_INCIDENT_ANGLE = true;//false;//

// ==== (3rd paper)酸素が吸着することによる形状不安定を防ぐため、
//      酸化の影響（収率の低下）を分散させる
const bool FLAG_DISPERSE_OXIDATION = true;

#include <math.h>

// ===================================
// ==  物理定数               ========
// ===================================

const  double  PI    = 3.1415926535897932384626433 ; // 円周率
const  double  Q_ELEMENTAL    =  1.602e-19  ;    // 素電荷 ( c )
const  double  MASS_ELEC      =  9.1095e-31 ;   // 電子質量 ( kg )
const  double  E_PERMITTIVITY =  8.8542e-12  ;  // 真空の誘電率( C2 / J m )
const  double  K_BOLTZMANN    =  1.3806503e-23; // Boltzmann定数(m2 kg s-2 K-1)
const  double  K_BOLTZMANN_kcal = 1.987e-3  ;  // Boltzmann定数(kcal mol-1 K-1)

const  double  BOHR_RADIUS    = 5.29177e-11 ;  // Bohr半径 (m)



// --SCOFF を構造体定義する
struct Atom_struct
{
  int atomic_n ;   // 原子番号
  
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

// 式 L = N^(-1/3)  (p-118)による最も単純な自由行程の設定(m)
//                   ( 2.718592533e-10 m )
const double FREE_FLIGHT_PATH_Si = pow( (Si.density_atom * 1.0e+28),
                                        - 1.0 / 3.0 ) ;

// 原子間距離：自由行程と同じとする
const double L_INTER_ATOMIC      = FREE_FLIGHT_PATH_Si ;

//const struct Atom_struct Cl = {0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}  ;
const struct Atom_struct Cl = {17, 35,  34.969,  35.453, 
			       1.8956,  3.22,  0.70827,  0.97 }  ;

const double CL_MASS = Cl.atomic_weight * 1.6726231e-27  ;

//====================================
//==  modeling geometry     ==========
//====================================
// ion stopping の座標に合わせて、x-z 平面を扱うことにする(zは下向き)
// 2.718592533e-10 m / 1 cell
// 30nm  ゲートで、110.3512191 cells
// 50nm  ゲートで、183.9186983 cells
// 70nm  ゲートで、257.4861777 cells
// 100nm ゲートで、367.8373967 cells
// 200nm : 735 cells

const double CELL_SIZE = L_INTER_ATOMIC ;
const int    N_CELL_X  = 368 + 368 ;
//120;//184;//257;//110 111;//735;//368 上記のセル数に70足す @@@ 重要 @@@ ：
//マスクなしで実行する場合は全部で 184 とする
const int    N_CELL_Z  = 1200;//2570;//1000;// 

const int    N_TOP_MASK   = 10 + 184 ;// マスク上部の位置
const int    N_TOP_Si  = 184 + N_TOP_MASK;//735;// 368
//Si基板上部の位置（セル数：マスク厚さではなく、絶対的なZ座標）
const int    N_BOTTOM_Si = N_TOP_Si + 735;//N_CELL_Z ;//N_TOP_Si+368;//
//Si の膜厚を考慮する場合

const int    N_LEFT_MASK  = 184;//60;//35;//70;// マスク左壁面の位置
const int    N_RIGHT_MASK = N_CELL_X - 184;//60;//35;//70;// マスク右壁面の位置

const double SYSTEM_WIDTH_X  = L_INTER_ATOMIC * N_CELL_X ;
const double SYSTEM_HEIGHT_Z = L_INTER_ATOMIC * N_CELL_Z ;

//== マスク側壁の傾斜角 ==
const double SLOPE_ANGLE_SIDEWALL = 0.0 * PI / 180.0 ;//2.5 


//----------------------------------
//-- flux についての考察
// Sano's thesis では、イオン：Γi  = 1×10^16 /cm2 s (1×10^20 /m2 s)
//                   ラジカル：ΓCl = 2×10^18 /cm2 s (2×10^22 /m2 s)
// としている。ここでは、N_ION_INJECTがどれだけの時間に
// 相当するか計算しておく。

// このシミュレーションでの(x-y平面の)表面積は、
// SYSTEM_WIDTH_X * L_INTER_ATOMIC

// したがって、相当する時間は
// N_ION_INJECT / (Γi × SYSTEM_WIDTH_X × L_INTER_ATOMIC )
// 
// となる。

const double ION_FLUX = 1.0e+20 ;//(/m2 s)
//const double REAL_ELAPSED_TIME = 135.305 ;
//33.82625 ; //67.6525; //13.5305 ; //N_ION_INJECT 
// 135.305 としておくと、50nm, 1.0e+20 ion flux でちょうど になる




// @@@@@@  Cl -> Si の吸着確率  @@@@@@@
// @@@@@@  シミュレーション上での吸着確率は、これに
// @@@@@@  1 / Ymax をかけたものとする。
const  double  ADSORPTION_PROBABILITY = 1.0;//0.01;//0.1;//0.55 ;


// == １つのサイトに許されるCl adatom の数
const  int     NMAX_ADATOM  = 4 ;

// @@@@@@  Etch Yield の最大値：
// @@@@@@  エネルギー依存性を導入するため、粒子の吸着、
// @@@@@@  脱離などの反応は、1 / Ymax をかけた確率で起るとする
const  double  YIELD_MAX     = 5.0;//8.0;//3.0;//7.0;//4.0;
const  int     YIELD_MAX_INT = int (YIELD_MAX) ;

const  double  ENERGY_ION_STOPPING    = 4.0  ; // ( eV )

//=================================================
//==== Si -> SiO, SiO2 と酸化すると、etch yield は著しく
// === 小さくなる [Tuda, 1996]によると 50倍
// ここでは、単純に SiO/SiO2のyield が Si の x倍になると仮定する。
// この倍数を設定して、desorption の関数に導入する：一般の選択比の
// 逆数で表示する。

// おそらく50倍というのはbulkのSiO2で測定した選択比だと考えられる
// が、当面はそれを採用する （逆数は 0.02 ）。
// SiO では25倍とする（逆数は 0.04）。

const double SELECTIVITY_SiO  =  1.0/7.0;//0.05;//0.2;//0.1;//0.02 ; //0.04 ;
const double SELECTIVITY_SiO2 =  1.0/7.0;//0.05;//0.2;//0.1;//0.02 ;

//== ハードマスクの選択性：上記のSiO2を基準に決定
const double SELECTIVITY_HARDMASK =  0.05 /  SELECTIVITY_SiO2 ;// 0.05 / SELECTIVITY_SiO2 ;

//============================================================
// --- ファイルに記録する頻度（シミュレーション上の時間で計るように変更）
// int( 5.0  <- この数字が秒数を表す）
const int INTERVAL_FILE_OUTPUT 
=         int( 5.0 *   ION_FLUX * SYSTEM_WIDTH_X * L_INTER_ATOMIC ) ;
//=         int( 30.0 *   ION_FLUX * SYSTEM_WIDTH_X * L_INTER_ATOMIC ) ;
//const int INTERVAL_FILE_OUTPUT = 40000 ;//int(N_ION_INJECT / 1) ; //10) ; //


const int N_ION_INJECT = INTERVAL_FILE_OUTPUT * 10;//5;//15;//20;//2;//1;//3;//
// 10000000 ;//-- 注入イオン数
//=int(ION_FLUX * SYSTEM_WIDTH_X * 
//     L_INTER_ATOMIC * REAL_ELAPSED_TIME) ; 
const double REAL_ELAPSED_TIME 
=     N_ION_INJECT /  (ION_FLUX * SYSTEM_WIDTH_X * L_INTER_ATOMIC ) ;





/* ---------------------------------- 
  座標の取り方（速度には球面座標を用いる）
     
   \  |  /
    \ | /
     \|/) ψ
      \--------------------->z
      |\
      | \  z 軸回り:θ
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
  同じdirectory でshellによる連続実行ができるように
  引数(argc, argv)を導入

  archive: 03Feb2004
  Yield の入射エネルギー依存性を導入するために、
  Si原子の脱離を確率的にする (Y / Ymax) 前に
  アーカイブ
*/
