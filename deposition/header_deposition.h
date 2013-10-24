
#ifndef _HEADER_DEPOSITION_INCLUDED_

#include "atom_struct.h"

//== マクロ設定（etch products/酸素をモデルに入れるか）

const bool FLAG_SiCl4_INCLUDED = true;//false;// <- basically true

const bool FLAG_SiCl2_INCLUDED = false;//true;//

const bool FLAG_OXYGEN_INCLUDED = false;//true;//

const bool FLAG_SiOCl2_INCLUDED = false;//true;//


//=== SiCl2 のflux : Tuda [1996] によると 
// Γp/Γi = 0, 0.1, 0.5
// としている。

// -> イオンを回数のカウンタにするのはそのままで、
// i_ion % [flux ratio] == 0 の時に SiCl2 を入射する
// （flux比は上の逆数 Γi/Γp で設定）

//const int ION_SiCl2_RATIO = 2;//1;//10 ;


// == 気相中に出たエッチング生成物がSiCl2に解離して
//    基板に戻ってくる確率

const  double  SiCl2_REINCIDENCE_PROBABILITY = 0.0;//0.5 ;//0.2;//0.05;//

// @@@@@@  SiCl4 の吸着確率  @@@@@@@
// Tuda, 1996 JVST B では、実際の吸着確率は
//   Sp <= 0.002 と記載されているが、
// より大きな値も仮定した場合の解析結果も提示されている
const  double  SiCl4_DEPOSITION_PROBABILITY = 0.02;//0.05;//0.1;//0.5;//0.2;//0.002;//0.01;//0.005;//0.001;

// @@@@@@  SiCl2 の吸着確率  @@@@@@@
// Tuda, 1996 JVST B では、吸着確率Spは 0〜1 の
// の中で変化させる、と記述されている。
//  Fig. 7 では Sp = 0.1 として、 Γp/Γi = 0, 0.1, 0.5
// としている。 
const  double  SiCl2_DEPOSITION_PROBABILITY = 0.1;//0.0;//0.5;//0.2;//


//=== O （酸素）のflux: Tuda [1996] によると 
// Γo/Γi = 0, 0.02, 0.05

const  double OXYGEN_ION_RATIO = 0.0;//1.0;//2.0;//0.5;//5.0;//0.1;//0.05;
//20.0;//50.0;//100.0;//30.0;//

// ΓSiO/Γi = 0, SiOCl2 の吸着確率
//const  double SiO_ION_RATIO = 0.0;//0.05 ;

const  double SiOCl2_REINCIDENCE_PROBABILITY = 0.0;//0.5;//0.2;//0.05;//
const  double SiOCl2_DEPOSITION_PROBABILITY  = 0.1;//0.5;//0.2; 

#define _HEADER_DEPOSITION_INCLUDED_
#endif
