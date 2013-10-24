// neutral_particle.h

// 中性粒子のクラス。Particle class を継承する。

#ifndef  _NEUTRAL_PARTICLE_H_DEFINED_

#include "particle.h"

class Neutral_class : public Particle_class
{
protected:
  bool flag_adsorption ; //== 吸着を示すためのflag ;
  int  i_exception_reflect; //== 乱反射における例外処理のためのカウンタ

public:
  //-- constructor 
  Neutral_class() ;
  Neutral_class(double mass_input,  int flag_boundary_input);
  ~Neutral_class();

  //-- x座標は乱数、速度はx,z方向に与える
  // （y 方向に速度を与える必要はなし）

  // -- 向きはcosine distributionで決定
  //  -> このコードでは x軸方向が0°なので、sinθ
  // distribution を求めると、[-cosθ]θ〜0°/ [-cosθ]180°〜 0°
  // 0〜1の乱数R に対して
  // acos(1-2R)

  // 入力：並進速度
  void inject_from_top( double v ) ;

  //-- 右側からの入射
  void inject_from_right_side(double v) ;
  
  //==================================
  // 06 Aug 2004
  // ***  substrate surface における乱反射
  // monte_carlo.cc のmain 内で行っていたが、このクラス内で
  // 行うように refactoring 
 

  //=================================
  // 中性粒子の入射から放出、吸着まで
  // 入力：入射方向：領域上部(false)／右側(true)
  void all_process(class Shape_trim_class    *Shape,
		   class Shape_counter_class *Neutral_counter,
		   bool   flag_inject_from_side ,
		   bool   flag_flux_count ,
		   double incident_energy ,
		   bool   flag_chemical_etch,
		   double yield_chemical) ;
  

};


#define  _NEUTRAL_PARTICLE_H_DEFINED_
#endif
