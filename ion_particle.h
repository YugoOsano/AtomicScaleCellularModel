// ion_particle.h

// イオンのクラス：Particle class を継承

// 19Oct2004
// SiCl4 の脱離を考慮するためのフラグをここで導入


#ifndef  _ION_PARTICLE_H_DEFINED_

#include <iostream>
#include <fstream>
#include "particle.h"
#include "shape_trim.h"

class Ion_class : public Particle_class
{
  
public:
  //== イオンの衝突によって表面のSi原子が脱離したかを示すカウンタ ==
  // （イオンの入射のたびにfalseに戻す）
  bool  flag_desorb_Si ;
  //== 脱離した生成物中の酸素の数 (SiClxOyのy)
  int   n_oxy_desorb_Si ; 
  
  //== 脱離が起こった時のイオンの座標: 構造体で記録しておく（速度は不使用）
  // （SiCl4のflightがここから始まるとする）
  struct Particle_location_velocity_struct position_at_desorption ;

protected:
  //========================================================
  //    イオンの前方散乱に関する記録を行う
  bool flag_reflection ; // 前方散乱が起こったかどうかのflag
  int  ctr_reflection  ; // カウンタ
  std::ofstream velocity_file ; //== 散乱後の速度を記録するためのファイル

protected:
  //========================================================
  //    イオンが最表面に衝突した後、反射するか／基板上に侵入してstopするか
  //    判断し、yieldを区別する場合に用いる変数
  //    あらかじめ衝突時の位置等を記録し、事後的にエッチングの処理を行う
  //==  衝突（エッチング）位置、エネルギーを記録する変数  ==
  int    i_x_etch, i_z_etch ;
  double energy_etch ;
  struct Particle_location_velocity_struct pos_v_etch ;

public:
  void  record_etch_position(int i_x, int i_z, double input_energy)
  {
    i_x_etch    = i_x ;          i_z_etch = i_z ;
    energy_etch = input_energy ; pos_v_etch = pos_v ; 
  }

public:
  //-- constructor 
  Ion_class() ;
  Ion_class(double mass_input, int flag_boundary_input,
	    int n_adf_input );
  ~Ion_class();

  //== 脱離カウンタのリセット ==
  //void reset_ctr_desorb_Si(){ ctr_desorb_Si = 0 ;}
  
  //== 前方散乱関連 ==
  void set_flag_reflection(bool input_flag) //== flagの設定 ==
  { flag_reflection = input_flag ; }
  
  void add_ctr_reflection()//== カウンタの加算 ==
  { if(flag_reflection == true) ctr_reflection++ ;  }

  int  put_ctr_reflection()//== カウンタの出力 ==
  { return  ctr_reflection ; }

  //=============================================
  //*********************************************************
  // -- microstructure/charged_particle::inject_from_top 関数
  // を、改良して使用する。ここでは、まずangular distribution のみ
  // 導入する。また、投入する粒子は１つとする。
  void inject_iadf( //int n_injection, 
		   //	 double energy_array[],double energy_threshold,
		   //	 double energy_df[],
		   double v_z) ;  
  // double angle_array[],  double angle_df[] );


  //==============================================
  //== 脱離が起こった時にフラグを立てて位置を記録
  //   入力：（エッチング時点の）位置／速度, 酸素の数
  void record_desorption(Particle_location_velocity_struct pos_v_recorded,
			 int input_n_oxygen) ;


  //==============================================
  //== イオンの衝突によるエッチングの処理
  //   イオンが反射した場合のエッチング収率は Er < Eth であれば C(√Ei - √Eth),
  //                                          Er > Eth であれば C(√Ei - √Er)とする
  //   反射で無い場合は Er = 0 を入力
  //   入力：Shape_trim class 
  //   （エッチング時点の）入射エネルギー Ei, 反射エネルギー Er,
  //   位置／速度、 格子点のindex (i_x, i_z)
  void ion_enhanced_etch(class  Shape_trim_class *Shape_trim,
			 double incident_energy, double reflected_energy,
			 Particle_location_velocity_struct pos_v_recorded,
			 int i_x, int i_z ) ;

  //== ハードマスクのスパッタリング
  //   収率は酸化されたSi基板のものにさらに係数をかける
  void hardmask_sputter(class  Shape_trim_class *Shape_trim,
			double incident_energy, double reflected_energy,
			Particle_location_velocity_struct pos_v_recorded,
			int i_x, int i_z ) ;

  //==============================================
  //==   パターン表面に到達したかどうかの判定 
  //   ->エッチング、衝突の処理
  //  返り値：true: 領域を出る、もしくはイオンがstopする
  //   （while loop を出る）
  bool impact_on_surface(class  Shape_trim_class *Shape_trim) ;

  //==============================================
  //==   上と同じく衝突の処理：前方散乱のも含める

  //     1 イオンのいる格子点がSiであるか？ 
  //        YES -> 従来通りの衝突の処理 
  //  NO -> 2  隣接格子点がSi/Hard mask であるか？
  //        YES -> impact parameterを正確に計算する
  //               (collision_accutate)
   bool impact_scattering(class  Shape_trim_class *Shape_trim,
			  bool  flag_mask_erosion) ;

};


#define  _ION_PARTICLE_H_DEFINED_
#endif
