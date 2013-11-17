
// particle.h

// 粒子もクラスとして扱う。

#ifndef _PARTICLE_H_DEFINED_

#include "atom_struct.h"

//== 境界条件のためのflag == 
const int NO_PERIODIC      = 1;//== 周期境界なし ==
const int SPECULAR_REFLECT = 2;//== 境界で鏡面反射 ==

//
class Particle_class
{
public:
  struct Particle_location_velocity_struct pos_v ;
  //--- 進ませた場合の仮の位置
  struct Particle_location_velocity_struct pos_v_togo ;

  //===== 2004/08/04 追加
  //  その粒子が位置するセル番号、セル内の相対的な位置
  // Shape_trim_class::put_shape_matrix で行っていた計算
  // の一部分である。
  /*     dx
       +--------+
       |   |    |
     dz|   |    |
       |---*    |
       |        |
       +--------+
  */
  int    i_cell_x, i_cell_z ;
  double dx_cell , dz_cell  ;


  //--- 質量
  double mass ;
protected:
  //---運動エネルギー
  double energy_k ;
public:


  //=========================================
  //===  入射角度分布 =======================
  int n_adf_array ;  // 角度分布関数の階級数→ 1 を指定すると等方入射とする
  double *angular_df ;
  double *angle_array_degree ;
  
  //--- 表面に侵入しているかどうかの flag
  // これをもとに脱離を行うかどうか判断する。
  // -- 初期値は FALSE, inject でFALSE に戻す。
  bool   flag_contact ; 

protected: 
  //=============================
  //===  入射粒子数のカウンタ ===
  int n_injection ;

  //== 例外処理のためのカウンタ
  int i_exception_move ;

  //== 領域内部にあるか判定するためのフラグ（内部：true）
  bool flag_inside ;
  bool flag_inside_togo ;

  //== x方向の境界条件を表すフラグ === 
  // （0：周期境界条件／1：周期境界なし／2：境界で鏡面反射）
  int  flag_boundary ;

public:
  //-- constructor
  Particle_class() ;
  Particle_class(double mass_input, int flag_boundary_input ) ;
  ~Particle_class() ;

  int  put_n_injection(){ return n_injection ;}
  void add_n_injection(){ n_injection++ ;}
  
  bool put_flag_inside(){ return flag_inside ;}

  //== flag_boundary の変更 ==
  void set_flag_boundary(int input_flag) 
  { flag_boundary = input_flag ; }

  //=======================================
  //== 角度分布をファイルから取得する =====
  // from "flux_model" directory
  // input : file name
  void read_angular_df(char filename[]);

  //-- 粒子入射：座標はx 軸の中心、速度はz方向のみ与える
  void inject_from_center( double v_z ) ;

  //-- 粒子入射：x座標は乱数、速度はz方向のみ与える
  void inject_from_top(double v_z) ;

  //-- 粒子入射：斜め方向に初期速度を与える
  // θの範囲は -π/2 〜 π/2 
  //  入力：位置、入射角、速度
  void inject_oblique(double x, double psi, double v) ;


  //===== 2004/08/04 追加
  //  その粒子が位置するセル番号、セル内の相対的な位置
  // ( i_cell_x, i_cell_z , dx_cell, dz_cell ) を計算する
  void get_position() ;

  //=================================
  //-- 長さ l だけ並進運動させる
  //     領域から出た場合 flag_inside をfalseにする
  //
  // -- 周期境界条件を入力
  //================================
  void  move_trans(double  l) ;
    
  // -- pos_v_togo を求める
  void  move_trans_togo(double  l ) ;

  //=== 06 Aug 2004 modified 
  // 基板表面における thermal reflection
  //-- 
  // flag FALSE -> 単純な速度の変更 
  //      TRUE  -> cosine distribution の導入
  // （法線方向も入力）
  // 法線方向を０°としてdistribution を求めると、
  // 0〜1 の乱数R に対して
  // asin(2R-1)

  void random_reflection(bool   flag_cosine_dist ,
			 double normal_x, double normal_z ) ;

  //-------------------------
  //-- nuclear energy loss (T) を出力する
  // 入力：入射イオン、ターゲット固体、入射エネルギー(eV)、
  //                   impact parameter (p)
  // 出力：scattering angle

  double put_energy_loss(Atom_struct  Ion,   Atom_struct  Solid,
			 double  energy_i ,
			 double  b_impact_para,
			 double  *theta_c ,    double  *psi_scattering ) ;

  //--------------------------
  // 一回のイオン／原子の衝突
  // 返り値：跳ね返った原子のエネルギー: (tmp_energy)
  void  collision_with_solid_atom(double p_impact_parameter ) ;

  //==========================
  // イオン／原子の衝突
  // イオンの前方散乱を表すため、impact parameterを
  // 粒子の位置、速度から計算して求める。
  // 簡略化のため、２次元の計算のみとする。
  // （奥行き方向の速度は一定とする）

  // 入力値：固体原子の位置
  void  collision_accurate(double x_solid, double z_solid,
			   double *tmp_energy ) ;

  // イオン／原子の衝突：3次元計算

  /* イオンの位置を(x, y, z), 速度を(vx, vy, vz), 
     Si原子の位置を(X, Y, Z) として、イオン-原子の最接近の位置を
     求める。パラメータλ、ベクトルr = (rx, ry, rz) を利用して、
     
     (x, y, z)+λ(vx, vy, vz) = (X, Y, Z)+(rx, ry, rz) 
     を解く。
     λ = [(X-x)vx + (Y-y)vy + (Z-z)vz]/(vx^2 + vy^2 + vz^2)

     rx = x - X + λ vx
     ry = y - Y + λ vy
     rz = z - Z + λ vz 

     impact parameter は |r|

     散乱後の速度ベクトルは v = (vx, vy, vz) と r
     の作る平面上にあるので、これらの一次結合。
     それぞれ単位ベクトル[v],[r]を求め、散乱角がψであるとすると、
     v' = |v'|([v] cosψ + [r] sinψ)
     となる。

     位置ベクトルはSi原子より、v'と垂直の方向に |r|だけ離れた
     地点 (X, Y, Z) + |r|([v](-sinψ) + [r]cosψ) とする。
  */
  void  collision_accurate3D(double x_solid, double y_solid,
			     double z_solid, double *tmp_energy ) ;

} ;



#define _PARTICLE_H_DEFINED_
#endif

