// etch_product.h 
// 
// 2004/07/29 -
//  エッチング生成物のクラス：Neutral_particle 
//  クラスを継承する。deposition をメソッドとして実装する

#ifndef _ETCH_PRODUCT_H_DEFINED_

#include "neutral_particle.h"

//== 形状のセルの状態を判定するため、shapeクラスもincludeする。
#include "shape_trim.h"

class Etch_product_class : public Neutral_class
{
protected:
  //=== SiClx における塩素化の数 x ====
  int    n_chlorination ;

  //=== SiOy におけるoxidationの数 y ==== 
  int    n_oxidation ;

  //=== 吸着係数 Sp ===================
  double sticking_probability ;
  double sticking_probability_1st_contact ;

  //=== 現在のセルに入射した一つ前のセル番号
  int  i_previous_x, i_previous_z ;

public:
  //===実効的な入射確率
  //     パターン幅によって入射確率が変化する：
  // 基板全体のSiの露出がxパーセントとして、simulation domainの
  // Siの露出をyパーセントとすると、その領域では全体のy/x倍の放出がある。
  //  したがって、Pr × (x/y)が実効的な入射確率となる （ x = 0.5 とする）
  
  double reincidence_probability ;

public:
  //-- constructor 
  Etch_product_class() ;
  Etch_product_class(double mass_input, int flag_boundary_input,
		     bool   flag_open_space, 
		     int    n_chlorination_input ,
		     int    n_oxidation_input ,
		     double reincident_probability_whole,
		     double sticking_probability_input   );
  ~Etch_product_class();


  //-- deposition のメソッド
  // 単純な ballistic deposition model を実装する
  // --> 山崎、南部(2002, 2003) 機械学会誌を参照

  // 1. 固体部分(SHAPE_Si)に入射したか判定

  // =============================================
  // 現在の位置、速度ベクトルより、どのセルから
  // 現在のセルに入射したか計算する。

  //(1) v_x の正負によって場合分けを行う
  //(2) 速度ベクトルの逆方向を延長し、切片の位置から判定を行う

  // v_x > 0 ならば、0切片を取る
  //  -> 切片の求め方：(Δx, Δz)を通り、方向ベクトルが
  //   ( - v_x , - v_z) の直線の0切片は、
  //    Δz - (v_z/v_x)Δx 

  // v_x < 0 ならば、x = CELL_SIZE の切片を取る
  // 上と同じく、切片は
  //    Δz + (v_z/v_x)( CELL_SIZE - Δx)

  void get_previous_cell() ;


  //=========================================
  // （SiCl4の脱離を前提として）Si表面に粒子を配置し
  //初速を与える：下の whole_flight 関数の拡散反射の部分を抽出

  // 入力：Shape_trim_classのポインタ、
  //       位置／速度の構造体（位置情報のみ）、速度
  // 返り値：配置が行われた場合 TRUE そうでなければ FALSE

  bool allocate_on_surface(Shape_trim_class *Shape_pointer,
			   Particle_location_velocity_struct  pos_v_input,
			   double v ) ;

  // ===================================
  // クラス（Shape_trim_class）のポインタを利用して、
  // 粒子の入射以降の処理を関数化する

  // （吸着せずに）領域の外に出た場合、true を返す
  //  :  脱離した SiCl4 の数に応じて入射するSiCl2の数を決定するため
  void whole_flight(Shape_trim_class *Shape_pointer,
		    bool flag_inject_from_side ) ;

};




#define _ETCH_PRODUCT_H_DEFINED_
#endif

