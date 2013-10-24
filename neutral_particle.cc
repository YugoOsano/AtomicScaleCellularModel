
#include <iostream>
#include <math.h>
#include "neutral_particle.h"
#include "random_number.h"
#include "atom_struct.h"
#include "shape_trim.h"

Neutral_class::Neutral_class() 
{
}

Neutral_class::Neutral_class(double mass_input , int flag_boundary_input) 
  : Particle_class(mass_input, flag_boundary_input) 
{
  //mass = mass_input ; 
}
Neutral_class::~Neutral_class() 
{}


void Neutral_class::inject_from_top( double v )
{
  pos_v.x = RAN0() * SYSTEM_WIDTH_X ; 
  pos_v.y = 0.0 ;
  pos_v.z = 0.0 ;

  // acos(1-2R)
  //double theta = RAN0() * PI ;// obsolete: 
  double theta = acos(1.0 - 2.0 * RAN0()) ;

  pos_v.v_x = v * cos(theta) ;
  pos_v.v_y = 0.0 ;
  pos_v.v_z = v * sin(theta) ;

  pos_v.v_r     = v ;
  pos_v.v_theta = 0.0 ;
  pos_v.v_psi   = sin(theta) ;

	//	cout << pos_v.v_x << "\t" << pos_v.v_y << "\t"
     	//	<< pos_v.v_z << "\n" ;
  i_exception_move = 0 ;//== 例外処理用カウンタの初期化 ==
}
void Neutral_class::
inject_from_right_side(double v)
{
  pos_v.x = SYSTEM_WIDTH_X ; 
  pos_v.y = 0.0 ;
  pos_v.z = RAN0() * SYSTEM_HEIGHT_Z ;

  // acos(1-2R)
  //double theta = RAN0() * PI ;// obsolete: 
  double theta = acos(1.0 - 2.0 * RAN0()) ;

  pos_v.v_x = - v * sin(theta) ;
  pos_v.v_y =   0.0 ;
  pos_v.v_z =   v * cos(theta) ;

  pos_v.v_r     = v ;
  pos_v.v_theta = 0.0 ;
  pos_v.v_psi   = sin(theta) ;
}
//=================================
// 中性粒子の入射から放出、吸着まで
void Neutral_class::
all_process(class Shape_trim_class *Shape,
	    class Shape_counter_class *Neutral_counter,
	    bool   flag_inject_from_side ,
	    bool   flag_flux_count ,
	    double incident_energy,
	    bool   flag_chemical_etch,
	    double yield_chemical )
{
  //=== セルインデックス ===
  int i_x, i_z, i_x_particle, i_z_particle ;
  //--- 入射 ---
  inject_from_top
    ( sqrt(2.0 * incident_energy * 
	   Q_ELEMENTAL / mass) ) ; 

  i_exception_move = 0 ;
  
  //== 横方向からの入射 ==
  if(flag_inject_from_side == true)
    {
      inject_from_right_side
	( sqrt(2.0 * incident_energy * 
	       Q_ELEMENTAL / mass) ) ; 
      //== 入射の位置が空間で無ければ i_exception_move に大きい
      //   値を代入して while loop をskip する
      //   put_shape.. 関数は領域右端で0を返すので、L/2を引く
      if(Shape->put_shape_matrix(pos_v.x - 0.5 * CELL_SIZE, pos_v.z ,
				     &i_x, &i_z ) != SHAPE_SPACE )
	i_exception_move = N_CELL_Z * 200 + 1 ;
    }
 
  while(1) 
    {
      //-- 領域内部にあまり長くとどまっている場合、
      //   例外処理（粒子を消去）
      i_exception_move++ ;
      if( i_exception_move > N_CELL_Z * 200 )
	break ;
		
      // -- まず並進運動 ----
      move_trans(FREE_FLIGHT_PATH_Si) ;

      if (flag_inside == false) // -- 外に出た？
	break ;
	     
      //--- セルが空白でない？ ---
      //if(Shape.put_shape_matrix(pos_v.x, pos_v.z ,
      //			&i_x, &i_z ) != SHAPE_SPACE)
      if( Shape->find_solid_nearest_neighbor
	  (pos_v.x,    pos_v.z ,
	   pos_v.v_x,  pos_v.v_z ,
	   &i_x, &i_z, &i_x_particle, &i_z_particle ) != SHAPE_SPACE )
	{
	  //== flux counting 
	  if(flag_flux_count == true )
	    Neutral_counter->count(i_x, i_z ) ;
	  
	  //---- 化学的エッチング ----
	  // YIELD_CHEMICALの確率で脱離が起こる 
	  if( RAN0() < yield_chemical &&
	      flag_chemical_etch == true  )
	    {
	      Shape->desorb_Si_chemical_etch(i_x, i_z);
	    }
	  else //---- 吸着 -----
	    {
	      flag_adsorption = 
		Shape->settle_Cl_into_bareSi(i_x, i_z,  // 吸着した？
					    ADSORPTION_PROBABILITY );
	      if(flag_adsorption == true )
		break;   // 吸着したらループを出る
	    }
	  // 乱反射 : 
	  // 無限ループをさけるために、カウンタを用意する。
	  // 粒子を仮に進ませてみて、(pos_v_togo)
	  // もし進んだ先が空白部分でなければ、元に戻って
	  //     やり直し(50回になったらもうやめる)
	  i_exception_reflect = 0;
	  Shape->get_surfacenormal(i_x, i_z) ;
	reflection:
	  if( i_exception_reflect > 50 )
	    break ;
	  
	  random_reflection //(FALSE,0.0,0.0);//乱反射
	    (true,//<-cosine diffusion を導入すると左右非対称になる<- this bug was fixed (29Aug2004)
	     Shape->surfacenormal_x[i_x][i_z],
	     Shape->surfacenormal_z[i_x][i_z]) ;
		
	  move_trans_togo(FREE_FLIGHT_PATH_Si) ;
	  
	  if(Shape->put_shape_matrix(pos_v_togo.x, pos_v_togo.z ,
				     &i_x, &i_z ) != SHAPE_SPACE )
	    {
	      i_exception_reflect++ ;
	      goto reflection ;
	    }
	  // cout << Shape.surfacenormal_x[i_x][i_z] * 10000.0 
	  //	<< "\t" << Shape.surfacenormal_z[i_x][i_z] * 10000.0
	  //	<<  "\t" << Cl_neutral.pos_v.v_x << "\t" << Cl_neutral.pos_v.v_z << endl ;
	}
    }//=== while loop 終了
  
}
