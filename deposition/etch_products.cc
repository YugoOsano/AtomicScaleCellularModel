
// etch_products.cc
// 2004/07/29

#include <iostream>
#include <math.h>
#include "atom_struct.h"

#include "etch_products.h"


Etch_product_class::
Etch_product_class(double mass_input,  int flag_boundary_input,
		   bool   flag_open_space, 
		   int    n_chlorination_input ,
		   int    n_oxidation_input ,
		   double reincident_probability_whole,
		   double sticking_probability_input  ) 
  : Neutral_class(mass_input, flag_boundary_input)
{
  n_chlorination = n_chlorination_input  ;
  n_oxidation    = n_oxidation_input  ; 
  
  sticking_probability = sticking_probability_input ;
  
  if(flag_open_space == true)
    reincidence_probability = reincident_probability_whole
      *           ( 0.5 / ((double(N_CELL_X) - 184.0)/double(N_CELL_X)) );
  else
    reincidence_probability = reincident_probability_whole
      *           ( 0.5 / ((double(N_CELL_X) - 368.0)/double(N_CELL_X)) );

  std::cout << "effective Pr: " << reincidence_probability  << "\n" ;
}
Etch_product_class::~Etch_product_class()
{
}

void Etch_product_class::get_previous_cell()
{
  //(1) v_x > 0 
  if(pos_v.v_x > 0.0 )
    {
      //-- Δz - (v_z/v_x)Δx(0切片) < 0 ならば 
      //      上（zの負の方向）からの入射
      if( dz_cell - ( pos_v.v_z / pos_v.v_x ) * dx_cell < 0.0 )
	{
	  i_previous_x = i_cell_x ;
	  i_previous_z = i_cell_z - 1 ;
	}
      //-- Δz - (v_z/v_x)Δx(0切片) > cell size ならば 
      //      下（zの正の方向）からの入射
      else if(dz_cell - ( pos_v.v_z / pos_v.v_x ) * dx_cell >
	      CELL_SIZE)
	{
	  i_previous_x = i_cell_x ;
	  i_previous_z = i_cell_z + 1 ;
	}
      else // 横（xの負の方向）からの入射
	{
	  i_previous_x = i_cell_x - 1 ;
	  i_previous_z = i_cell_z  ;
	}
    }
  // v_x < 0 では、x = CELL_SIZE の切片
  // Δz + (v_z/v_x)( CELL_SIZE - Δx)
  else if(pos_v.v_x < 0.0) 
    {
      //-- (切片) < 0 ならば 
      //      上（zの負の方向）からの入射
      if( dz_cell + 
	  ( pos_v.v_z / pos_v.v_x ) * (CELL_SIZE - dx_cell) < 0.0 )
	{
	  i_previous_x = i_cell_x ;
	  i_previous_z = i_cell_z - 1 ;
	}
      //-- (切片) > cell size ならば 
      //      下（zの正の方向）からの入射
      else if(dz_cell +
	      ( pos_v.v_z / pos_v.v_x ) * (CELL_SIZE - dx_cell) >
	      CELL_SIZE)
	{
	  i_previous_x = i_cell_x ;
	  i_previous_z = i_cell_z + 1 ;
	}
       else // 横（xの正の方向）からの入射
	{
	  i_previous_x = i_cell_x + 1 ;
	  i_previous_z = i_cell_z  ;
	}


    }

  //=== 例外処理（domain から出ていないか？）
  if     ( i_previous_x < 0 )
    i_previous_x = 0 ;
  else if( i_previous_x > N_CELL_X - 1)
    i_previous_x = N_CELL_X - 1 ;
  
  if     ( i_previous_z < 0 )
    i_previous_z = 0 ;
  else if( i_previous_z > N_CELL_Z - 1)
    i_previous_z = N_CELL_Z - 1 ;
	
}

//========================================
// Si表面に粒子を配置し初速を与える
bool Etch_product_class::
allocate_on_surface(Shape_trim_class *Shape_pointer,
		    Particle_location_velocity_struct  pos_v_input,
		    double v )
{
  // random_reflection 関数は反射しても同じ速度を
  // 保つので最初に適切な速度を与える （後で位置のみ配置）

  inject_from_top( v ) ;

  pos_v.x = pos_v_input.x ;
  pos_v.y = pos_v_input.y ;
  pos_v.z = pos_v_input.z ;

  // -- 粒子のいる場所のphaseを取得
  //（セル番号 i_cell_x, i_cell_zも取得）
  get_position() ;

  //表面からの拡散 : 
  // 無限ループをさけるために、カウンタを用意する。
  // 粒子を仮に進ませてみて、(pos_v_togo)
  // もし進んだ先が空白部分でなければ、元に戻って
  //     やり直し(50回になったらもうやめる)
 
  int i_x, i_z ;
 
  int i_exception_reflect = 0;
  Shape_pointer->get_surfacenormal(i_cell_x, 
				   i_cell_z) ;
 
 diffusion_SiClx:
  if( i_exception_reflect > 50 )
    return false ;
  
  random_reflection 
    (true,  // 乱反射
     Shape_pointer->surfacenormal_x[i_cell_x][i_cell_z],
     Shape_pointer->surfacenormal_z[i_cell_x][i_cell_z]) ;
  
  move_trans_togo(FREE_FLIGHT_PATH_Si) ;
  
  if( Shape_pointer->put_shape_matrix(pos_v_togo.x, 
				      pos_v_togo.z ,
				      &i_x, &i_z ) != SHAPE_SPACE )
    {
      i_exception_reflect++ ;
      goto diffusion_SiClx ;
    }
  return true ;

}

//=======================================================
// 粒子の入射以降の処理を関数化する
void Etch_product_class::
whole_flight(Shape_trim_class *Shape_pointer,
	     bool flag_inject_from_side) 
{
  bool flag_deposition    ;
  int  i_exception_reflect ;
  int  i_x, i_z ;

  i_exception_move = 0 ;

  //== 横方向からの入射 : neutral_particle.cc のall_process関数とほぼ同じ
  if(flag_inject_from_side == true)
    {
      inject_from_right_side(pos_v.v_r);
      //== 入射の位置が空間で無ければ i_exception_move に大きい
      //   値を代入して while loop をskip する
      if(Shape_pointer->put_shape_matrix(pos_v.x - 0.5 * CELL_SIZE, pos_v.z ,
					 &i_x, &i_z ) != SHAPE_SPACE )
	i_exception_move = N_CELL_Z * 200 + 1 ;
    }
  
  while(1) 
    {
      //--- 領域内部にあまり長くとどまっている場合、
      //    例外処理（粒子を消去）
      i_exception_move++ ;
      if( i_exception_move > N_CELL_Z * 200 )
	break ;
      
      // -- 並進運動
      move_trans(FREE_FLIGHT_PATH_Si) ;
      
      if (flag_inside == false) // -- 外に出た？
	break ;
      
      // -- 粒子のいる場所のphaseを取得
      get_position() ;
	      
      // -- 固相であれば、一つ前のセル番号を取得 -> deposit
      if(Shape_pointer->shape_matrix
	 [i_cell_x][i_cell_z] != SHAPE_SPACE )
	{
	  get_previous_cell() ; // セル番号取得
	  flag_deposition = 
	    Shape_pointer->deposit_Si(i_previous_x,     // deposit
				      i_previous_z,
				      i_cell_x,    i_cell_z , 
				      n_chlorination, n_oxidation,
				      sticking_probability ) ;	  
	  if(flag_deposition == true )
	  {
	    // cout << "SiClx deposited: \t" << i_previous_x 
	     //<< "\t" << i_previous_z << endl ;
	    break ;  // get out of the loop if deposited
	  }
	  //乱反射 : 
	  // 無限ループをさけるために、カウンタを用意する。
	  // 粒子を仮に進ませてみて、(pos_v_togo)
	  // もし進んだ先が空白部分でなければ、元に戻って
	  //     やり直し(50回になったらもうやめる)
	  i_exception_reflect = 0;
	  Shape_pointer->get_surfacenormal(i_cell_x, 
					   i_cell_z) ;
	reflection_SiClx:
	  if( i_exception_reflect > 50 )
	    break ;
	  
	  random_reflection //(FALSE, 0.0, 0.0 ) ;
	    (true,  // 乱反射
	     Shape_pointer->surfacenormal_x[i_cell_x][i_cell_z],
	     Shape_pointer->surfacenormal_z[i_cell_x][i_cell_z]) ;
	  
	  move_trans_togo(FREE_FLIGHT_PATH_Si) ;
	  
	  if( Shape_pointer->put_shape_matrix(pos_v_togo.x, 
					      pos_v_togo.z ,
					      &i_x, &i_z ) != SHAPE_SPACE )
	    {
	      i_exception_reflect++ ;
	      goto reflection_SiClx ;
	    }
	} 
    } //=== while loop 終了（１粒子のtrajectory）
}
