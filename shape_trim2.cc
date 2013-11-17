// shape_trim2.cc
//  
//   shape_trim.cc が大きくなったので分割

#include <iostream>
#include <fstream>
#include <stddef.h>
#include <math.h>


#include "shape_trim.h"
#include "random_number.h"
#include "fileio_particle.h"
#include "common_utils.h"

//=====================================================
//== chemical etching
void Shape_trim_class::desorb_Si_chemical_etch(int i_x, int i_z)
{
  if(shape_matrix[i_x][i_z] == SHAPE_Si &&
     n_Clbond[i_x][i_z]     >= 4        &&//== 吸着Clが4以上
     n_oxygen[i_x][i_z]     == 0  )  //== 酸素の吸着がなければ
    {
      shape_matrix[i_x][i_z] = SHAPE_SPACE ;
      n_Clbond[i_x][i_z]     = 0 ;
      n_oxygen[i_x][i_z]     = 0 ;

      // 脱離Si のカウンタ
      cntr_desorbed_Si++ ;
      //  std::cout << i_x << "\t" << i_z << "\tchemical_etch!\n" ;
    }
}

//=====================================================
//==   SiClx(x = 1〜4)の脱離の関数

void Shape_trim_class::
desorb_SiClx(int i_x, int i_z) 
{
  //==  SiClx の脱離の確率はSiCl4の x/4 とする
  if(shape_matrix[i_x][i_z] == SHAPE_Si &&
     RAN0() < double(n_Clbond[i_x][i_z]) / double(NMAX_ADATOM))
    {
      shape_matrix[i_x][i_z] = SHAPE_SPACE ;
      n_Clbond[i_x][i_z]     = 0 ;
      n_oxygen[i_x][i_z]     = 0 ;

      // 脱離Si のカウンタ
      cntr_desorbed_Si++ ;
           
      if(FLAG_INCIDENT_ANGLE == true)
	get_surfacenormal_around(i_x, i_z) ;

    }
  // -- domain の下から２番目まで到達したら、flag_get_bottom を
  //    TRUE にする（一番下は hard mask にする）
  if( i_z >= N_CELL_Z - 1)
    flag_get_bottom  =  true  ;
  
}

double etch_yield_disperse_oxidation(double theta, int n_oxygen )
{
  double tmp_coef = 0.0; 

  //===酸化膜の場合：Mahorowala, Sawin JVST 20, 1064 (2002)
  //   paper では最大値を元に正規化されているが、
  //   ここでは θ=0 を基準にするため、0.24 を割る
  if(n_oxygen >= 1)
    {
      tmp_coef = 
	0.4 * (18.7  * cos(theta)          -  64.7 * pow(cos(theta),2.0) +
	       145.2 * pow(cos(theta),3.0) - 206.0 * pow(cos(theta),4.0) + 
	       147.3 * pow(cos(theta),5.0) -  39.9 * pow(cos(theta),6.0))/0.24;

      if(tmp_coef < 0.0)
	tmp_coef = 0.0 ;
      if(n_oxygen > 6)
	n_oxygen = 6 ;

      tmp_coef = ( n_oxygen / 6.0 ) * tmp_coef * SELECTIVITY_SiO2 
	+        ( 1.0 - n_oxygen / 6.0 ) ;

      //std::std::cout << theta << "\tCoef: " << tmp_coef << "\n" ;
    }
  else //== poly-Siの場合 ==
    {
      tmp_coef = 1.0 ;
    }
  return  tmp_coef ;
}



//***************************************
// Cl density distribution - depthの統計をとる。

void Shape_trim_class::output_Cl_density_depth(char file_name[],
					       int  time_step)
{
  //-- 統計をとるためのカウンタ（浅い順に 0 - CL_DEPTH_RANGE）
  int *counter_Cl ;
  counter_Cl = new int[CL_DEPTH_RANGE] ;
 
  for (int i = 0 ; i < CL_DEPTH_RANGE ; i++)
    counter_Cl[i] = 0 ; //-- 初期化

  bool flag_below_surface ; // 最表面よりも下であるかどうかのflag
  int  tmp_depth ;          // 最表面からの深さ

  for(int i_x = 0; i_x < N_CELL_X ; i_x++)
    {
      flag_below_surface = false ;
      tmp_depth          = 0     ;

      for(int i_z = 0; i_z < N_CELL_Z ; i_z++)
	{
	  if(shape_matrix[i_x][i_z] == SHAPE_Si)
	    flag_below_surface = true ;

	  if(flag_below_surface == true &&
	     tmp_depth < CL_DEPTH_RANGE )
	    {
	      counter_Cl[tmp_depth] += n_Clbond[i_x][i_z] ;
	      tmp_depth++ ;
	    }
	}
    }
  // -- ファイル出力
  char OUT1[50] ; 
  sprintf( OUT1, "%s%d.dat", file_name, time_step);
  std::ofstream  cl_density_depth_file(OUT1) ;
  for(int i = 0 ; i < CL_DEPTH_RANGE ; i++)
    cl_density_depth_file << i << "\t" 
			  << counter_Cl[i] << std::endl ;

  delete [] counter_Cl ;
}

//--- 形状ファイルの読み込み
//    入力：file name(Si, Cl)
void Shape_trim_class::
input_profile(char file_Si[],char file_Cl[]) 
{
  input_array_splot(file_Si, shape_matrix ,
		     N_CELL_X, N_CELL_Z ) ;

  input_cellinfo(n_oxygen, n_Clbond , file_Cl, 
		 N_CELL_X, N_CELL_Z ) ;
}
	  
//----------------------------
void  Shape_trim_class::remove_isolated_Si() 
{
  // -- 
  int i_x, i_z ;

  for( i_z = 2; i_z < N_CELL_Z - 2 ; i_z++ )
    {
      for( i_x = 1; i_x < N_CELL_X - 1; i_x++ )
	{
	  //-- そのセルはSi 原子があるか？
	  if(shape_matrix[i_x][i_z] == SHAPE_Si)
	    {
	      //-- 周囲の８セルがすべて空白であれば、
	      //    Si原子を取り除く
	      if(shape_matrix[i_x - 1][i_z - 1] == SHAPE_SPACE &&
		 shape_matrix[i_x - 1][i_z    ] == SHAPE_SPACE &&
		 shape_matrix[i_x - 1][i_z + 1] == SHAPE_SPACE &&
		 shape_matrix[i_x    ][i_z - 1] == SHAPE_SPACE &&
		 shape_matrix[i_x    ][i_z + 1] == SHAPE_SPACE &&
		 shape_matrix[i_x + 1][i_z - 1] == SHAPE_SPACE &&
		 shape_matrix[i_x + 1][i_z    ] == SHAPE_SPACE &&
		 shape_matrix[i_x + 1][i_z + 1] == SHAPE_SPACE )
		{
		  shape_matrix[i_x][i_z] = SHAPE_SPACE ;
		  n_Clbond[i_x][i_z]     = 0 ;
		  std::cout << i_x << "\t" << i_z << "\n" ;
		}
	    }
	}
    }
}
