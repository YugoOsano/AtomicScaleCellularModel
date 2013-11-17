// ***  main_deposition.cc  *** 
//

#include <stdio.h>
#include <time.h>
#include <iostream.h>
#include <fstream.h>

#include "atom_struct.h"
#include "header_main.h"
#include "header_deposition.h"

#include "fileio_particle.h"

#include "particle.h"
#include "neutral_particle.h"
#include "etch_products.h"

#include "shape_trim.h"

int main ()
{
  time_t	computing_start,  computing_current;
  time(&computing_start);
  

  //== クラスの生成 ===

  class Particle_class      Cl_ion(CL_MASS, N_IADF_ARRAY ) ;
  class Neutral_class       Cl_neutral(CL_MASS ) ;
  class Shape_trim_class    Shape  ;		
  class Etch_product_class  SiCl2_neutral(CL_MASS) ;  

  if( ION_INJECT_FLAG == 3 ) //== IADF 読み込み
    Cl_ion.read_angular_df(IADF_FILE) ; 

  //== 一時的に使用する変数
  
  int    flag_out_of_domain , flag_deposition, flag_adsorption;
  int    i_exception_move ,i_exception_reflect ;
  int    i_x, i_z ;

  //== 計算条件の記録／出力
  ofstream  condition_file(CONDITION_FILE) ;
  ofstream  log_file(YIELD_FILE) ;
  cout  << "Number of ions to be injected: " << N_ION_INJECT 
	<< "\nReal time of evolution: "      << REAL_ELAPSED_TIME 
	<< "\nInterval for file output: "    << INTERVAL_FILE_OUTPUT 
	<< "\nYield beta: "    << YIELD_BETA
	<< "\n\narguments: neutral to ion flux ratio: " << NEUTRAL_ION_RATIO
	<< "\narguments: ion incident energy: "         << INCIDENT_ENERGY_ION
	<< endl << endl ;
  condition_file << "Number of ions to be injected: " << N_ION_INJECT 
		 << "\tReal time of evolution: "      << REAL_ELAPSED_TIME 
		 << "\tInterval for file output: "    
		 << INTERVAL_FILE_OUTPUT 
		 << "\nYield beta: "    << YIELD_BETA 
		 << "\n\narguments: neutral to ion flux ratio: " << NEUTRAL_ION_RATIO
		 << "\narguments: ion incident energy: "         << INCIDENT_ENERGY_ION<< endl ;
  
  //*****************************************
  //----    ループ開始         --------------
  //*****************************************

  int i_ion = 0 ;
  do
    {
      //== ループ回数、脱離Siの数の出力、記録
      if(i_ion % N_COUT_INTERVAL == 0)
	{
	  time(&computing_current);
	  
	  cout << "iteration: " << i_ion 
	       << "\tdesorbed Si: " << Shape.cntr_desorbed_Si 
	       << "\tcomputing: " 
	       << difftime(computing_current, computing_start) << endl ;

	  log_file << i_ion << "\t" 
		   << Shape.cntr_desorbed_Si << endl ;	  
	}
      //*****************************************
      // --        Cl ラジカルの入射
      //*****************************************
      
#ifdef _INJECT_CL_RADICAL_
      for (int i_neutral = 0; 
	   i_neutral < NEUTRAL_ION_RATIO; i_neutral++)
	{
	  //--- 入射 ---
	  Cl_neutral.inject_from_top( sqrt(2.0 * INCIDENT_ENERGY_NEUTRAL * 
					   Q_ELEMENTAL / Cl_neutral.mass) ) ; 
	  
	  i_exception_move = 0 ;
	  while(1) 
	    {
	      //--- 領域内部にあまり長くとどまっている場合、例外処理（粒子を消去）
	      i_exception_move++ ;
	      if( i_exception_move > N_CELL_Z * 200 )
		{
		  //cout << "i_ion: " << i_ion 
		  //   << "\ti_neutral: " << i_neutral << endl ; 
		  break ;
		}
	      /*cout << flag_out_of_domain << "\t"
		<< Cl_neutral.pos_v.x << "\t"
		   << Cl_neutral.pos_v.z 
		   << endl ;
	      */
	      flag_out_of_domain            // -- まず並進運動
		= Cl_neutral.move_trans(FREE_FLIGHT_PATH_Si) ;
	      if (flag_out_of_domain == FALSE) // -- 外に出た？
		break ;
	     
	      //--- セルが空白でない？ ---
	      if(Shape.put_shape_matrix(Cl_neutral.pos_v.x, 
					Cl_neutral.pos_v.z ,
					&i_x, &i_z ) != SHAPE_SPACE)
		{
		  flag_adsorption = 
		    Shape.settle_Cl_into_bareSi(i_x, i_z,  // 吸着した？
						ADSORPTION_PROBABILITY );
		  if(flag_adsorption == TRUE )
		    break;                          // 吸着したらループを出る
		  
		  // 乱反射 : 
		  // 無限ループをさけるために、カウンタを用意する。
		  // 粒子を仮に進ませてみて、(pos_v_togo)
		  // もし進んだ先が空白部分でなければ、元に戻って
		  //     やり直し(50回になったらもうやめる)
		  i_exception_reflect = 0;
		  Shape.get_surfacenormal(i_x, i_z) ;
		reflection:
		  if( i_exception_reflect > 50 )
		    break ;

		  Cl_neutral.random_reflection(FALSE, 0.0, 0.0 ) ; // 乱反射
		  // (TRUE,//<-cosine diffusion を導入すると左右非対称になる：当面は使用せず
		  // Shape.surfacenormal_x[i_x][i_z],
		  // Shape.surfacenormal_x[i_x][i_z]) ;
		
		  flag_out_of_domain 
		    = Cl_neutral.move_trans_togo(FREE_FLIGHT_PATH_Si) ;

		  if(Shape.put_shape_matrix(Cl_neutral.pos_v_togo.x, 
					    Cl_neutral.pos_v_togo.z ,
					    &i_x, &i_z ) != SHAPE_SPACE )
		    {
		      i_exception_reflect++ ;
		      goto reflection ;
		    }
		}
	    }
	}  // === END: Cl ラジカルの処理 
#endif


      //*************************************************
      //---   SiCl2  入射 （ION_SiCl2_RATIO回に１回）
      //*************************************************
      if( i_ion % ION_SiCl2_RATIO == 0 )
	{
	  SiCl2_neutral.inject_from_top( sqrt(2.0 * INCIDENT_ENERGY_NEUTRAL * 
					      Q_ELEMENTAL / SiCl2_neutral.mass) ) ; 
	  
	  i_exception_move = 0 ;
	  while(1) 
	    {
	      //--- 領域内部にあまり長くとどまっている場合、
	      //    例外処理（粒子を消去）
	      i_exception_move++ ;
	      if( i_exception_move > N_CELL_Z * 200 )
		break ;
	      
	      // -- 並進運動
	      flag_out_of_domain          
		= SiCl2_neutral.move_trans(FREE_FLIGHT_PATH_Si) ;
	      if (flag_out_of_domain == FALSE) // -- 外に出た？
		break ;
	      
	      // -- 粒子のいる場所のphaseを取得
	      SiCl2_neutral.get_position() ;
	      
	      // -- 固相であれば、一つ前のセル番号を取得 -> deposit
	      if(Shape.shape_matrix
		 [SiCl2_neutral.i_cell_x][SiCl2_neutral.i_cell_z] == SHAPE_Si )
		{
		  SiCl2_neutral.get_previous_cell() ; // セル番号取得
		  flag_deposition = 
		    Shape.deposit_Si( SiCl2_neutral.i_previous_x,     // deposit
				      SiCl2_neutral.i_previous_z,
				      2  ,SiCl2_DEPOSITION_PROBABILITY ) ;
		  
		  if(flag_deposition == TRUE )
		    break ;  // get out of the loop if deposited
		  
		  //乱反射 : 
		  // 無限ループをさけるために、カウンタを用意する。
		  // 粒子を仮に進ませてみて、(pos_v_togo)
		  // もし進んだ先が空白部分でなければ、元に戻って
		  //     やり直し(50回になったらもうやめる)
		  i_exception_reflect = 0;
		  Shape.get_surfacenormal(SiCl2_neutral.i_cell_x, 
					  SiCl2_neutral.i_cell_z) ;
		reflection_SiCl2:
		  if( i_exception_reflect > 50 )
		    break ;
		  
		  SiCl2_neutral.random_reflection (FALSE, 0.0, 0.0 ) ;
		  //(TRUE,  // 乱反射
		  //Shape.surfacenormal_x
		  //[SiCl2_neutral.i_cell_x][SiCl2_neutral.i_cell_z],
		  //Shape.surfacenormal_x
		  //[SiCl2_neutral.i_cell_x][SiCl2_neutral.i_cell_z]) ;
		  
		  flag_out_of_domain 
		    = SiCl2_neutral.move_trans_togo(FREE_FLIGHT_PATH_Si) ;
		  
		  if( Shape.put_shape_matrix(SiCl2_neutral.pos_v_togo.x, 
					     SiCl2_neutral.pos_v_togo.z ,
					     &i_x, &i_z ) != SHAPE_SPACE )
		    {
		      i_exception_reflect++ ;
		      goto reflection_SiCl2 ;
		    }
		  
		} 
	      
	    } //=== while loop 終了（１粒子のtrajectory）
	} // SiCl2 の処理終了	


	// == record shape
	if(i_ion % INTERVAL_FILE_DEPOSITION == 0 ||
	 i_ion == N_ION_INJECT ||
	 Shape.flag_get_bottom == TRUE ) 
	{
	  int tmp; 
	  char OUT1[50] ; 
		  
	  tmp = sprintf( OUT1, "%s%d.dat", SI_MATRIX_FILE, i_ion);
	  output_array_splot(OUT1 , 
			     Shape.shape_matrix ,
			     N_CELL_X, N_CELL_Z , 
			     CELL_SIZE ) ;
	}

	
	
	
	    i_ion++ ;
    }// for loop N;
  while( i_ion <= N_ION_INJECT &&
	 Shape.flag_get_bottom != TRUE ) ; // End of the do loop 


  return 0 ;
}

