// 

// -- simplified Monte Carlo simulation of 
//    low energy ion penetration into solids
// 
// p-109〜 を参照

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <time.h>

#include "atom_struct.h"
#include "header_main.h"

#include "random_number.h"
#include "fileio_particle.h"

#include "particle.h"
#include "neutral_particle.h"
#include "ion_particle.h"
#include "shape_trim.h"
#include "rotating_matrix.h"
#include "integration_angle.h"
#include "common_utils.h"

//------------
//-- reduced energy を出力する
// 入力：入射イオン、ターゲット固体、入射エネルギー(eV)

double put_reduced_energy(Atom_struct  Ion,   Atom_struct  Solid,
			  double energy_i ,   double *a_universal );

//--------------------------
int main(int argc, char *argv[])
{
  if(argc >= 2 ) // --- 引数はあるか？
    {
      NEUTRAL_ION_RATIO   = atof(argv[1]) ;
      INCIDENT_ENERGY_ION = atof(argv[2]) ;
      //YIELD_BETA=0.90537*sqrt(INCIDENT_ENERGY_ION)-4.219450374;
      //YIELD_BETA=0.36*(sqrt( INCIDENT_ENERGY_ION ) - sqrt(10.0));
      YIELD_BETA = 0.77*(sqrt( INCIDENT_ENERGY_ION ) - sqrt(20.0));
      std::cout << "YIELD_BETA : " << YIELD_BETA << std::endl ;
      //-- ファイル名
      /*sprintf( CL_BOND_FILE,   "%s%seV%sratio",
	       CL_BOND_FILE ,   argv[2], argv[1]);
      sprintf( SI_MATRIX_FILE, "%s%seV%sratio", 
	       SI_MATRIX_FILE , argv[2], argv[1]);
      sprintf( CONDITION_FILE, "%s%seV%sratio.dat", 
	       CONDITION_FILE , argv[2], argv[1]);
      */
      sprintf( YIELD_FILE,     "%s%seV%sratio.dat", 
	       YIELD_FILE ,     argv[2], argv[1]);
    }

  time_t	computing_start,  computing_current;
  time(&computing_start);

  sgenrand(RANDOM_NUMBER_SEED) ;
  srand( time(NULL) );
 
  //== クラスの生成 ===
  class Ion_class         Cl_ion(CL_MASS, 0, N_IADF_ARRAY ) ;
  class Neutral_class     Cl_neutral(CL_MASS, 0 ) ;
  class Shape_trim_class  Shape(FLAG_OPEN_SPACE)  ;

  if(FLAG_OPEN_SPACE == true)//== 片側 open spaceの場合
    {
      Cl_ion.set_flag_boundary(SPECULAR_REFLECT) ;
      Cl_neutral.set_flag_boundary(SPECULAR_REFLECT);//NO_PERIODIC) ;
    }
  //== flux count ==
  class Shape_counter_class Ion_counter ;
  class Shape_counter_class Neutral_counter ;
  if(FLAG_FLUX_COUNT == true)
    {
      YIELD_BETA = 0.0 ; //== エッチングは進展しないものとする
      Ion_counter.allocate_memory(); //== メモリ確保 ==
      Neutral_counter.allocate_memory();
    }

  if( ION_INJECT_FLAG == 3 ) //== IADF 読み込み
    Cl_ion.read_angular_df(IADF_FILE) ; 

  if(FLAG_INPUT_PROFILE == true) //== 形状読み込み ===
    Shape.input_profile(SI_INPUT_FILE, CL_INPUT_FILE ) ;

  //== 一時的に使用する変数
  
  bool   tmp_flag ;
  //  bool   flag_out_of_domain;
  int    i_x, i_z ;

  //int    i_counter ; // -- debug 用カウンタ

  //== 計算条件の記録／出力
  std::ofstream  condition_file(CONDITION_FILE) ;
  std::ofstream  log_file(YIELD_FILE) ;
  std::cout << "Pattern width(nm): " 
       << L_INTER_ATOMIC * (N_RIGHT_MASK - N_LEFT_MASK) * 1.0e+9
       << "\nNumber of ions to be injected: " << N_ION_INJECT 
       << "\nReal time of evolution: "        << REAL_ELAPSED_TIME 
       << "\nInterval for file output: "      << INTERVAL_FILE_OUTPUT 
       << "\nYield beta: "                    << YIELD_BETA
       << "\n\narguments: neutral to ion flux ratio: " << NEUTRAL_ION_RATIO
       << "\narguments: ion incident energy: "         << INCIDENT_ENERGY_ION
       << "\nflag of chemical etching: "      << FLAG_CHEMICAL_ETCH
       << "\nsubstrate temperature: "         << T_SUBSTRATE
       << "\nYield for chemical etching: "    << YIELD_CHEMICAL
       << "\nflag of forward scattering: "    << FLAG_FORWARD_SCATTER
       << "\nflag of mask erosion: "          << FLAG_MASK_EROSION
       << "\nslope angle of mask: "           << SLOPE_ANGLE_SIDEWALL
       << "\nflag of one sided open space: "  << FLAG_OPEN_SPACE
       << "\n\n" ;   
  condition_file 
    << "Pattern width(nm): " 
    << L_INTER_ATOMIC * (N_RIGHT_MASK - N_LEFT_MASK) * 1.0e+9
    << "\nNumber of ions to be injected: " << N_ION_INJECT 
    << "\tReal time of evolution: "      << REAL_ELAPSED_TIME 
    << "\tInterval for file output: "    
    << INTERVAL_FILE_OUTPUT 
    << "\nYield beta: "    << YIELD_BETA 
    << "\n\narguments: neutral to ion flux ratio: " << NEUTRAL_ION_RATIO
    << "\narguments: ion incident energy: "         << INCIDENT_ENERGY_ION
    << "\nflag of chemical etching: " << FLAG_CHEMICAL_ETCH
    << "\nsubstrate temperature: "    << T_SUBSTRATE
    << "\nYield for chemical etching: " << YIELD_CHEMICAL
    << "\nflag of forward scattering: "    << FLAG_FORWARD_SCATTER
    << "\nflag of mask erosion: "          << FLAG_MASK_EROSION
    << "\nslope angle of mask: "           << SLOPE_ANGLE_SIDEWALL
    << "\nflag of one sided open space: "  << FLAG_OPEN_SPACE
    << "\n\n" ;

#ifdef _DEPOSITION_PROGRAM_
#include "deposition/included_to_main1.cc"  
#endif


  //*****************************************
  //----    ループ開始         --------------
  //*****************************************
  
  int i_ion     = 0 ; 
  int i_neutral = 0 ; //int i_neutral_side = 0 ; 
  do
    //for(int i_ion = 1; i_ion <= N_ION_INJECT ; i_ion++)
    {
      //*****************************************
      //== ループ回数、脱離Siの数の出力、記録
      //*****************************************
      if(i_ion % N_COUT_INTERVAL == 0) 
	{
	  time(&computing_current);
	  
	  std::cout << "iteration: " << i_ion 
	       << "\tdesorbed Si: "  << Shape.cntr_desorbed_Si 
	       << "\tdeposited Si: " << Shape.cntr_deposited_Si
	       << "\tdesorbed mask: "<< Shape.cntr_desorbed_mask
	       << "\tcomputing: " 
	       << difftime(computing_current, computing_start) << std::endl ;
			 
	  log_file << i_ion << "\t" << Shape.cntr_desorbed_Si 
		   << "\t" << Shape.cntr_deposited_Si << std::endl ;
	}
      //*****************************************
      // --        Cl ラジカルの入射
      //*****************************************
    
#ifdef _INJECT_CL_RADICAL_
      while(double(i_ion) * NEUTRAL_ION_RATIO > double(i_neutral))
	{
	  i_neutral++ ;
	  Cl_neutral.all_process(&Shape, &Neutral_counter, false,
				 FLAG_FLUX_COUNT,
				 INCIDENT_ENERGY_NEUTRAL,
				 FLAG_CHEMICAL_ETCH, YIELD_CHEMICAL) ;
	  
	  //== 片側 open space の場合の処理 
	  // neutral/ion ratio と同様に、top/side ratioを元にして
	  // 横方向からの入射を行う
	  // if(FLAG_OPEN_SPACE == true)
	  //{
	  //  while(double(i_neutral) * SYSTEM_HEIGHT_Z / SYSTEM_WIDTH_X >
	  //    i_neutral_side)
	  //{
	  //  i_neutral_side++ ;
	  //  Cl_neutral.all_process(&Shape, &Neutral_counter, true,
	  //			 FLAG_FLUX_COUNT,
	  //		 INCIDENT_ENERGY_NEUTRAL,
	  //			 FLAG_CHEMICAL_ETCH, YIELD_CHEMICAL) ;
	  //	}
	  // }
	}  //  === END: Cl ラジカルの処理 
#endif
      
      //*****************************************
      // -- イオン入射：あるエネルギーで与えてやって、
      //    衝突を繰り返す。エネルギーがあるしきい値を下回ったら
      //　　ストップする。
      // もしくは跳ね返って z < 0 となった場合もストップ。
      //*****************************************
      // flag によって入射の条件を変える
      // 150eV -> m/s に変換
      
      // === 代表粒子数を 1/Ymax 倍するので、Ymax回繰り返す ===
      //i_counter = 0 ; //--debug
      for (int i_virtual_ion = 0; 
	   i_virtual_ion  < YIELD_MAX_INT; i_virtual_ion++)
	{
	  //== イオン入射 ==
	  if(ION_INJECT_FLAG == 2 )
	    Cl_ion.inject_from_center( sqrt(2.0 * INCIDENT_ENERGY_ION * 
					    Q_ELEMENTAL / Cl_ion.mass) ) ;
	  else if(ION_INJECT_FLAG == 3 )
	    Cl_ion.inject_iadf( sqrt(2.0 * INCIDENT_ENERGY_ION * 
				     Q_ELEMENTAL / Cl_ion.mass)) ;
	  else
	    Cl_ion.inject_from_top( sqrt(2.0 * INCIDENT_ENERGY_ION * 
					 Q_ELEMENTAL / Cl_ion.mass) ) ;
	  //std::cout << Cl_ion.pos_v.x << "\t" << Cl_ion.pos_v.y << "\t"
	  //  << Cl_ion.pos_v.z << std::endl ;
	  Cl_ion.set_flag_reflection(false);//== 前方散乱のflagをリセット 
	  
	  while(1) 
	    {
	      //-- 位置座標を求める（粒子の移動）
	      // -- free flight path(L) は一定と見なす
	      
	      Cl_ion.move_trans(FREE_FLIGHT_PATH_Si) ;
	      
	      if(Cl_ion.put_flag_inside() == false)
		break ; // 領域を飛び出した場合
	      
	      //==マスクに到達した場合も break ==
	      if(Shape.put_shape_matrix(Cl_ion.pos_v.x, 
					Cl_ion.pos_v.z ,
					&i_x, &i_z ) == HARD_MASK )
		{
		  if(FLAG_FLUX_COUNT == true)
		    Ion_counter.count(i_x, i_z ) ;
		  break ;
		}

	      if(FLAG_FORWARD_SCATTER == true)
		tmp_flag = Cl_ion.impact_scattering(&Shape, FLAG_MASK_EROSION);
	      else
		tmp_flag = Cl_ion.impact_on_surface(&Shape);

	      if(tmp_flag == true)
		break ;

	    } // <- end of while loop

	  Cl_ion.add_ctr_reflection(); //== 前方散乱のカウンタ加算 ==

	  // ==== エネルギーが小さくなったところがSi 原子であればそこに定着
	  if(Cl_ion.put_flag_inside() == true )
	    {
	      //***** 02Dec2003 **************
	      if(RAN0() < 1.0 / YIELD_MAX) //<- 吸着にかんしても仮想粒子を確率的に扱う
		Shape.adsorb_Cl(i_x, i_z) ;

	      // 座標出力（領域を飛び出たものではないか？）
	      //std::cout << Cl_ion.pos_v.x << "\t"
	      //<< Cl_ion.pos_v.y << "\t" << Cl_ion.pos_v.z << std::endl ;
	    }
	  
	  //*************************************************
	  //---   SiCl2  入射 （ION_SiCl2_RATIO回に１回）
	  //*************************************************
#ifdef _DEPOSITION_PROGRAM_
#include "deposition/included_to_main2.cc"  
#endif
	}// === for loop  END : イオン入射の処理

      //-- isolation check
      //Shape.remove_isolated_Si() ; <- コメントアウトしておく

      //****************************
      //--     ファイル出力
      //  -> 一つのセルにおけるClの数と酸素の数を同一のファイルで
      //     記録したい。-> 10進数で十の位を酸素の数、一の位をClの数とする
      //     -> プロットする際にmodをとってそれぞれの数を抽出する
      //****************************
      if(i_ion % INTERVAL_FILE_OUTPUT == 0 ||
	 i_ion == N_ION_INJECT ||
	 Shape.flag_get_bottom == true ) 
	//-- 底面に到達したらファイル出力してループ終了
	{
	  int tmp;   char OUT1[50] ; 
	  tmp = sprintf( OUT1, "%s%d.dat", CL_BOND_FILE, i_ion);
	  
	  //output_array_splot(OUT1 ,Shape.n_Clbond ,
	  //	     N_CELL_X, N_CELL_Z , CELL_SIZE ) ;
	  output_cellinfo(Shape.n_oxygen, Shape.n_Clbond ,
			  i_ion,  OUT1, 
			  N_CELL_X, N_CELL_Z , CELL_SIZE ) ;
			
	  
	  tmp = sprintf( OUT1, "%s%d.dat", SI_MATRIX_FILE, i_ion);
	  output_array_splot(OUT1 , Shape.shape_matrix ,
			     N_CELL_X, N_CELL_Z ) ;  
			    // CELL_SIZE ) ;

	  //Shape.output_Cl_density_depth
	  //  (CL_DENSITY_DEPTH_FILE ,i_ion );

	  if(FLAG_FLUX_COUNT == true && i_ion != 0)
	    {
	      tmp = sprintf( OUT1, "%s%d.dat", ION_COUNTER_FILE, i_ion);
	      output_array_splot(OUT1 , Ion_counter.ctr ,
				 N_CELL_X, N_CELL_Z ) ;  
	      tmp = sprintf( OUT1, "%s%d.dat", NEUTRAL_COUNTER_FILE, i_ion);
	      output_array_splot(OUT1 , Neutral_counter.ctr ,
				 N_CELL_X, N_CELL_Z ) ;  
	    }
	}

      i_ion++ ;
    }// do loop 終了
  while( i_ion <= N_ION_INJECT &&
	 Shape.flag_get_bottom != true ) ;

  
 
  // --- For check
  /*
  double a_u ; // universal radius
  double epsilon, t_loss, theta_c, psi ;

  for(int i = 1 ; i <= 200; i++)
    {
      epsilon = put_reduced_energy(Cl, Si, double(i)  , &a_u ) ;
      
      t_loss  = put_energy_loss(Cl, Si, double(i) ,
				2.0e-10 ,   // impact parameter 
				&theta_c, &psi ) ;
      std::cout << i << "\t" 
	   << "ε(reduced energy): " << epsilon   << "\t"
             << "Θ(theta_c): " << theta_c   << "\t" 
	   << "ψ(psi): " << psi       << "\t"
	   << "T(eV): " << t_loss << "\t" 
	   << std::endl ;
    }
   */
  return 0;

}

//------------
//-- reduced energy を出力する
// 入力：入射イオン、ターゲット固体、入射エネルギー(eV)

double put_reduced_energy(Atom_struct  Ion,   Atom_struct  Solid,
			  double energy_i ,   double *a_universal )
{
  // unversal screening length の決定
  *a_universal = 0.8853 * BOHR_RADIUS / (pow(Ion.atomic_n,   0.23) +
					 pow(Solid.atomic_n, 0.23) ) ;
  
  // 入射エネルギー(eV)を無次元化するためには、
  // 4πε a_u Ec / (Z1 Z2 e)  とする。
  // (Ec = E * M2 / (M1 + M2))
  return   
    4.0 * PI * E_PERMITTIVITY * (*a_universal) 
    * energy_i * (Solid.atomic_weight_solid / 
		  (Ion.atomic_weight + Solid.atomic_weight_solid ))
    / (Ion.atomic_n * Solid.atomic_n * Q_ELEMENTAL ) ;
  
}
