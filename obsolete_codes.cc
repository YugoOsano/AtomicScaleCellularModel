
//-- obsolete
 // sin^2(θ/2) の値をまず求める(4-13)
  /*
  double tmp_sin2 = 1.0 / (1.0 + 4.0 * reduced_energy * reduced_energy *
			   b_impact_para * b_impact_para) ;

  // scattering angle ψ (4-15)
  // ψ = arctan{sinθ / [cosθ + (M1/M2)]}

  // cosθ = 1 - 2 sin^2(θ/2), sinθ = √(1 - cos^2θ) を利用
  double tmp_cos = 1.0 - 2.0 * tmp_sin2 * tmp_sin2 ;

  *psi_scattering = atan( sqrt(1.0 - tmp_cos * tmp_cos) /
			  (tmp_cos + Ion.atomic_weight / Solid.atomic_weight_solid)) ;

  // T : (4-14)
  return 
    4.0 * Ion.atomic_weight * Solid.atomic_weight_solid * energy_i
    * tmp_sin2 / ((Ion.atomic_weight + Solid.atomic_weight_solid) *
		  (Ion.atomic_weight + Solid.atomic_weight_solid) ) ;
  */

 //-- 座標変換
/*
          v_x = tmp_v_x * cos(Cl_ion.pos_v.v_theta) * cos(Cl_ion.pos_v.v_psi) +
            tmp_v_y * (- sin(Cl_ion.pos_v.v_theta) ) +
              tmp_v_z * (- cos(Cl_ion.pos_v.v_theta)) * sin(Cl_ion.pos_v.v_psi) ;
      
          v_y = tmp_v_x * sin(Cl_ion.pos_v.v_theta) * cos(Cl_ion.pos_v.v_psi) +
            tmp_v_y * cos(Cl_ion.pos_v.v_theta)  +
              tmp_v_z * sin(Cl_ion.pos_v.v_theta) * sin(Cl_ion.pos_v.v_psi) ;
          
          v_z = tmp_v_x * sin(Cl_ion.pos_v.v_psi) +
            tmp_v_z * cos(Cl_ion.pos_v.v_psi) ;
  
  // x-z 平面でψ回転、x-y平面でθ回転を合成すると
          // ┌                               ┐
          // | cosθcosψ  -sinθ  -cosθsinψ |  
          // | sinθcosψ   cosθ   sinθsinψ |
          // |      sinψ      0         cosψ |
          // └                               ┘
          // ┌                               ┐
          // | cosψcosθ  -cosψsinθ  -sinψ |  
          // | sinθ             cosθ      0  |
          // | sinψcosθ  -sinψsinθ   cosψ |
          // └                               ┘
  
*/




//  tmp_x = Cl_neutral.pos_v.x ;// 前のステップの位置を記録  
//tmp_z = Cl_neutral.pos_v.z ;
// Cl_neutral.pos_v.x = tmp_x ;
// Cl_neutral.pos_v.z = tmp_z ;


//=== included_to_main2.cc 中のSiCl2の処理

 /*
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
	 [SiCl2_neutral.i_cell_x][SiCl2_neutral.i_cell_z] != SHAPE_SPACE )
	{
	  SiCl2_neutral.get_previous_cell() ; // セル番号取得
	  flag_deposition = 
	    Shape.deposit_Si( SiCl2_neutral.i_previous_x,     // deposit
			      SiCl2_neutral.i_previous_z,
				SiCl2_neutral.i_cell_x,    SiCl2_neutral.i_cell_z , 
				SiCl2_neutral.n_chlorination, 
				SiCl2_neutral.sticking_probability ) ;	  
	  if(flag_deposition == TRUE )
	  {
	    // cout << "SiCl2 deposited: \t" << SiCl2_neutral.i_previous_x 
	    // << "\t" << SiCl2_neutral.i_previous_z << endl ;
	    break ;  // get out of the loop if deposited
	  }
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
	  
	  SiCl2_neutral.random_reflection //(FALSE, 0.0, 0.0 ) ;
	  (TRUE,  // 乱反射
	  Shape.surfacenormal_x[SiCl2_neutral.i_cell_x][SiCl2_neutral.i_cell_z],
	  Shape.surfacenormal_z[SiCl2_neutral.i_cell_x][SiCl2_neutral.i_cell_z]) ;
	  
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
  */


//=================================
// monte_carlo.cc 中のエッチングの処理
// （Ion_class::ion_enhanced_etch 関数に置き換え）
//***** 02Dec2003 **************
//--- 飽和状態（吸着Cl が４つある）であれば、脱離させる
//	   --（かつ、今まで脱離反応にかかわっていなければ）
/*
if(Cl_ion.flag_contact == false)
{
  double tmp_beta ;
#ifndef  _INCIDENT_ANGLE_
  tmp_beta = yield_beta(tmp_energy);
#endif
#ifdef  _INCIDENT_ANGLE_
  tmp_beta = yield_beta(tmp_energy) * etch_yield_angle
    (Shape.get_incident_angle(i_x, i_z,
			      Cl_ion.pos_v.v_x,//<- pos_v.x としていたので、形状進展が変だった！29Aug2004 
			      Cl_ion.pos_v.v_z ),
     Shape.n_oxygen[i_x][i_z] ) ;
#endif
  if(RAN0() < tmp_beta / YIELD_MAX) // <- 仮想粒子を確率的に扱う
    {
      Shape.desorb_Si(i_x, i_z);
      Cl_ion.record_desorption() ; // flag for treatment of SiCl4
    }
  Cl_ion.flag_contact = true ;
  
  //== flux counting 
  if(FLAG_FLUX_COUNT == true)
    Ion_counter.count(i_x, i_z ) ;
}    
*/
//=================================
// monte_carlo.cc 中のイオンの表面衝突の処理
// （Ion_class::impact_on_surface 関数に置き換え）
  // -- そのセルに Si があれば、衝突
	      //  if( Shape.put_shape_matrix(Cl_ion.pos_v.x, 
	      //			 Cl_ion.pos_v.z ,
	      //			 &i_x, &i_z ) == SHAPE_Si )
	      /*
	      if( Shape.find_solid_nearest_neighbor
		  (Cl_ion.pos_v.x,    Cl_ion.pos_v.z ,
		   Cl_ion.pos_v.v_x,  Cl_ion.pos_v.v_z ,
		   &i_x, &i_z ) == SHAPE_Si )
		{
		  double p_impact_parameter ;
		  double tmp_energy ;
		  
		  tmp_energy = Cl_ion.mass * 
		    Cl_ion.pos_v.v_r * Cl_ion.pos_v.v_r 	      
		    /  (2.0 * Q_ELEMENTAL) ;  

		  //----  break if the energy is 
		  //      under the predetermined value
		  if(tmp_energy <= ENERGY_ION_STOPPING)
		    break ;
		  
		  Cl_ion.ion_enhanced_etch(&Shape,
					   tmp_energy, i_x, i_z ) ;
		  
		  // impact parameter は乱数Xで決定：
		  // p = (L/2)√X  where  L は原子間距離
		  p_impact_parameter = (L_INTER_ATOMIC / 2.0) * sqrt(RAN0());
		  
		  // -- collision
		  Cl_ion.collision_with_solid_atom( p_impact_parameter,
						    &tmp_energy );
		}
	      // == マスクであれば、break
	      else if(Shape.put_shape_matrix(Cl_ion.pos_v.x, 
					     Cl_ion.pos_v.z ,
					     &i_x, &i_z ) == HARD_MASK )
		{
		  flag_out_of_domain = false ;
		  break ;
		}
	      */
	      //=== 出力  =============
	      //cout << Cl_ion.pos_v.x << "\t"
	      //<< Cl_ion.pos_v.y << "\t"
	      //<< Cl_ion.pos_v.z << endl ;


//=================================
// monte_carlo.cc 中の中性粒子の処理
 //--- 入射 ---
	  /*
	  Cl_neutral.inject_from_top
	    ( sqrt(2.0 * INCIDENT_ENERGY_NEUTRAL * 
		   Q_ELEMENTAL / Cl_neutral.mass) ) ; 
	  
	  i_exception_move = 0 ;
	  while(1) 
	    {
	      //-- 領域内部にあまり長くとどまっている場合、
	      //   例外処理（粒子を消去）
	      i_exception_move++ ;
	      if( i_exception_move > N_CELL_Z * 200 )
		break ;
		
	      // -- まず並進運動 ----
	      flag_out_of_domain 
		= Cl_neutral.move_trans(FREE_FLIGHT_PATH_Si) ;
	      if (flag_out_of_domain == false) // -- 外に出た？
		break ;
	     
	      //--- セルが空白でない？ ---
	      //if(Shape.put_shape_matrix(Cl_neutral.pos_v.x, 
	      //			Cl_neutral.pos_v.z ,
	      //			&i_x, &i_z ) != SHAPE_SPACE)
	      if( Shape.find_solid_nearest_neighbor
		  (Cl_neutral.pos_v.x,    Cl_neutral.pos_v.z ,
		   Cl_neutral.pos_v.v_x,  Cl_neutral.pos_v.v_z ,
		   &i_x, &i_z, &i_x_particle, &i_z_particle ) != SHAPE_SPACE )
		{
		  //== flux counting 
		  if(FLAG_FLUX_COUNT == true )
		    Neutral_counter.count(i_x, i_z ) ;

		  //---- 化学的エッチング ----
		  // YIELD_CHEMICALの確率で脱離が起こる 
		  if( RAN0() < YIELD_CHEMICAL &&
		      FLAG_CHEMICAL_ETCH == true  )
		    {
		      Shape.desorb_Si_chemical_etch(i_x, i_z);
		    }
		  else //---- 吸着 -----
		    {
		      flag_adsorption = 
			Shape.settle_Cl_into_bareSi(i_x, i_z,  // 吸着した？
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
		  Shape.get_surfacenormal(i_x, i_z) ;
		reflection:
		  if( i_exception_reflect > 50 )
		    break ;

		  Cl_neutral.random_reflection //(FALSE,0.0,0.0);//乱反射
		   (true,//<-cosine diffusion を導入すると左右非対称になる<- this bug was fixed (29Aug2004)
		   Shape.surfacenormal_x[i_x][i_z],
		   Shape.surfacenormal_z[i_x][i_z]) ;
		
		  flag_out_of_domain 
		    = Cl_neutral.move_trans_togo(FREE_FLIGHT_PATH_Si) ;

		  if(Shape.put_shape_matrix(Cl_neutral.pos_v_togo.x, 
					    Cl_neutral.pos_v_togo.z ,
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
	  */

const double GATE_WIDTH   = 100.0 * 1.0e-9 ;//== ゲート幅(nm) ==
const double MARGIN_WIDTH =  50.0 * 1.0e-9 ;//== 領域左右の余白(nm)
=     int(floor(GATE_WIDTH / CELL_SIZE) + 0.5)  
  + 2*int(floor(MARGIN_WIDTH/CELL_SIZE) + 0.5) ; //floor関数による四捨五入

const int    N_LEFT_MASK  = int(floor(MARGIN_WIDTH/CELL_SIZE) + 0.5);//60;//35;//70;// マスク左壁面の位置
const int    N_RIGHT_MASK = N_CELL_X - int(floor(MARGIN_WIDTH/CELL_SIZE) + 0.5);//60;//35;//70;// マスク右壁面の位置
