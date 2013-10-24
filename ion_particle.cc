
#include <iostream>
#include "particle.h"
#include "random_number.h"
#include "iedf_to_distribution.h"

#include "ion_particle.h"

//==  入射エネルギーに対して etch yield を返す ==
//    yield の係数β = C(√E - √Eth): C = 0.77, Eth = 20.0 
//    *: 0.77 * sqrt(20.0) = 3.44354468534968

inline double yield_beta(double incident_energy)
{
  double  tmp ;
  tmp = 0.77 * sqrt(incident_energy) - 3.44354468534968 ;

  if(tmp > 0.0)
    return tmp ;
  else 
    return 0.0 ;
}
//== イオンが反射した場合のエッチング収率は Er < Eth であれば C(√Ei - √Eth),
//                                      Er > Eth であれば C(√Ei - √Er)とする
inline double yield_beta(double incident_energy, double reflected_energy)
{
  double  tmp ;
  if(reflected_energy > 20.0)
    tmp = 0.77 *(sqrt(incident_energy) - sqrt(reflected_energy));
  else
    tmp = 0.77 * sqrt(incident_energy) - 3.44354468534968 ;
  
  if(tmp > 0.0)
    return tmp ;
  else 
    return 0.0 ;
}

//======================================

Ion_class::Ion_class() 
{
}

Ion_class::Ion_class(double mass_input,  int flag_boundary_input,
		     int n_adf_input ) 
  : Particle_class(mass_input, flag_boundary_input) 
{
  n_adf_array = n_adf_input ;

  angular_df           =  new double[n_adf_array] ;
  angle_array_degree   =  new double[n_adf_array] ;

  flag_contact    = false ;
  flag_desorb_Si  = false ;
  flag_reflection = false ;
  ctr_reflection  = 0 ;
  velocity_file.open(NULL);//== 散乱後の速度を記録 ==

  i_x_etch = 0; i_z_etch = 0 ;
  energy_etch = 0.0 ;
}
Ion_class::~Ion_class() 
{
  delete [] angular_df  ;
  delete [] angle_array_degree ;
  velocity_file.close();
}

//--==================================================
void Ion_class::inject_iadf(double v )
{
  flag_contact   = false ; 
  flag_desorb_Si = false ;

  // randomized_distribution   (long int 乱数の種, int 作りだす乱数の数,
  //   double[] エネルギー  ,double[] 分布,
  //   -> 出力  double[] (元の乱数), double[] 得られた乱数) 
  
  // -- 動的確保！
  double *dummy_array_angle ;

  double *randomized_angle_array ;
  
  dummy_array_angle      = new double[1] ;
  randomized_angle_array = new double[1] ;

  randomized_distribution (RANDOM_NUMBER_SEED,
			   1 ,
			   angle_array_degree , 0.0,
			   angular_df, 
			   dummy_array_angle,
			   randomized_angle_array) ;

  //cout << "randomized_angle_array[0]: " << randomized_angle_array[0] << "\n";

  // X座標は乱数によって決定する
  pos_v.x = RAN0() * SYSTEM_WIDTH_X ; 
  pos_v.y = 0.0 ;
  pos_v.z = 0.0 ;
  
  // --- 入射角はこの時点で決定している(randomized_angle_array)
  // が、y 軸周りの回転角は乱数で作る 
  // 入射角ψ、回転角θとすると z = v cosψ, x = v sinψcosθ, y = v sinψsinθ
      
  // randomized_energy は eV の単位なので Q_ELEMENTAL をかけてVolt 単位にする

  pos_v.v_r     =  v ; // <--- Absolute velocity

  pos_v.v_theta =  2.0 * PI * RAN0() ; 
     
  pos_v.v_psi   = 
    PI * randomized_angle_array[0] / 180.0 ; // degree -> radへ換算 
  
  //----
  pos_v.v_x =  //0.0 ;
    pos_v.v_r * sin(pos_v.v_psi) * cos(pos_v.v_theta) ;
  pos_v.v_y = 
    pos_v.v_r * sin(pos_v.v_psi) * sin(pos_v.v_theta) ;

  pos_v.v_z = 
    pos_v.v_r * cos(pos_v.v_psi) ;
 
   //cout << pos_v.v_x << "\t" << pos_v.v_y << "\t"
     //<< pos_v.v_z << "\n" ;

  delete [] dummy_array_angle ;
  delete [] randomized_angle_array ;

  //== エッチング用変数初期化 ==
  i_x_etch = 0; i_z_etch = 0 ;   energy_etch = 0.0 ;
  pos_v_etch = pos_v ; 

}// --- End of inject_ion


//==  
void Ion_class::record_desorption(Particle_location_velocity_struct pos_v_recorded,
				  int input_n_oxygen) 
{
  flag_desorb_Si  = true ;
  n_oxy_desorb_Si = input_n_oxygen ;

  position_at_desorption  = pos_v_recorded ;
}
	
//== イオンの衝突によるエッチングの処理
//bug fix:<- pos_v.x としていたので、形状進展が変だった:29Aug2004 
//bug fix(*): 28/Nov/2005
//  ローカルな入射角の決定を「現時点での速度(pos_v)」で行っていた
void Ion_class::
ion_enhanced_etch(class  Shape_trim_class *Shape_trim,
		  double incident_energy, double reflected_energy,
		  Particle_location_velocity_struct pos_v_recorded,
		  int i_x, int i_z ) 
{
  double tmp_beta = 0.0;
  if(FLAG_INCIDENT_ANGLE == true)
    {
      //== (3rd paper)酸素が吸着することによる形状不安定を防ぐため、
      //   エッチング収率を分散させる。nearest-neighborのセルに含まれる酸素の数を合計し、
      //   6 で割ってcoverageΘとする。
      //   収率は Y = Θ Y(SiO2) + (1 - Θ)Y(Si)
      
      if(FLAG_DISPERSE_OXIDATION == true)
	{
	  int ctr_oxygen = 0;
	  if(i_x > 0  &&  i_x < N_CELL_X - 1 &&
	     i_z > 0  &&  i_z < N_CELL_Z - 1 )
	    ctr_oxygen = Shape_trim->n_oxygen[i_x][i_z] 
	      +          Shape_trim->n_oxygen[i_x + 1][i_z] + Shape_trim->n_oxygen[i_x - 1][i_z] 
	      +          Shape_trim->n_oxygen[i_x][i_z + 1] + Shape_trim->n_oxygen[i_x][i_z - 1] ;
	  else
	    ctr_oxygen = Shape_trim->n_oxygen[i_x][i_z] ;

	  tmp_beta = yield_beta(incident_energy, reflected_energy) * etch_yield_disperse_oxidation
	    (Shape_trim->get_incident_angle
	     (i_x, i_z,	    pos_v_recorded.v_x, pos_v_recorded.v_z ), //<-*
	     ctr_oxygen );
	}
      else //== (2nd paper) ==
	{
	   tmp_beta = yield_beta(incident_energy, reflected_energy) * etch_yield_angle
	    (Shape_trim->get_incident_angle
	     (i_x, i_z,	    pos_v_recorded.v_x, pos_v_recorded.v_z ), //<-*
	     Shape_trim->n_oxygen[i_x][i_z] ) ;
	}
    }
  else
    {
      tmp_beta = yield_beta(incident_energy, reflected_energy);
    }

  if(RAN0() < tmp_beta / YIELD_MAX) // <- 仮想粒子を確率的に扱う
    {
      Shape_trim->desorb_Si(i_x, i_z);
      record_desorption(pos_v_recorded, Shape_trim->n_oxygen[i_x][i_z]) ; 
    }
  //== flux counting 
  //    if(FLAG_FLUX_COUNT == true)
  //	Ion_counter.count(i_x, i_z ) ;
}

//==
void Ion_class::
hardmask_sputter(class  Shape_trim_class *Shape_trim,
		 double incident_energy, double reflected_energy,
		 Particle_location_velocity_struct pos_v_recorded,
		 int i_x, int i_z ) 
{
  double tmp_beta ;
 
  tmp_beta = SELECTIVITY_HARDMASK * 
    yield_beta(incident_energy, reflected_energy) * 
    etch_yield_angle(Shape_trim->get_incident_angle
		     (i_x, i_z, pos_v_recorded.v_x, pos_v_recorded.v_z ), 2 ) ; 
  //etch_yield_angle関数の最後の引数を2 とすることで酸化膜の角度依存性が得られる
 
  if(RAN0() < tmp_beta / YIELD_MAX) // <- 仮想粒子を確率的に扱う
    {
      Shape_trim->desorb_mask(i_x, i_z);
      record_desorption(pos_v_recorded, Shape_trim->n_oxygen[i_x][i_z]) ; 
    }
}

//==   パターン表面に到達したかどうかの判定 
bool Ion_class::
impact_on_surface(class  Shape_trim_class *Shape_trim)
{
  int i_x, i_z , i_x_particle, i_z_particle ;
  // -- そのセルに Si があれば、衝突
  //  if( Shape.put_shape_matrix(Cl_ion.pos_v.x, 
  //			 Cl_ion.pos_v.z ,
  //			 &i_x, &i_z ) == SHAPE_Si )
  if( Shape_trim->find_solid_nearest_neighbor
      (pos_v.x,    pos_v.z ,
       pos_v.v_x,  pos_v.v_z ,
       &i_x, &i_z, &i_x_particle, &i_z_particle ) == SHAPE_Si )
    {
      double p_impact_parameter ;
      double incident_energy ;
		  
      incident_energy = mass * pos_v.v_r * pos_v.v_r /(2.0 * Q_ELEMENTAL) ; 

      //----  break if the energy is 
      //      under the predetermined value
      if(incident_energy <= ENERGY_ION_STOPPING)
	return true ;
      
      //== 最表面でのみエッチング処理を行うため、接触しているか
      //   フラグで判定
      if(flag_contact == false)
	{
	  ion_enhanced_etch(Shape_trim, incident_energy , 0.0,  pos_v, i_x, i_z ) ;
	  flag_contact = true ; //== 固体内部でのエッチングはなし ==
	}
		  
      // impact parameter は乱数Xで決定：
      // p = (L/2)√X  where  L は原子間距離
      p_impact_parameter = (L_INTER_ATOMIC / 2.0) * sqrt(RAN0());
		  
      // -- collision
      collision_with_solid_atom( p_impact_parameter );
    }
  // == マスクであれば、break
  else if(Shape_trim->put_shape_matrix(pos_v.x, pos_v.z ,
				       &i_x, &i_z ) == HARD_MASK )
    {
      flag_inside = false ;
      return true ;
    }
  return false ;
}

//==   上と同じく衝突の処理：前方散乱のも含める
bool Ion_class::
impact_scattering(class  Shape_trim_class *Shape_trim,
		  bool  flag_mask_erosion) 
{
  int i_x, i_z, i_x_particle, i_z_particle ;
  double p_impact_parameter ;
  double incident_energy ;
  incident_energy = mass * pos_v.v_r * pos_v.v_r /(2.0 * Q_ELEMENTAL) ; 

  // == マスクであれば、break
  if(Shape_trim->put_shape_matrix(pos_v.x, pos_v.z ,
				  &i_x, &i_z ) == HARD_MASK )
    {
      flag_inside = false ;
      return true ;
    }  
  // -- そのセルに Si があれば、衝突
  else if(Shape_trim->put_shape_matrix(pos_v.x, pos_v.z ,
				       &i_x, &i_z ) == SHAPE_Si )
    {
      //エネルギーがある値を下回っていたら stop
      if(incident_energy <= ENERGY_ION_STOPPING)
	{
	  //== stopした時点で事後的にエッチングの処理を行う
	  ion_enhanced_etch(Shape_trim, energy_etch, 0.0, pos_v_etch,
			    i_x_etch , i_z_etch ) ;
	  return true ;
	}
      if(flag_contact == false)
	{
	  record_etch_position(i_x, i_z, incident_energy);
	  flag_contact = true ; 
	}
      
      // impact parameter は乱数Xで決定：
      // p = (L/2)√X  where  L は原子間距離
      p_impact_parameter = (L_INTER_ATOMIC / 2.0) * sqrt(RAN0());
      
      // -- collision --
      //collision_with_solid_atom( p_impact_parameter );
      collision_accurate3D( (i_x + 0.5) * L_INTER_ATOMIC,
			    0.5         * L_INTER_ATOMIC,
			    (i_z + 0.5) * L_INTER_ATOMIC,
			    &incident_energy );
    }
  // -- 隣接セルが Si or hardmask の場合
  else if( Shape_trim->find_solid_nearest_neighbor
	   (pos_v.x,    pos_v.z ,
	    pos_v.v_x,  pos_v.v_z ,
	    &i_x, &i_z, &i_x_particle, &i_z_particle) == SHAPE_Si )
    {
      if(flag_contact == false)
	{
	  record_etch_position(i_x, i_z, incident_energy);
	  flag_contact = true ; 
	}
      // -- collision （周期境界条件を考慮）
      if(i_x_particle == 0 && i_x == N_CELL_X - 1)
	collision_accurate3D( (- 0.5)     * L_INTER_ATOMIC,
			      0.5         * L_INTER_ATOMIC,
			      (i_z + 0.5) * L_INTER_ATOMIC,
	 		    &incident_energy );
      else if(i_x_particle == N_CELL_X - 1 && i_x == 0)
	collision_accurate3D( (N_CELL_X + 0.5) * L_INTER_ATOMIC,
			      0.5         * L_INTER_ATOMIC,
			      (i_z + 0.5) * L_INTER_ATOMIC,
			    &incident_energy );
      else
	collision_accurate3D( (i_x + 0.5) * L_INTER_ATOMIC,
			      0.5         * L_INTER_ATOMIC,
			      (i_z + 0.5) * L_INTER_ATOMIC,
			      &incident_energy );

      //if(i_z < 360 && i_z > 11)
      //cout << i_x << "\t" << i_z << "\t" 
      //     << pos_v.v_x << "\t" << pos_v.v_z << "\n" ;

    }
  else if(Shape_trim->find_solid_nearest_neighbor
	  (pos_v.x,    pos_v.z ,
	   pos_v.v_x,  pos_v.v_z ,
	   &i_x, &i_z, &i_x_particle, &i_z_particle) == HARD_MASK )
    {
      if(flag_contact == false)
	{
	  if(flag_mask_erosion == true)
	    hardmask_sputter(Shape_trim, incident_energy, 
			     0.0, pos_v, i_x, i_z);//<-added
	  flag_contact = true ; 
	}
       // -- collision （周期境界条件を考慮）
       if(i_x_particle == 0 && i_x == N_CELL_X - 1)
	 collision_accurate3D( (- 0.5)     * L_INTER_ATOMIC,
			       0.5         * L_INTER_ATOMIC,
			       (i_z + 0.5) * L_INTER_ATOMIC,
			       &incident_energy );
       else if(i_x_particle == N_CELL_X - 1 && i_x == 0)
	 collision_accurate3D( (N_CELL_X + 0.5) * L_INTER_ATOMIC,
			       0.5         * L_INTER_ATOMIC,
			       (i_z + 0.5) * L_INTER_ATOMIC,
			        &incident_energy );
       else
	 collision_accurate3D( (i_x + 0.5) * L_INTER_ATOMIC,
			       0.5         * L_INTER_ATOMIC,
			       (i_z + 0.5) * L_INTER_ATOMIC,
			       &incident_energy );
    }
  else if( Shape_trim->find_solid_nearest_neighbor
	   (pos_v.x,    pos_v.z ,
	    pos_v.v_x,  pos_v.v_z ,
	    &i_x, &i_z, &i_x_particle, &i_z_particle) == SHAPE_SPACE )
    {
      if(flag_contact    == true)//== 固体部分に接触状態 -> detachしたことを意味する
	{
	  flag_reflection = true ; //== 前方散乱が起こったことを記録
	  //velocity_file << (pos_v.v_psi * 180 / PI - 90.0 ) << "\t" 
	  //	<< incident_energy / 50.0 << "\n" ;

	  //== イオンがreflectionを起こした場合も事後的にエッチング処理する
	  //   収率はエネルギー変化分
	  //  （エッチング衝突時−現在のEi）のみによるものとする
	  //ion_enhanced_etch(Shape_trim, energy_etch, 0.0, pos_v_etch,
	  //	    i_x_etch , i_z_etch ) ;
	  ion_enhanced_etch(Shape_trim, energy_etch, incident_energy, pos_v_etch,
			    i_x_etch , i_z_etch ) ;
	}
      flag_contact    = false ;
    }
  
  return false ;
}
