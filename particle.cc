
#include <iostream>
#include <math.h>
#include "particle.h"
#include "random_number.h"
#include "integration_angle.h"
#include "rotating_matrix.h"
#include "iedf_to_distribution.h"
#include "fileio_particle.h"
#include "common_utils.h"

Particle_class::Particle_class()
{
  flag_contact  = false ; 
}
Particle_class::Particle_class(double mass_input, int flag_boundary_input ) 
{
  mass          = mass_input  ;
  flag_boundary = flag_boundary_input;
  n_injection   = 0 ;
} 

Particle_class::~Particle_class()
{
  delete [] angular_df  ;
  delete [] angle_array_degree ;
}

//== 角度分布をファイルから取得する =====
//   下記の理由により、flux_model から修正が必要
void Particle_class::read_angular_df(char filename[])
{
  //== microstructure と同じく、
  // distribution に変換する必要がある。
  iedf_to_distribution(filename , n_adf_array ,
		       angle_array_degree , 
		       angular_df  	) ;
  /*
    input_array(filename ,angle_array_degree , angular_df,
    n_adf_array ) ;
    // (正規化は必要!)
    double tmp_summation = 0.0 ;
    for(int i = 0; i < n_adf_array; i++)
    {
    tmp_summation += angular_df[i] ;
    }
    for(int i = 0; i < n_adf_array; i++)
    {
    angular_df[i] = angular_df[i] / tmp_summation ;
    //cout << "angular_df[" << i << "]" << angular_df[i] << "\n" ;
    }
  */
  //cout << "n_adf_array: " << n_adf_array << endl ;
}

void Particle_class::inject_from_center( double v_z )
{
  pos_v.x = SYSTEM_WIDTH_X / 2.0  ; 
  pos_v.y = 0.0 ;
  pos_v.z = 0.0 ;

  pos_v.v_x = 0.0 ;
  pos_v.v_y = 0.0 ;
  pos_v.v_z = v_z ;

  pos_v.v_r     = v_z ;
  pos_v.v_theta = 0.0 ;
  pos_v.v_psi   = 0.0 ;

  flag_contact  = false ; 
}


void Particle_class::inject_from_top(double v_z)
{
  double tmp_random = RAN0() ;
  pos_v.x = tmp_random * SYSTEM_WIDTH_X ; 

  pos_v.y = 0.0 ;
  pos_v.z = CELL_SIZE * RAN0() ; 
  // <- 0.0 にすると、あるセルを「飛びこえて」しまう

  pos_v.v_x = 0.0 ;
  pos_v.v_y = 0.0 ;
  pos_v.v_z = v_z ;

  pos_v.v_r     = v_z ;
  pos_v.v_theta = 0.0 ;
  pos_v.v_psi   = 0.0 ;

  flag_contact  = false ; 
}

//-- 粒子入射：斜め方向に初期速度を与える
void Particle_class::
inject_oblique(double x, double psi, double v ) 
{
  //  double tmp_random = RAN0() ;
  pos_v.x = x ; //0.5 * SYSTEM_WIDTH_X ; 
  pos_v.y = RAN0() * CELL_SIZE ;//0.0 ;
  pos_v.z = 0.0 ;

  pos_v.v_r     = v   ;
  pos_v.v_theta = 2.0 * PI * RAN0() ; ;
  pos_v.v_psi   = psi ;
  
  pos_v.v_x = pos_v.v_r * sin(pos_v.v_psi) * cos(pos_v.v_theta) ;
  pos_v.v_y = pos_v.v_r * sin(pos_v.v_psi) * sin(pos_v.v_theta) ;
  pos_v.v_z = pos_v.v_r * cos(pos_v.v_psi) ;

  flag_contact  = false ; 
}


//-------------------------------------
void Particle_class::get_position()
{
  if(pos_v.x >= 0.0 && pos_v.x < SYSTEM_WIDTH_X  &&
     pos_v.z >= 0.0 && pos_v.z < SYSTEM_HEIGHT_Z  )
    {
      //cout << shape_matrix[ int( x / CELL_SIZE ) ][ int( z / CELL_SIZE ) ]
      //   << endl ;

      i_cell_x = int( pos_v.x / CELL_SIZE ) ;
      i_cell_z = int( pos_v.z / CELL_SIZE ) ;

      dx_cell  = pos_v.x - i_cell_x * CELL_SIZE ;
      dz_cell  = pos_v.z - i_cell_z * CELL_SIZE ;
    }
  else
    {
      i_cell_x = -1 ;
      i_cell_z = -1 ;
      dx_cell  =  0.0 ;
      dz_cell  =  0.0 ;
    }
}

//-------------------------------------------
void  Particle_class::move_trans(double  l )
{
  pos_v.x += l * pos_v.v_x
    /  sqrt(pos_v.v_x * pos_v.v_x + 
	    pos_v.v_y * pos_v.v_y + 
	    pos_v.v_z * pos_v.v_z) ;
          
  pos_v.y += l * pos_v.v_y
    /  sqrt(pos_v.v_x * pos_v.v_x + 
	    pos_v.v_y * pos_v.v_y + 
	    pos_v.v_z * pos_v.v_z) ;
  
  pos_v.z += l * pos_v.v_z
    /  sqrt(pos_v.v_x * pos_v.v_x + 
	    pos_v.v_y * pos_v.v_y + 
	    pos_v.v_z * pos_v.v_z) ;
  //--------------------
  //== 境界条件の処理 ==
  if(flag_boundary == NO_PERIODIC)//--境界条件なし
    {
      if(pos_v.x < 0 || pos_v.x > SYSTEM_WIDTH_X)
	flag_inside = false ;
    }
  else if(flag_boundary == SPECULAR_REFLECT)//--境界で鏡面反射
    {
      if(pos_v.x < 0) // 左から出た場合
	{
	  pos_v.x   = - pos_v.x ;
	  pos_v.v_x = - pos_v.v_x ;
	}
      if(pos_v.x > SYSTEM_WIDTH_X ) // 右から出た場合
	{
	  pos_v.x   = 2.0 * SYSTEM_WIDTH_X - pos_v.x ;
	  pos_v.v_x = - pos_v.v_x ;
	}
    }
  else // -- 周期境界条件
    {
      if(pos_v.x < 0) // 左から出た場合
	pos_v.x = SYSTEM_WIDTH_X - fabs(pos_v.x) ;
      
      if(pos_v.x > SYSTEM_WIDTH_X ) // 右から出た場合
	pos_v.x = pos_v.x - SYSTEM_WIDTH_X ;
    }
  //-- 奥行き方向
  if(pos_v.y < 0)
    pos_v.y = L_INTER_ATOMIC + pos_v.y ;

  if(pos_v.y > L_INTER_ATOMIC)
    pos_v.y = pos_v.y - L_INTER_ATOMIC ;
    
  // -- 上下に出たらflagをfalseにする
  if(pos_v.z < 0 || 
     pos_v.z > SYSTEM_HEIGHT_Z )
    flag_inside = false ;
  else
    flag_inside = true ;
}

// -- pos_v_togo を求める
void  Particle_class::move_trans_togo(double  l )
{
  pos_v_togo.x = pos_v.x + l * pos_v.v_x
    /  sqrt(pos_v.v_x * pos_v.v_x + 
	    pos_v.v_y * pos_v.v_y + 
	    pos_v.v_z * pos_v.v_z) ;
  
  pos_v_togo.y = pos_v.y + l * pos_v.v_y
    /  sqrt(pos_v.v_x * pos_v.v_x + 
	    pos_v.v_y * pos_v.v_y + 
	    pos_v.v_z * pos_v.v_z) ;
  
  pos_v_togo.z = pos_v.z + l * pos_v.v_z
    /  sqrt(pos_v.v_x * pos_v.v_x + 
	    pos_v.v_y * pos_v.v_y + 
	    pos_v.v_z * pos_v.v_z) ;
  //------------------
  // -- 周期境界条件を考慮
  if(pos_v_togo.x < 0) // 左から出た場合
    pos_v_togo.x = SYSTEM_WIDTH_X - fabs(pos_v_togo.x) ;

  if(pos_v_togo.x > SYSTEM_WIDTH_X ) // 右から出た場合
    pos_v_togo.x = pos_v_togo.x - SYSTEM_WIDTH_X ;
    
  // -- 上に出たらFALSEを返す
  if(pos_v_togo.z < 0 || 
     pos_v_togo.z > SYSTEM_HEIGHT_Z )
    flag_inside_togo = false ;
  else
    flag_inside_togo =  true ;
}

//-- 速度の変更 (x, z 方向成分のみ）
void Particle_class::random_reflection(bool flag_cosine_dist ,
			 double normal_x, double normal_z)
{
  double theta ;
  if (flag_cosine_dist  == true )
    {
      // θn + α ：θnは法線方向、αは cosine distribution
      // normal: 空間 -> 固体内部の方向 であることに注意
      theta = rotating_angle(- normal_x, - normal_z) 
	+ asin(2.0 * RAN0() - 1.0) ;
      /*	cout << rotating_angle(- normal_x, - normal_z)   << "\t" 
		<<  asin(2.0 * RAN0() - 1.0)  << "\t" 
		<<  "reflection angle: " << theta << endl ;
      */
    } 
  else
    {
      theta = RAN0() * 2.0 * PI ;
    }
  double v     = sqrt(pos_v.v_x * pos_v.v_x + pos_v.v_z * pos_v.v_z ) ;
  
  pos_v.v_x = v * cos(theta) ;
  pos_v.v_z = v * sin(theta) ;
  
  //-- 速度：デカルト座標→球面座標
  
  pos_v.v_r = sqrt(pos_v.v_x * pos_v.v_x + 
		   pos_v.v_y * pos_v.v_y + 
		   pos_v.v_z * pos_v.v_z) ;
  pos_v.v_theta = atan(pos_v.v_y / pos_v.v_x) ;
  pos_v.v_psi   = acos(pos_v.v_z / pos_v.v_r) ;

}
//-----------------------------------------------------
double Particle_class::put_energy_loss(Atom_struct  Ion,   
				       Atom_struct  Solid,
				       double  energy_i ,
				       double  b_impact_para,
				       double  *theta_c ,    
				       double  *psi_scattering )
{
 
  //double theta_c ; // 重心座標における跳ね返り角

  *theta_c = theta_integrate_NC(b_impact_para,
				energy_i ) ;

  // scattering angle ψ (4-15)
  // ψ = arctan{sinθ / [cosθ + (M1/M2)]}
  
  *psi_scattering = atan( sin(*theta_c) / 
			  (cos(*theta_c) 
			   + Ion.atomic_weight / Solid.atomic_weight_solid )) ;
  
  // T : (4-14)
  return 
    4.0 * Ion.atomic_weight * Solid.atomic_weight_solid * energy_i
    * sin(*theta_c / 2.0) * sin(*theta_c / 2.0) 
    / ((Ion.atomic_weight + Solid.atomic_weight_solid) *
       (Ion.atomic_weight + Solid.atomic_weight_solid) ) ;
}

//-----------------------------------------------------
void  Particle_class::
collision_with_solid_atom(double p_impact_parameter )
{
  double psi_c, psi ,theta  ;
  double tmp_energy_loss ;

  // psi_c : 重心中心座標におけるscattering angle
  // psi   : laboratoryの座標におけるscattering angle
  // theta : scattering の回転角（ランダムに決定）

  // -- 衝突前の進行方向軸からの相対的な速度
  double tmp_v_x, tmp_v_y, tmp_v_z , tmp_v_r ;  

  // 変換のための直交行列
  double a11, a12, a13 ;
  double a21, a22, a23 ;
  double a31, a32, a33 ;
  bool    sign_flag ;  //  sqrtの正負

  // -- 衝突前の進行方向軸によって変換を施した速度
  //  double v_x, v_y, v_z ;

  energy_k = mass * pos_v.v_r * pos_v.v_r 	      
    /  (2.0 * Q_ELEMENTAL) ;  

  // -- collision
  tmp_energy_loss = put_energy_loss(Cl, Si, energy_k ,
				    p_impact_parameter,
				    &psi_c, &psi ) ;
          
  // エネルギー消失
  energy_k -= tmp_energy_loss ;
          
  //-- 速度変化：衝突前の速度方向の軸回りの
  //　　　　　　衝突角をランダムに決定(2 * PI * RAN0() )
  // これをもとに速度の x, y, z 成分を求めて、さらに
  // 衝突前の速度ベクトルを元に座標変換を行う。
  
  //  座標変換の変換行列の決定
  //  (sqrtの正負-> 乱数で、0.5以上なら正、0.5未満は負
  if(RAN0() >= 0.5)
    sign_flag = true  ;
  else
    sign_flag = false ;
  
  rotating_matrix(pos_v.v_theta, pos_v.v_psi,
		  sign_flag , 
		  &a11, &a12, &a13 ,
		  &a21, &a22, &a23 ,
		  &a31, &a32, &a33 ) ;
                            
  theta = 2.0 * PI * RAN0() ;
  tmp_v_r = sqrt(2.0 * energy_k * Q_ELEMENTAL / mass) ;
          
  //-- 速度：球面座標→デカルト座標
  tmp_v_x = tmp_v_r * cos(theta) * sin(psi) ;
  tmp_v_y = tmp_v_r * sin(theta) * sin(psi) ;
  tmp_v_z = tmp_v_r * cos(psi) ;

  /*
    cout << a11 << "\t" << a12 << "\t" << a13 << endl ;
    cout << a21 << "\t" << a22 << "\t" << a23 << endl ;
    cout << a31 << "\t" << a32 << "\t" << a33 << endl ;
  */

  //-- 座標変換
  pos_v.v_x = a11 * tmp_v_x + a12 * tmp_v_y  + a13 * tmp_v_z ;
  pos_v.v_y = a21 * tmp_v_x + a22 * tmp_v_y  + a23 * tmp_v_z ;
  pos_v.v_z = a31 * tmp_v_x + a32 * tmp_v_y  + a33 * tmp_v_z ;

  //-- 速度：デカルト座標→球面座標
  pos_v.v_r = sqrt(pos_v.v_x * pos_v.v_x + 
		   pos_v.v_y * pos_v.v_y + 
		   pos_v.v_z * pos_v.v_z) ;
  pos_v.v_theta = atan(pos_v.v_y / pos_v.v_x) ;
  pos_v.v_psi   = acos(pos_v.v_z / pos_v.v_r) ;

  //-- 例外処理（座標変換後、速度が以前のものを超えた場合、
  //   以前の速度とする）
  if( pos_v.v_r > tmp_v_r )
    {
      pos_v.v_r = tmp_v_r ;
      pos_v.v_x = tmp_v_r * cos(pos_v.v_theta) * sin(pos_v.v_psi) ;
      pos_v.v_y = tmp_v_r * sin(pos_v.v_theta) * sin(pos_v.v_psi) ;
      pos_v.v_z = tmp_v_r * cos(pos_v.v_psi) ;
    }
}


// イオン／原子の衝突: impact parameterを粒子の位置、速度から計算して求める。
void  Particle_class::
collision_accurate(double x_solid, double z_solid,
		   double *tmp_energy ) 
{
  double psi_c, psi ;
  double tmp_energy_loss ;
  struct Particle_location_velocity_struct tmp_pos_v ;

  // psi_c : 重心中心座標におけるscattering angle
  // psi   : laboratoryの座標におけるscattering angle
  // theta : scattering の回転角

  // -- 衝突前の進行方向軸からの相対的な速度
  //double tmp_v_x, tmp_v_y, tmp_v_z , tmp_v_r ;  

  // -- impact parameter の計算
  // 気体粒子の位置、速度を(x, z, vx, vz)、
  // 固体原子の位置を      (X, Z) とすると
  // impact parameter は
  // p = [(x - X)vz - (z - Z)vx]/√(vx^2 + vy^2)

  // ここでは最初に無次元のparameter を用意する（後で用いる）

  double p_no_dimension = 
    ((pos_v.x - x_solid) * pos_v.v_z - 
     (pos_v.z - z_solid) * pos_v.v_x) 
    / ( pos_v.v_x * pos_v.v_x + pos_v.v_z * pos_v.v_z)  ;

  double p_impact_parameter = 
    p_no_dimension * 
    sqrt( pos_v.v_x * pos_v.v_x + pos_v.v_z * pos_v.v_z)  ;

  // -- collision: S-W potential で、斥力が働く場合のみ
  if(fabs(p_impact_parameter) <= LENGTH_SW_MINIMUM)
    {
      tmp_energy_loss = put_energy_loss(Cl, Si, *tmp_energy ,
					fabs(p_impact_parameter),
					&psi_c, &psi ) ;
          
      // -- 上記の p の正負によって、気体粒子の位置／速度
      //    ベクトルのどちら側に固体原子が位置しているか判断できる
      //    （ r > 0 : 左側：この座標系で負の回転方向
      //       r < 0 : 右側：この座標系で正の回転方向）
      //    これに応じて回転変換を行う。
      
      tmp_pos_v = pos_v ;
      psi       = fabs(psi) ; //== 絶対値にしておく ==

      //=== エネルギー変化 -> 速度変化への変換 ===
      tmp_pos_v.v_x = 
 	sqrt(fabs(*tmp_energy - tmp_energy_loss)/(*tmp_energy)) * pos_v.v_x ;
      tmp_pos_v.v_z = 
 	sqrt(fabs(*tmp_energy - tmp_energy_loss)/(*tmp_energy)) * pos_v.v_z ;
        
      //--- エネルギー消失
      *tmp_energy -= tmp_energy_loss ;

      if(p_impact_parameter > 0)
	{
	  pos_v.v_x = tmp_pos_v.v_x *  cos(psi) + tmp_pos_v.v_z * sin(psi);
	  pos_v.v_z = tmp_pos_v.v_x *(-sin(psi))+ tmp_pos_v.v_z * cos(psi);
	}
      else 
	{
	  pos_v.v_x = tmp_pos_v.v_x *  cos(psi) + tmp_pos_v.v_z *(-sin(psi));
	  pos_v.v_z = tmp_pos_v.v_x *  sin(psi) + tmp_pos_v.v_z * cos(psi);
	}
      
      pos_v.x = x_solid +   pos_v.v_z * p_no_dimension ; 
      pos_v.z = z_solid +(- pos_v.v_x)* p_no_dimension ; 

      //-- 速度：デカルト座標→球面座標
      pos_v.v_r = sqrt(pos_v.v_x * pos_v.v_x + 
		       pos_v.v_y * pos_v.v_y + 
		       pos_v.v_z * pos_v.v_z) ;
      pos_v.v_theta = atan(pos_v.v_y / pos_v.v_x) ;
      pos_v.v_psi   = acos(pos_v.v_z / pos_v.v_r) ;
      
      //-- 新たに配置するのはSi原子からpだけ離れた地点 
      //   並進運動させないと衝突の処理がうまくいかないので
      //   CELL_SIZE だけ進ませる
      move_trans(CELL_SIZE);
    }
  
}

// イオン／原子の衝突：3次元計算
void  Particle_class::
collision_accurate3D(double x_solid, double y_solid,
		     double z_solid, double *tmp_energy ) 
{
  double psi_c, psi ;
  double tmp_energy_loss ;
 
  double uv_x, uv_y, uv_z ; //== 速度方向の単位ベクトル ==
  double lambda ;  // パラメータλ (particle.h参照)

  double  r_x,  r_y,  r_z ; //== Si原子 -> impact 方向のベクトル
  double ur_x, ur_y, ur_z ; //== Si原子 -> impact 方向の単位ベクトル

  // -- impact parameter の計算
  // 気体粒子の位置、速度を(x, y, z, vx, vy, vz)、
  // 固体原子の位置を      (X, Y, Z) とすると
  // impact parameter は
  // p = [(x - X)^2 + (y - Y)^2 + (z - Z)^2
  //     - [(X - x)vx + (Y - y)vy + (Z - z)vz]^2/(vx^2 + vy^2 + vz^2)]^(-1/2)

  double p_impact = 
    sqrt( (pos_v.x - x_solid) * (pos_v.x - x_solid) +
	  (pos_v.y - y_solid) * (pos_v.y - y_solid) +
	  (pos_v.z - z_solid) * (pos_v.z - z_solid) 
	  - ((x_solid - pos_v.x) * pos_v.v_x +
	     (y_solid - pos_v.y) * pos_v.v_y +
	     (z_solid - pos_v.z) * pos_v.v_z )  
	  * ((x_solid - pos_v.x) * pos_v.v_x +
	     (y_solid - pos_v.y) * pos_v.v_y +
	     (z_solid - pos_v.z) * pos_v.v_z )
	  / (pos_v.v_x * pos_v.v_x +
	     pos_v.v_y * pos_v.v_y +
	     pos_v.v_z * pos_v.v_z ) );

  // -- collision: S-W potential で、斥力が働く場合のみ
  if(fabs(p_impact) <= LENGTH_SW_MINIMUM)
    {
      //== 2005/12/11 debugged: 直交座標->極座標の変換に問題があった（例外処理の部分）
      double tmp_v_before, tmp_v_after ; //== 衝突前後の速度 ==
      tmp_v_before = sqrt(2.0 * (*tmp_energy) * Q_ELEMENTAL / mass) ; 

      tmp_energy_loss = put_energy_loss(Cl, Si, *tmp_energy ,
					fabs(p_impact),
					&psi_c, &psi ) ;
          
      //--- エネルギー消失
      *tmp_energy -= tmp_energy_loss ;

      tmp_v_after = sqrt(2.0 * (*tmp_energy) * Q_ELEMENTAL / mass) ; 
      psi         = fabs(psi) ; //== 絶対値にしておく ==

      //--- 単位ベクトル計算
      uv_x = pos_v.v_x / sqrt(pos_v.v_x * pos_v.v_x +
			      pos_v.v_y * pos_v.v_y +
			      pos_v.v_z * pos_v.v_z ) ;
      uv_y = pos_v.v_y / sqrt(pos_v.v_x * pos_v.v_x +
			      pos_v.v_y * pos_v.v_y +
			      pos_v.v_z * pos_v.v_z ) ;
      uv_z = pos_v.v_z / sqrt(pos_v.v_x * pos_v.v_x +
			      pos_v.v_y * pos_v.v_y +
			      pos_v.v_z * pos_v.v_z ) ;
      
      lambda = ((x_solid - pos_v.x) * pos_v.v_x +
		(y_solid - pos_v.y) * pos_v.v_y +
		(z_solid - pos_v.z) * pos_v.v_z ) 
	/  (pos_v.v_x * pos_v.v_x +
	    pos_v.v_y * pos_v.v_y +
	    pos_v.v_z * pos_v.v_z ) ;
      
      r_x = pos_v.x - x_solid + lambda * pos_v.v_x ;
      r_y = pos_v.y - y_solid + lambda * pos_v.v_y ;
      r_z = pos_v.z - z_solid + lambda * pos_v.v_z ;

      ur_x = r_x / sqrt(r_x * r_x + r_y * r_y + r_z * r_z);
      ur_y = r_y / sqrt(r_x * r_x + r_y * r_y + r_z * r_z);
      ur_z = r_z / sqrt(r_x * r_x + r_y * r_y + r_z * r_z);

      pos_v.v_x = tmp_v_after * (uv_x * cos(psi) + ur_x * sin(psi)) ;
      pos_v.v_y = tmp_v_after * (uv_y * cos(psi) + ur_y * sin(psi)) ;
      pos_v.v_z = tmp_v_after * (uv_z * cos(psi) + ur_z * sin(psi)) ;
     
      //== 新たな位置の決定 ===
      pos_v.x = x_solid 
	+   p_impact * (uv_x * (-sin(psi)) + ur_x * cos(psi)); 
      pos_v.y = y_solid 
	+   p_impact * (uv_y * (-sin(psi)) + ur_y * cos(psi)); 
      pos_v.z = z_solid 
	+   p_impact * (uv_z * (-sin(psi)) + ur_z * cos(psi));       

      //-- 速度：デカルト座標→球面座標 : 問題あり：atanは1対1の関数でない？
      pos_v.v_r = sqrt(pos_v.v_x * pos_v.v_x + 
		       pos_v.v_y * pos_v.v_y + 
		       pos_v.v_z * pos_v.v_z) ;
      pos_v.v_theta = atan(pos_v.v_y / pos_v.v_x) ;
      pos_v.v_psi   = acos(pos_v.v_z / pos_v.v_r) ;
      
      //-- 例外処理（座標変換後、速度が以前のものを超えた場合、
      //   以前の速度とする）
      if( pos_v.v_r > tmp_v_before )
	{
	  pos_v.v_x = tmp_v_before * (pos_v.v_x / pos_v.v_r);
	  pos_v.v_y = tmp_v_before * (pos_v.v_y / pos_v.v_r);
	  pos_v.v_z = tmp_v_before * (pos_v.v_z / pos_v.v_r);
	  pos_v.v_r = tmp_v_before ;
	  std::cout << "caution!\n";
	}
      //-- 新たに配置するのはSi原子からpだけ離れた地点 
      //   並進運動させないと衝突の処理がうまくいかないので
      //   CELL_SIZE だけ進ませる
      move_trans(CELL_SIZE);
    }
  
}
