
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

//== ����ʬ�ۤ�ե����뤫��������� =====
//   ��������ͳ�ˤ�ꡢflux_model ���齤����ɬ��
void Particle_class::read_angular_df(char filename[])
{
  //== microstructure ��Ʊ������
  // distribution ���Ѵ�����ɬ�פ����롣
  iedf_to_distribution(filename , n_adf_array ,
		       angle_array_degree , 
		       angular_df  	) ;
  /*
    input_array(filename ,angle_array_degree , angular_df,
    n_adf_array ) ;
    // (��������ɬ��!)
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
  // <- 0.0 �ˤ���ȡ����륻�������Ӥ����ơפ��ޤ�

  pos_v.v_x = 0.0 ;
  pos_v.v_y = 0.0 ;
  pos_v.v_z = v_z ;

  pos_v.v_r     = v_z ;
  pos_v.v_theta = 0.0 ;
  pos_v.v_psi   = 0.0 ;

  flag_contact  = false ; 
}

//-- γ�����͡��Ф������˽��®�٤�Ϳ����
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
  //== �������ν��� ==
  if(flag_boundary == NO_PERIODIC)//--�������ʤ�
    {
      if(pos_v.x < 0 || pos_v.x > SYSTEM_WIDTH_X)
	flag_inside = false ;
    }
  else if(flag_boundary == SPECULAR_REFLECT)//--�����Ƕ���ȿ��
    {
      if(pos_v.x < 0) // ������Ф����
	{
	  pos_v.x   = - pos_v.x ;
	  pos_v.v_x = - pos_v.v_x ;
	}
      if(pos_v.x > SYSTEM_WIDTH_X ) // ������Ф����
	{
	  pos_v.x   = 2.0 * SYSTEM_WIDTH_X - pos_v.x ;
	  pos_v.v_x = - pos_v.v_x ;
	}
    }
  else // -- �����������
    {
      if(pos_v.x < 0) // ������Ф����
	pos_v.x = SYSTEM_WIDTH_X - fabs(pos_v.x) ;
      
      if(pos_v.x > SYSTEM_WIDTH_X ) // ������Ф����
	pos_v.x = pos_v.x - SYSTEM_WIDTH_X ;
    }
  //-- ���Ԥ�����
  if(pos_v.y < 0)
    pos_v.y = L_INTER_ATOMIC + pos_v.y ;

  if(pos_v.y > L_INTER_ATOMIC)
    pos_v.y = pos_v.y - L_INTER_ATOMIC ;
    
  // -- �岼�˽Ф���flag��false�ˤ���
  if(pos_v.z < 0 || 
     pos_v.z > SYSTEM_HEIGHT_Z )
    flag_inside = false ;
  else
    flag_inside = true ;
}

// -- pos_v_togo �����
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
  // -- �������������θ
  if(pos_v_togo.x < 0) // ������Ф����
    pos_v_togo.x = SYSTEM_WIDTH_X - fabs(pos_v_togo.x) ;

  if(pos_v_togo.x > SYSTEM_WIDTH_X ) // ������Ф����
    pos_v_togo.x = pos_v_togo.x - SYSTEM_WIDTH_X ;
    
  // -- ��˽Ф���FALSE���֤�
  if(pos_v_togo.z < 0 || 
     pos_v_togo.z > SYSTEM_HEIGHT_Z )
    flag_inside_togo = false ;
  else
    flag_inside_togo =  true ;
}

//-- ®�٤��ѹ� (x, z ������ʬ�Τߡ�
void Particle_class::random_reflection(bool flag_cosine_dist ,
			 double normal_x, double normal_z)
{
  double theta ;
  if (flag_cosine_dist  == true )
    {
      // ��n + �� ����n��ˡ������������ cosine distribution
      // normal: ���� -> �������������� �Ǥ��뤳�Ȥ����
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
  
  //-- ®�١��ǥ���Ⱥ�ɸ�����̺�ɸ
  
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
 
  //double theta_c ; // �ſ���ɸ�ˤ�����ķ���֤��

  *theta_c = theta_integrate_NC(b_impact_para,
				energy_i ) ;

  // scattering angle �� (4-15)
  // �� = arctan{sin�� / [cos�� + (M1/M2)]}
  
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

  // psi_c : �ſ��濴��ɸ�ˤ�����scattering angle
  // psi   : laboratory�κ�ɸ�ˤ�����scattering angle
  // theta : scattering �β�ž�ѡʥ�����˷����

  // -- �������οʹ����������������Ū��®��
  double tmp_v_x, tmp_v_y, tmp_v_z , tmp_v_r ;  

  // �Ѵ��Τ����ľ�����
  double a11, a12, a13 ;
  double a21, a22, a23 ;
  double a31, a32, a33 ;
  bool    sign_flag ;  //  sqrt������

  // -- �������οʹ��������ˤ�ä��Ѵ���ܤ���®��
  //  double v_x, v_y, v_z ;

  energy_k = mass * pos_v.v_r * pos_v.v_r 	      
    /  (2.0 * Q_ELEMENTAL) ;  

  // -- collision
  tmp_energy_loss = put_energy_loss(Cl, Si, energy_k ,
				    p_impact_parameter,
				    &psi_c, &psi ) ;
          
  // ���ͥ륮���ü�
  energy_k -= tmp_energy_loss ;
          
  //-- ®���Ѳ�����������®�������μ�����
  //���������������ͳѤ������˷���(2 * PI * RAN0() )
  // ������Ȥ�®�٤� x, y, z ��ʬ����ơ������
  // ��������®�٥٥��ȥ�򸵤˺�ɸ�Ѵ���Ԥ���
  
  //  ��ɸ�Ѵ����Ѵ�����η���
  //  (sqrt������-> ����ǡ�0.5�ʾ�ʤ�����0.5̤������
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
          
  //-- ®�١����̺�ɸ���ǥ���Ⱥ�ɸ
  tmp_v_x = tmp_v_r * cos(theta) * sin(psi) ;
  tmp_v_y = tmp_v_r * sin(theta) * sin(psi) ;
  tmp_v_z = tmp_v_r * cos(psi) ;

  /*
    cout << a11 << "\t" << a12 << "\t" << a13 << endl ;
    cout << a21 << "\t" << a22 << "\t" << a23 << endl ;
    cout << a31 << "\t" << a32 << "\t" << a33 << endl ;
  */

  //-- ��ɸ�Ѵ�
  pos_v.v_x = a11 * tmp_v_x + a12 * tmp_v_y  + a13 * tmp_v_z ;
  pos_v.v_y = a21 * tmp_v_x + a22 * tmp_v_y  + a23 * tmp_v_z ;
  pos_v.v_z = a31 * tmp_v_x + a32 * tmp_v_y  + a33 * tmp_v_z ;

  //-- ®�١��ǥ���Ⱥ�ɸ�����̺�ɸ
  pos_v.v_r = sqrt(pos_v.v_x * pos_v.v_x + 
		   pos_v.v_y * pos_v.v_y + 
		   pos_v.v_z * pos_v.v_z) ;
  pos_v.v_theta = atan(pos_v.v_y / pos_v.v_x) ;
  pos_v.v_psi   = acos(pos_v.v_z / pos_v.v_r) ;

  //-- �㳰�����ʺ�ɸ�Ѵ��塢®�٤������Τ�Τ�Ķ������硢
  //   ������®�٤Ȥ����
  if( pos_v.v_r > tmp_v_r )
    {
      pos_v.v_r = tmp_v_r ;
      pos_v.v_x = tmp_v_r * cos(pos_v.v_theta) * sin(pos_v.v_psi) ;
      pos_v.v_y = tmp_v_r * sin(pos_v.v_theta) * sin(pos_v.v_psi) ;
      pos_v.v_z = tmp_v_r * cos(pos_v.v_psi) ;
    }
}


// �����󡿸��Ҥξ���: impact parameter��γ�Ҥΰ��֡�®�٤���׻����Ƶ��롣
void  Particle_class::
collision_accurate(double x_solid, double z_solid,
		   double *tmp_energy ) 
{
  double psi_c, psi ;
  double tmp_energy_loss ;
  struct Particle_location_velocity_struct tmp_pos_v ;

  // psi_c : �ſ��濴��ɸ�ˤ�����scattering angle
  // psi   : laboratory�κ�ɸ�ˤ�����scattering angle
  // theta : scattering �β�ž��

  // -- �������οʹ����������������Ū��®��
  //double tmp_v_x, tmp_v_y, tmp_v_z , tmp_v_r ;  

  // -- impact parameter �η׻�
  // ����γ�Ҥΰ��֡�®�٤�(x, z, vx, vz)��
  // ���θ��Ҥΰ��֤�      (X, Z) �Ȥ����
  // impact parameter ��
  // p = [(x - X)vz - (z - Z)vx]/��(vx^2 + vy^2)

  // �����ǤϺǽ��̵������parameter ���Ѱդ���ʸ���Ѥ����

  double p_no_dimension = 
    ((pos_v.x - x_solid) * pos_v.v_z - 
     (pos_v.z - z_solid) * pos_v.v_x) 
    / ( pos_v.v_x * pos_v.v_x + pos_v.v_z * pos_v.v_z)  ;

  double p_impact_parameter = 
    p_no_dimension * 
    sqrt( pos_v.v_x * pos_v.v_x + pos_v.v_z * pos_v.v_z)  ;

  // -- collision: S-W potential �ǡ����Ϥ�Ư�����Τ�
  if(fabs(p_impact_parameter) <= LENGTH_SW_MINIMUM)
    {
      tmp_energy_loss = put_energy_loss(Cl, Si, *tmp_energy ,
					fabs(p_impact_parameter),
					&psi_c, &psi ) ;
          
      // -- �嵭�� p ������ˤ�äơ�����γ�Ҥΰ��֡�®��
      //    �٥��ȥ�Τɤ���¦�˸��θ��Ҥ����֤��Ƥ��뤫Ƚ�ǤǤ���
      //    �� r > 0 : ��¦�����κ�ɸ�Ϥ���β�ž����
      //       r < 0 : ��¦�����κ�ɸ�Ϥ����β�ž������
      //    ����˱����Ʋ�ž�Ѵ���Ԥ���
      
      tmp_pos_v = pos_v ;
      psi       = fabs(psi) ; //== �����ͤˤ��Ƥ��� ==

      //=== ���ͥ륮���Ѳ� -> ®���Ѳ��ؤ��Ѵ� ===
      tmp_pos_v.v_x = 
 	sqrt(fabs(*tmp_energy - tmp_energy_loss)/(*tmp_energy)) * pos_v.v_x ;
      tmp_pos_v.v_z = 
 	sqrt(fabs(*tmp_energy - tmp_energy_loss)/(*tmp_energy)) * pos_v.v_z ;
        
      //--- ���ͥ륮���ü�
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

      //-- ®�١��ǥ���Ⱥ�ɸ�����̺�ɸ
      pos_v.v_r = sqrt(pos_v.v_x * pos_v.v_x + 
		       pos_v.v_y * pos_v.v_y + 
		       pos_v.v_z * pos_v.v_z) ;
      pos_v.v_theta = atan(pos_v.v_y / pos_v.v_x) ;
      pos_v.v_psi   = acos(pos_v.v_z / pos_v.v_r) ;
      
      //-- ���������֤���Τ�Si���Ҥ���p����Υ�줿���� 
      //   �¿ʱ�ư�����ʤ��Ⱦ��ͤν��������ޤ������ʤ��Τ�
      //   CELL_SIZE �����ʤޤ���
      move_trans(CELL_SIZE);
    }
  
}

// �����󡿸��Ҥξ��͡�3�����׻�
void  Particle_class::
collision_accurate3D(double x_solid, double y_solid,
		     double z_solid, double *tmp_energy ) 
{
  double psi_c, psi ;
  double tmp_energy_loss ;
 
  double uv_x, uv_y, uv_z ; //== ®��������ñ�̥٥��ȥ� ==
  double lambda ;  // �ѥ�᡼���� (particle.h����)

  double  r_x,  r_y,  r_z ; //== Si���� -> impact �����Υ٥��ȥ�
  double ur_x, ur_y, ur_z ; //== Si���� -> impact ������ñ�̥٥��ȥ�

  // -- impact parameter �η׻�
  // ����γ�Ҥΰ��֡�®�٤�(x, y, z, vx, vy, vz)��
  // ���θ��Ҥΰ��֤�      (X, Y, Z) �Ȥ����
  // impact parameter ��
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

  // -- collision: S-W potential �ǡ����Ϥ�Ư�����Τ�
  if(fabs(p_impact) <= LENGTH_SW_MINIMUM)
    {
      //== 2005/12/11 debugged: ľ���ɸ->�˺�ɸ���Ѵ������꤬���ä����㳰��������ʬ��
      double tmp_v_before, tmp_v_after ; //== ���������®�� ==
      tmp_v_before = sqrt(2.0 * (*tmp_energy) * Q_ELEMENTAL / mass) ; 

      tmp_energy_loss = put_energy_loss(Cl, Si, *tmp_energy ,
					fabs(p_impact),
					&psi_c, &psi ) ;
          
      //--- ���ͥ륮���ü�
      *tmp_energy -= tmp_energy_loss ;

      tmp_v_after = sqrt(2.0 * (*tmp_energy) * Q_ELEMENTAL / mass) ; 
      psi         = fabs(psi) ; //== �����ͤˤ��Ƥ��� ==

      //--- ñ�̥٥��ȥ�׻�
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
     
      //== �����ʰ��֤η��� ===
      pos_v.x = x_solid 
	+   p_impact * (uv_x * (-sin(psi)) + ur_x * cos(psi)); 
      pos_v.y = y_solid 
	+   p_impact * (uv_y * (-sin(psi)) + ur_y * cos(psi)); 
      pos_v.z = z_solid 
	+   p_impact * (uv_z * (-sin(psi)) + ur_z * cos(psi));       

      //-- ®�١��ǥ���Ⱥ�ɸ�����̺�ɸ : ���ꤢ�ꡧatan��1��1�δؿ��Ǥʤ���
      pos_v.v_r = sqrt(pos_v.v_x * pos_v.v_x + 
		       pos_v.v_y * pos_v.v_y + 
		       pos_v.v_z * pos_v.v_z) ;
      pos_v.v_theta = atan(pos_v.v_y / pos_v.v_x) ;
      pos_v.v_psi   = acos(pos_v.v_z / pos_v.v_r) ;
      
      //-- �㳰�����ʺ�ɸ�Ѵ��塢®�٤������Τ�Τ�Ķ������硢
      //   ������®�٤Ȥ����
      if( pos_v.v_r > tmp_v_before )
	{
	  pos_v.v_x = tmp_v_before * (pos_v.v_x / pos_v.v_r);
	  pos_v.v_y = tmp_v_before * (pos_v.v_y / pos_v.v_r);
	  pos_v.v_z = tmp_v_before * (pos_v.v_z / pos_v.v_r);
	  pos_v.v_r = tmp_v_before ;
	  std::cout << "caution!\n";
	}
      //-- ���������֤���Τ�Si���Ҥ���p����Υ�줿���� 
      //   �¿ʱ�ư�����ʤ��Ⱦ��ͤν��������ޤ������ʤ��Τ�
      //   CELL_SIZE �����ʤޤ���
      move_trans(CELL_SIZE);
    }
  
}
