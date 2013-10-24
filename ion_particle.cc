
#include <iostream>
#include "particle.h"
#include "random_number.h"
#include "iedf_to_distribution.h"

#include "ion_particle.h"

//==  ���ͥ��ͥ륮�����Ф��� etch yield ���֤� ==
//    yield �η����� = C(��E - ��Eth): C = 0.77, Eth = 20.0 
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
//== ������ȿ�ͤ������Υ��å��󥰼�Ψ�� Er < Eth �Ǥ���� C(��Ei - ��Eth),
//                                      Er > Eth �Ǥ���� C(��Ei - ��Er)�Ȥ���
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
  velocity_file.open(NULL);//== ������®�٤�Ͽ ==

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

  // randomized_distribution   (long int ����μ�, int ����������ο�,
  //   double[] ���ͥ륮��  ,double[] ʬ��,
  //   -> ����  double[] (�������), double[] ����줿���) 
  
  // -- ưŪ���ݡ�
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

  // X��ɸ������ˤ�äƷ��ꤹ��
  pos_v.x = RAN0() * SYSTEM_WIDTH_X ; 
  pos_v.y = 0.0 ;
  pos_v.z = 0.0 ;
  
  // --- ���ͳѤϤ��λ����Ƿ��ꤷ�Ƥ���(randomized_angle_array)
  // ����y ������β�ž�Ѥ�����Ǻ�� 
  // ���ͳѦס���ž�ѦȤȤ���� z = v cos��, x = v sin��cos��, y = v sin��sin��
      
  // randomized_energy �� eV ��ñ�̤ʤΤ� Q_ELEMENTAL �򤫤���Volt ñ�̤ˤ���

  pos_v.v_r     =  v ; // <--- Absolute velocity

  pos_v.v_theta =  2.0 * PI * RAN0() ; 
     
  pos_v.v_psi   = 
    PI * randomized_angle_array[0] / 180.0 ; // degree -> rad�ش��� 
  
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

  //== ���å������ѿ������ ==
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
	
//== ������ξ��ͤˤ�륨�å��󥰤ν���
//bug fix:<- pos_v.x �Ȥ��Ƥ����Τǡ�������Ÿ���Ѥ��ä�:29Aug2004 
//bug fix(*): 28/Nov/2005
//  ����������ͳѤη����ָ������Ǥ�®��(pos_v)�פǹԤäƤ���
void Ion_class::
ion_enhanced_etch(class  Shape_trim_class *Shape_trim,
		  double incident_energy, double reflected_energy,
		  Particle_location_velocity_struct pos_v_recorded,
		  int i_x, int i_z ) 
{
  double tmp_beta = 0.0;
  if(FLAG_INCIDENT_ANGLE == true)
    {
      //== (3rd paper)���Ǥ����夹�뤳�Ȥˤ������԰�����ɤ����ᡢ
      //   ���å��󥰼�Ψ��ʬ�������롣nearest-neighbor�Υ���˴ޤޤ����Ǥο����פ���
      //   6 �ǳ�ä�coverage���Ȥ��롣
      //   ��Ψ�� Y = �� Y(SiO2) + (1 - ��)Y(Si)
      
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

  if(RAN0() < tmp_beta / YIELD_MAX) // <- ����γ�Ҥ��ΨŪ�˰���
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
  //etch_yield_angle�ؿ��κǸ�ΰ�����2 �Ȥ��뤳�Ȥǻ�����γ��ٰ�¸����������
 
  if(RAN0() < tmp_beta / YIELD_MAX) // <- ����γ�Ҥ��ΨŪ�˰���
    {
      Shape_trim->desorb_mask(i_x, i_z);
      record_desorption(pos_v_recorded, Shape_trim->n_oxygen[i_x][i_z]) ; 
    }
}

//==   �ѥ�����ɽ�̤���ã�������ɤ�����Ƚ�� 
bool Ion_class::
impact_on_surface(class  Shape_trim_class *Shape_trim)
{
  int i_x, i_z , i_x_particle, i_z_particle ;
  // -- ���Υ���� Si ������С�����
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
      
      //== ��ɽ�̤ǤΤߥ��å��󥰽�����Ԥ����ᡢ�ܿ����Ƥ��뤫
      //   �ե饰��Ƚ��
      if(flag_contact == false)
	{
	  ion_enhanced_etch(Shape_trim, incident_energy , 0.0,  pos_v, i_x, i_z ) ;
	  flag_contact = true ; //== ���������ǤΥ��å��󥰤Ϥʤ� ==
	}
		  
      // impact parameter �����X�Ƿ��ꡧ
      // p = (L/2)��X  where  L �ϸ��Ҵֵ�Υ
      p_impact_parameter = (L_INTER_ATOMIC / 2.0) * sqrt(RAN0());
		  
      // -- collision
      collision_with_solid_atom( p_impact_parameter );
    }
  // == �ޥ����Ǥ���С�break
  else if(Shape_trim->put_shape_matrix(pos_v.x, pos_v.z ,
				       &i_x, &i_z ) == HARD_MASK )
    {
      flag_inside = false ;
      return true ;
    }
  return false ;
}

//==   ���Ʊ�������ͤν�������������Τ�ޤ��
bool Ion_class::
impact_scattering(class  Shape_trim_class *Shape_trim,
		  bool  flag_mask_erosion) 
{
  int i_x, i_z, i_x_particle, i_z_particle ;
  double p_impact_parameter ;
  double incident_energy ;
  incident_energy = mass * pos_v.v_r * pos_v.v_r /(2.0 * Q_ELEMENTAL) ; 

  // == �ޥ����Ǥ���С�break
  if(Shape_trim->put_shape_matrix(pos_v.x, pos_v.z ,
				  &i_x, &i_z ) == HARD_MASK )
    {
      flag_inside = false ;
      return true ;
    }  
  // -- ���Υ���� Si ������С�����
  else if(Shape_trim->put_shape_matrix(pos_v.x, pos_v.z ,
				       &i_x, &i_z ) == SHAPE_Si )
    {
      //���ͥ륮���������ͤ򲼲�äƤ����� stop
      if(incident_energy <= ENERGY_ION_STOPPING)
	{
	  //== stop���������ǻ���Ū�˥��å��󥰤ν�����Ԥ�
	  ion_enhanced_etch(Shape_trim, energy_etch, 0.0, pos_v_etch,
			    i_x_etch , i_z_etch ) ;
	  return true ;
	}
      if(flag_contact == false)
	{
	  record_etch_position(i_x, i_z, incident_energy);
	  flag_contact = true ; 
	}
      
      // impact parameter �����X�Ƿ��ꡧ
      // p = (L/2)��X  where  L �ϸ��Ҵֵ�Υ
      p_impact_parameter = (L_INTER_ATOMIC / 2.0) * sqrt(RAN0());
      
      // -- collision --
      //collision_with_solid_atom( p_impact_parameter );
      collision_accurate3D( (i_x + 0.5) * L_INTER_ATOMIC,
			    0.5         * L_INTER_ATOMIC,
			    (i_z + 0.5) * L_INTER_ATOMIC,
			    &incident_energy );
    }
  // -- ���ܥ��뤬 Si or hardmask �ξ��
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
      // -- collision �ʼ������������θ��
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
       // -- collision �ʼ������������θ��
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
      if(flag_contact    == true)//== ������ʬ���ܿ����� -> detach�������Ȥ��̣����
	{
	  flag_reflection = true ; //== �������𤬵����ä����Ȥ�Ͽ
	  //velocity_file << (pos_v.v_psi * 180 / PI - 90.0 ) << "\t" 
	  //	<< incident_energy / 50.0 << "\n" ;

	  //== ������reflection�򵯤������������Ū�˥��å��󥰽�������
	  //   ��Ψ�ϥ��ͥ륮���Ѳ�ʬ
	  //  �ʥ��å��󥰾��ͻ��ݸ��ߤ�Ei�ˤΤߤˤ���ΤȤ���
	  //ion_enhanced_etch(Shape_trim, energy_etch, 0.0, pos_v_etch,
	  //	    i_x_etch , i_z_etch ) ;
	  ion_enhanced_etch(Shape_trim, energy_etch, incident_energy, pos_v_etch,
			    i_x_etch , i_z_etch ) ;
	}
      flag_contact    = false ;
    }
  
  return false ;
}
