
//-- obsolete
 // sin^2(��/2) ���ͤ�ޤ�����(4-13)
  /*
  double tmp_sin2 = 1.0 / (1.0 + 4.0 * reduced_energy * reduced_energy *
			   b_impact_para * b_impact_para) ;

  // scattering angle �� (4-15)
  // �� = arctan{sin�� / [cos�� + (M1/M2)]}

  // cos�� = 1 - 2 sin^2(��/2), sin�� = ��(1 - cos^2��) ������
  double tmp_cos = 1.0 - 2.0 * tmp_sin2 * tmp_sin2 ;

  *psi_scattering = atan( sqrt(1.0 - tmp_cos * tmp_cos) /
			  (tmp_cos + Ion.atomic_weight / Solid.atomic_weight_solid)) ;

  // T : (4-14)
  return 
    4.0 * Ion.atomic_weight * Solid.atomic_weight_solid * energy_i
    * tmp_sin2 / ((Ion.atomic_weight + Solid.atomic_weight_solid) *
		  (Ion.atomic_weight + Solid.atomic_weight_solid) ) ;
  */

 //-- ��ɸ�Ѵ�
/*
          v_x = tmp_v_x * cos(Cl_ion.pos_v.v_theta) * cos(Cl_ion.pos_v.v_psi) +
            tmp_v_y * (- sin(Cl_ion.pos_v.v_theta) ) +
              tmp_v_z * (- cos(Cl_ion.pos_v.v_theta)) * sin(Cl_ion.pos_v.v_psi) ;
      
          v_y = tmp_v_x * sin(Cl_ion.pos_v.v_theta) * cos(Cl_ion.pos_v.v_psi) +
            tmp_v_y * cos(Cl_ion.pos_v.v_theta)  +
              tmp_v_z * sin(Cl_ion.pos_v.v_theta) * sin(Cl_ion.pos_v.v_psi) ;
          
          v_z = tmp_v_x * sin(Cl_ion.pos_v.v_psi) +
            tmp_v_z * cos(Cl_ion.pos_v.v_psi) ;
  
  // x-z ʿ�̤Ǧײ�ž��x-yʿ�̤ǦȲ�ž����������
          // ��                               ��
          // | cos��cos��  -sin��  -cos��sin�� |  
          // | sin��cos��   cos��   sin��sin�� |
          // |      sin��      0         cos�� |
          // ��                               ��
          // ��                               ��
          // | cos��cos��  -cos��sin��  -sin�� |  
          // | sin��             cos��      0  |
          // | sin��cos��  -sin��sin��   cos�� |
          // ��                               ��
  
*/




//  tmp_x = Cl_neutral.pos_v.x ;// ���Υ��ƥåפΰ��֤�Ͽ  
//tmp_z = Cl_neutral.pos_v.z ;
// Cl_neutral.pos_v.x = tmp_x ;
// Cl_neutral.pos_v.z = tmp_z ;


//=== included_to_main2.cc ���SiCl2�ν���

 /*
  i_exception_move = 0 ;
  while(1) 
    {
      //--- �ΰ������ˤ��ޤ�Ĺ���ȤɤޤäƤ����硢
      //    �㳰������γ�Ҥ�õ��
      i_exception_move++ ;
      if( i_exception_move > N_CELL_Z * 200 )
	break ;
      
      // -- �¿ʱ�ư
      flag_out_of_domain          
	= SiCl2_neutral.move_trans(FREE_FLIGHT_PATH_Si) ;
      if (flag_out_of_domain == FALSE) // -- ���˽Ф���
	break ;
      
      // -- γ�ҤΤ������phase�����
      SiCl2_neutral.get_position() ;
	      
      // -- ����Ǥ���С�������Υ����ֹ����� -> deposit
      if(Shape.shape_matrix
	 [SiCl2_neutral.i_cell_x][SiCl2_neutral.i_cell_z] != SHAPE_SPACE )
	{
	  SiCl2_neutral.get_previous_cell() ; // �����ֹ����
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
	  //��ȿ�� : 
	  // ̵�¥롼�פ򤵤��뤿��ˡ������󥿤��Ѱդ��롣
	  // γ�Ҥ򲾤˿ʤޤ��Ƥߤơ�(pos_v_togo)
	  // �⤷�ʤ���褬������ʬ�Ǥʤ���С�������ä�
	  //     ���ľ��(50��ˤʤä���⤦����)
	  i_exception_reflect = 0;
	  Shape.get_surfacenormal(SiCl2_neutral.i_cell_x, 
				  SiCl2_neutral.i_cell_z) ;
	reflection_SiCl2:
	  if( i_exception_reflect > 50 )
	    break ;
	  
	  SiCl2_neutral.random_reflection //(FALSE, 0.0, 0.0 ) ;
	  (TRUE,  // ��ȿ��
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
      
    } //=== while loop ��λ�ʣ�γ�Ҥ�trajectory��
  */


//=================================
// monte_carlo.cc ��Υ��å��󥰤ν���
// ��Ion_class::ion_enhanced_etch �ؿ����֤�������
//***** 02Dec2003 **************
//--- ˰�¾��֡ʵ���Cl �����Ĥ���ˤǤ���С�æΥ������
//	   --�ʤ��ġ����ޤ�æΥȿ���ˤ�����äƤ��ʤ���С�
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
			      Cl_ion.pos_v.v_x,//<- pos_v.x �Ȥ��Ƥ����Τǡ�������Ÿ���Ѥ��ä���29Aug2004 
			      Cl_ion.pos_v.v_z ),
     Shape.n_oxygen[i_x][i_z] ) ;
#endif
  if(RAN0() < tmp_beta / YIELD_MAX) // <- ����γ�Ҥ��ΨŪ�˰���
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
// monte_carlo.cc ��Υ������ɽ�̾��ͤν���
// ��Ion_class::impact_on_surface �ؿ����֤�������
  // -- ���Υ���� Si ������С�����
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
		  
		  // impact parameter �����X�Ƿ��ꡧ
		  // p = (L/2)��X  where  L �ϸ��Ҵֵ�Υ
		  p_impact_parameter = (L_INTER_ATOMIC / 2.0) * sqrt(RAN0());
		  
		  // -- collision
		  Cl_ion.collision_with_solid_atom( p_impact_parameter,
						    &tmp_energy );
		}
	      // == �ޥ����Ǥ���С�break
	      else if(Shape.put_shape_matrix(Cl_ion.pos_v.x, 
					     Cl_ion.pos_v.z ,
					     &i_x, &i_z ) == HARD_MASK )
		{
		  flag_out_of_domain = false ;
		  break ;
		}
	      */
	      //=== ����  =============
	      //cout << Cl_ion.pos_v.x << "\t"
	      //<< Cl_ion.pos_v.y << "\t"
	      //<< Cl_ion.pos_v.z << endl ;


//=================================
// monte_carlo.cc �������γ�Ҥν���
 //--- ���� ---
	  /*
	  Cl_neutral.inject_from_top
	    ( sqrt(2.0 * INCIDENT_ENERGY_NEUTRAL * 
		   Q_ELEMENTAL / Cl_neutral.mass) ) ; 
	  
	  i_exception_move = 0 ;
	  while(1) 
	    {
	      //-- �ΰ������ˤ��ޤ�Ĺ���ȤɤޤäƤ����硢
	      //   �㳰������γ�Ҥ�õ��
	      i_exception_move++ ;
	      if( i_exception_move > N_CELL_Z * 200 )
		break ;
		
	      // -- �ޤ��¿ʱ�ư ----
	      flag_out_of_domain 
		= Cl_neutral.move_trans(FREE_FLIGHT_PATH_Si) ;
	      if (flag_out_of_domain == false) // -- ���˽Ф���
		break ;
	     
	      //--- ���뤬����Ǥʤ��� ---
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

		  //---- ����Ū���å��� ----
		  // YIELD_CHEMICAL�γ�Ψ��æΥ�������� 
		  if( RAN0() < YIELD_CHEMICAL &&
		      FLAG_CHEMICAL_ETCH == true  )
		    {
		      Shape.desorb_Si_chemical_etch(i_x, i_z);
		    }
		  else //---- ���� -----
		    {
		      flag_adsorption = 
			Shape.settle_Cl_into_bareSi(i_x, i_z,  // ���夷����
						    ADSORPTION_PROBABILITY );
		      if(flag_adsorption == true )
			break;   // ���夷����롼�פ�Ф�
		    }
		  // ��ȿ�� : 
		  // ̵�¥롼�פ򤵤��뤿��ˡ������󥿤��Ѱդ��롣
		  // γ�Ҥ򲾤˿ʤޤ��Ƥߤơ�(pos_v_togo)
		  // �⤷�ʤ���褬������ʬ�Ǥʤ���С�������ä�
		  //     ���ľ��(50��ˤʤä���⤦����)
		  i_exception_reflect = 0;
		  Shape.get_surfacenormal(i_x, i_z) ;
		reflection:
		  if( i_exception_reflect > 50 )
		    break ;

		  Cl_neutral.random_reflection //(FALSE,0.0,0.0);//��ȿ��
		   (true,//<-cosine diffusion ��Ƴ������Ⱥ������оΤˤʤ�<- this bug was fixed (29Aug2004)
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
	    }//=== while loop ��λ
	  */

const double GATE_WIDTH   = 100.0 * 1.0e-9 ;//== ��������(nm) ==
const double MARGIN_WIDTH =  50.0 * 1.0e-9 ;//== �ΰ躸����;��(nm)
=     int(floor(GATE_WIDTH / CELL_SIZE) + 0.5)  
  + 2*int(floor(MARGIN_WIDTH/CELL_SIZE) + 0.5) ; //floor�ؿ��ˤ��ͼθ���

const int    N_LEFT_MASK  = int(floor(MARGIN_WIDTH/CELL_SIZE) + 0.5);//60;//35;//70;// �ޥ��������̤ΰ���
const int    N_RIGHT_MASK = N_CELL_X - int(floor(MARGIN_WIDTH/CELL_SIZE) + 0.5);//60;//35;//70;// �ޥ��������̤ΰ���
