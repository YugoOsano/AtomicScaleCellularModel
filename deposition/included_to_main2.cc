// included_to_main2.cc 
//=== monte_carlo.cc ��Ǥ� main �ˤ�����

//*************************************************
//---   etch product  æΥ ��Cl_ion ��flag �ˤ���
//*************************************************

if(FLAG_SiCl4_INCLUDED == true &&
   Cl_ion.flag_desorb_Si == true )
{
  bool flag_product_diffusion  ;

  //== ���Ǥ�ޤ���: ��ñ�Τ��� SiOCl2�ΤߤȤ���==
  /*
  if(Cl_ion.n_oxy_desorb_Si >= 1)
    {
      flag_product_diffusion =
	SiOCl2_neutral.allocate_on_surface
	(&Shape,  Cl_ion.position_at_desorption,
	 sqrt(2.0 * INCIDENT_ENERGY_NEUTRAL * 
	      Q_ELEMENTAL / SiCl4_neutral.mass)) ;
  
      if( flag_product_diffusion == true )
	SiOCl2_neutral.whole_flight(&Shape, false) ;
    }
  */
  //else //== SiCl4�ξ��
  // {
      flag_product_diffusion =
	SiCl4_neutral.allocate_on_surface
	(&Shape,  Cl_ion.position_at_desorption,
	 sqrt(2.0 * INCIDENT_ENERGY_NEUTRAL * 
	      Q_ELEMENTAL / SiCl4_neutral.mass)) ;
      
      if( flag_product_diffusion == true )
	SiCl4_neutral.whole_flight(&Shape, false) ;
      
      //*************************************************
      //---   SiCl2  ���� 
      //*************************************************
      if(FLAG_SiCl2_INCLUDED == true )
	{
	  //std::cout << i_ion << "\t" << SiCl4_neutral.put_flag_inside() << "\n" ;
	  //if( i_ion % ION_SiCl2_RATIO == 0 )
	  if(SiCl4_neutral.put_flag_inside() == false &&
	     RAN0() < SiCl2_neutral.reincidence_probability )
	    {
	      //  std::cout << "SiCl2 incidence\n" ;
	      SiCl2_neutral.inject_from_top(sqrt(2.0 * INCIDENT_ENERGY_NEUTRAL * 
						 Q_ELEMENTAL/SiCl2_neutral.mass)); 
	      SiCl2_neutral.add_n_injection();
	      SiCl2_neutral.whole_flight(&Shape, false) ;
	      
	    } // SiCl2 �ν�����λ
	}
      //*************************************************
      //---   SiCl2O  ���� 
      //*************************************************
      if(FLAG_SiOCl2_INCLUDED == true)
	{
	  //while(double(i_ion) * SiO_ION_RATIO > interval_step_SiO )
	  if(SiCl4_neutral.put_flag_inside() == false &&
	     RAN0() < SiOCl2_neutral.reincidence_probability )
	    {
	      //interval_step_SiO++ ;
	      
	      SiOCl2_neutral.inject_from_top
		( sqrt(2.0 * INCIDENT_ENERGY_NEUTRAL * 
		       Q_ELEMENTAL / SiOCl2_neutral.mass) ) ; 
	      SiOCl2_neutral.add_n_injection();
	      SiOCl2_neutral.whole_flight(&Shape, false) ;
	    } // SiOCl2 �ν�����λ
	  //=== ���Ϳ������󥿤ν��� ===
	  if(i_ion % INTERVAL_FILE_OUTPUT == 0 ||
	     i_ion == N_ION_INJECT )
	    file_SiOCl2 << SiOCl2_neutral.put_n_injection() << "\n" ;
	  
	}
      //    }
}  
//=== ���Ϳ������󥿤ν��� ===
if(i_ion % INTERVAL_FILE_OUTPUT == 0 ||
   i_ion == N_ION_INJECT )
{
  if(FLAG_SiCl2_INCLUDED == true)
    file_SiCl2 << i_ion << "\t" << SiCl2_neutral.put_n_injection() << "\n" ;
  if(FLAG_SiOCl2_INCLUDED == true)
    file_SiOCl2 << i_ion << "\t"<< SiOCl2_neutral.put_n_injection() << "\n" ;
}

//*************************************************
//---   ���Ǹ���  ���� ��1/OXYGEN_ION_RATIO��ˣ����
//*************************************************

if(FLAG_OXYGEN_INCLUDED == true)
{
  //if( i_ion % ION_OXYGEN_RATIO == 0 )
  // ==== f���1���f�Ͼ����˽�����Ԥ����ᡢ�����Τ褦�ʼ�³����Ԥ�
  while(double(i_ion) * OXYGEN_ION_RATIO > interval_step_oxygen )
    {
      interval_step_oxygen++ ;

      Oxygen_atom.inject_from_top( sqrt(2.0 * INCIDENT_ENERGY_NEUTRAL * 
					Q_ELEMENTAL / Oxygen_atom.mass) ) ; 
      
      i_exception_move = 0 ;
      while(1) 
	{
	  //--- �ΰ������ˤ��ޤ�Ĺ���ȤɤޤäƤ����硢
	  //    �㳰������γ�Ҥ�õ��
	  i_exception_move++ ;
	  if( i_exception_move > N_CELL_Z * 200 )
	    break ;
      
	  // -- �¿ʱ�ư
	  
	  Oxygen_atom.move_trans(FREE_FLIGHT_PATH_Si) ;
	  if (Oxygen_atom.put_flag_inside() == false) // -- ���˽Ф���
	    break ;
      
	  // -- γ�ҤΤ������phase�����
	  Oxygen_atom.get_position() ;

	      
	  // -- ���֤Ǥʤ���С������ʤ⤷����ȿ�͡ˤν�����Ԥ�
	  if(Shape.shape_matrix
	     [Oxygen_atom.i_cell_x][Oxygen_atom.i_cell_z] != SHAPE_SPACE )
	    {
	      flag_deposition = 
		Shape.oxidation_Si(Oxygen_atom.i_cell_x, 
				   Oxygen_atom.i_cell_z,   1.0 ) ;

	      if(flag_deposition == true )
		{
		  // cout << "Oxidized: \t" <<  Oxygen_atom.i_previous_x 
		  //  << "\t" << Oxygen_atom.i_previous_z << endl ;
		  break ;  // get out of the loop if oxidized
		}
	      //��ȿ�� : 
	      // ̵�¥롼�פ򤵤��뤿��ˡ������󥿤��Ѱդ��롣
	      // γ�Ҥ򲾤˿ʤޤ��Ƥߤơ�(pos_v_togo)
	      // �⤷�ʤ���褬������ʬ�Ǥʤ���С�������ä�
	      //     ���ľ��(50��ˤʤä���⤦����)
	      i_exception_reflect = 0;
	      Shape.get_surfacenormal(Oxygen_atom.i_cell_x, 
				      Oxygen_atom.i_cell_z) ;
	    reflection_Oxygen:
	      if( i_exception_reflect > 50 )
		break ;
	  
	      Oxygen_atom.random_reflection //(false, 0.0, 0.0 ) ;
		(true,  // ��ȿ��
		 Shape.surfacenormal_x
		 [Oxygen_atom.i_cell_x][Oxygen_atom.i_cell_z],
		 Shape.surfacenormal_z
		 [Oxygen_atom.i_cell_x][Oxygen_atom.i_cell_z]);
	  
	      Oxygen_atom.move_trans_togo(FREE_FLIGHT_PATH_Si) ;
	  
	      if( Shape.put_shape_matrix(Oxygen_atom.pos_v_togo.x, 
					 Oxygen_atom.pos_v_togo.z ,
					 &i_x, &i_z ) != SHAPE_SPACE )
		{
		  i_exception_reflect++ ;
		  goto reflection_Oxygen ;
		}
	    } 
      	} //=== while loop ��λ�ʣ�γ�Ҥ�trajectory��
    } // ���Ǥν�����λ
}




  
