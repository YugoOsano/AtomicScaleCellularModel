
// etch_products.cc
// 2004/07/29

#include <iostream>
#include <math.h>
#include "atom_struct.h"

#include "etch_products.h"


Etch_product_class::
Etch_product_class(double mass_input,  int flag_boundary_input,
		   bool   flag_open_space, 
		   int    n_chlorination_input ,
		   int    n_oxidation_input ,
		   double reincident_probability_whole,
		   double sticking_probability_input  ) 
  : Neutral_class(mass_input, flag_boundary_input)
{
  n_chlorination = n_chlorination_input  ;
  n_oxidation    = n_oxidation_input  ; 
  
  sticking_probability = sticking_probability_input ;
  
  if(flag_open_space == true)
    reincidence_probability = reincident_probability_whole
      *           ( 0.5 / ((double(N_CELL_X) - 184.0)/double(N_CELL_X)) );
  else
    reincidence_probability = reincident_probability_whole
      *           ( 0.5 / ((double(N_CELL_X) - 368.0)/double(N_CELL_X)) );

  std::cout << "effective Pr: " << reincidence_probability  << "\n" ;
}
Etch_product_class::~Etch_product_class()
{
}

void Etch_product_class::get_previous_cell()
{
  //(1) v_x > 0 
  if(pos_v.v_x > 0.0 )
    {
      //-- ��z - (v_z/v_x)��x(0����) < 0 �ʤ�� 
      //      ���z����������ˤ��������
      if( dz_cell - ( pos_v.v_z / pos_v.v_x ) * dx_cell < 0.0 )
	{
	  i_previous_x = i_cell_x ;
	  i_previous_z = i_cell_z - 1 ;
	}
      //-- ��z - (v_z/v_x)��x(0����) > cell size �ʤ�� 
      //      ����z�����������ˤ��������
      else if(dz_cell - ( pos_v.v_z / pos_v.v_x ) * dx_cell >
	      CELL_SIZE)
	{
	  i_previous_x = i_cell_x ;
	  i_previous_z = i_cell_z + 1 ;
	}
      else // ����x����������ˤ��������
	{
	  i_previous_x = i_cell_x - 1 ;
	  i_previous_z = i_cell_z  ;
	}
    }
  // v_x < 0 �Ǥϡ�x = CELL_SIZE ������
  // ��z + (v_z/v_x)( CELL_SIZE - ��x)
  else if(pos_v.v_x < 0.0) 
    {
      //-- (����) < 0 �ʤ�� 
      //      ���z����������ˤ��������
      if( dz_cell + 
	  ( pos_v.v_z / pos_v.v_x ) * (CELL_SIZE - dx_cell) < 0.0 )
	{
	  i_previous_x = i_cell_x ;
	  i_previous_z = i_cell_z - 1 ;
	}
      //-- (����) > cell size �ʤ�� 
      //      ����z�����������ˤ��������
      else if(dz_cell +
	      ( pos_v.v_z / pos_v.v_x ) * (CELL_SIZE - dx_cell) >
	      CELL_SIZE)
	{
	  i_previous_x = i_cell_x ;
	  i_previous_z = i_cell_z + 1 ;
	}
       else // ����x�����������ˤ��������
	{
	  i_previous_x = i_cell_x + 1 ;
	  i_previous_z = i_cell_z  ;
	}


    }

  //=== �㳰������domain ����ФƤ��ʤ�������
  if     ( i_previous_x < 0 )
    i_previous_x = 0 ;
  else if( i_previous_x > N_CELL_X - 1)
    i_previous_x = N_CELL_X - 1 ;
  
  if     ( i_previous_z < 0 )
    i_previous_z = 0 ;
  else if( i_previous_z > N_CELL_Z - 1)
    i_previous_z = N_CELL_Z - 1 ;
	
}

//========================================
// Siɽ�̤�γ�Ҥ����֤���®��Ϳ����
bool Etch_product_class::
allocate_on_surface(Shape_trim_class *Shape_pointer,
		    Particle_location_velocity_struct  pos_v_input,
		    double v )
{
  // random_reflection �ؿ���ȿ�ͤ��Ƥ�Ʊ��®�٤�
  // �ݤĤΤǺǽ��Ŭ�ڤ�®�٤�Ϳ���� �ʸ�ǰ��֤Τ����֡�

  inject_from_top( v ) ;

  pos_v.x = pos_v_input.x ;
  pos_v.y = pos_v_input.y ;
  pos_v.z = pos_v_input.z ;

  // -- γ�ҤΤ������phase�����
  //�ʥ����ֹ� i_cell_x, i_cell_z�������
  get_position() ;

  //ɽ�̤���γȻ� : 
  // ̵�¥롼�פ򤵤��뤿��ˡ������󥿤��Ѱդ��롣
  // γ�Ҥ򲾤˿ʤޤ��Ƥߤơ�(pos_v_togo)
  // �⤷�ʤ���褬������ʬ�Ǥʤ���С�������ä�
  //     ���ľ��(50��ˤʤä���⤦����)
 
  int i_x, i_z ;
 
  int i_exception_reflect = 0;
  Shape_pointer->get_surfacenormal(i_cell_x, 
				   i_cell_z) ;
 
 diffusion_SiClx:
  if( i_exception_reflect > 50 )
    return false ;
  
  random_reflection 
    (true,  // ��ȿ��
     Shape_pointer->surfacenormal_x[i_cell_x][i_cell_z],
     Shape_pointer->surfacenormal_z[i_cell_x][i_cell_z]) ;
  
  move_trans_togo(FREE_FLIGHT_PATH_Si) ;
  
  if( Shape_pointer->put_shape_matrix(pos_v_togo.x, 
				      pos_v_togo.z ,
				      &i_x, &i_z ) != SHAPE_SPACE )
    {
      i_exception_reflect++ ;
      goto diffusion_SiClx ;
    }
  return true ;

}

//=======================================================
// γ�Ҥ����Ͱʹߤν�����ؿ�������
void Etch_product_class::
whole_flight(Shape_trim_class *Shape_pointer,
	     bool flag_inject_from_side) 
{
  bool flag_deposition    ;
  int  i_exception_reflect ;
  int  i_x, i_z ;

  i_exception_move = 0 ;

  //== ��������������� : neutral_particle.cc ��all_process�ؿ��Ȥۤ�Ʊ��
  if(flag_inject_from_side == true)
    {
      inject_from_right_side(pos_v.v_r);
      //== ���ͤΰ��֤����֤�̵����� i_exception_move ���礭��
      //   �ͤ��������� while loop ��skip ����
      if(Shape_pointer->put_shape_matrix(pos_v.x - 0.5 * CELL_SIZE, pos_v.z ,
					 &i_x, &i_z ) != SHAPE_SPACE )
	i_exception_move = N_CELL_Z * 200 + 1 ;
    }
  
  while(1) 
    {
      //--- �ΰ������ˤ��ޤ�Ĺ���ȤɤޤäƤ����硢
      //    �㳰������γ�Ҥ�õ��
      i_exception_move++ ;
      if( i_exception_move > N_CELL_Z * 200 )
	break ;
      
      // -- �¿ʱ�ư
      move_trans(FREE_FLIGHT_PATH_Si) ;
      
      if (flag_inside == false) // -- ���˽Ф���
	break ;
      
      // -- γ�ҤΤ������phase�����
      get_position() ;
	      
      // -- ����Ǥ���С�������Υ����ֹ����� -> deposit
      if(Shape_pointer->shape_matrix
	 [i_cell_x][i_cell_z] != SHAPE_SPACE )
	{
	  get_previous_cell() ; // �����ֹ����
	  flag_deposition = 
	    Shape_pointer->deposit_Si(i_previous_x,     // deposit
				      i_previous_z,
				      i_cell_x,    i_cell_z , 
				      n_chlorination, n_oxidation,
				      sticking_probability ) ;	  
	  if(flag_deposition == true )
	  {
	    // cout << "SiClx deposited: \t" << i_previous_x 
	     //<< "\t" << i_previous_z << endl ;
	    break ;  // get out of the loop if deposited
	  }
	  //��ȿ�� : 
	  // ̵�¥롼�פ򤵤��뤿��ˡ������󥿤��Ѱդ��롣
	  // γ�Ҥ򲾤˿ʤޤ��Ƥߤơ�(pos_v_togo)
	  // �⤷�ʤ���褬������ʬ�Ǥʤ���С�������ä�
	  //     ���ľ��(50��ˤʤä���⤦����)
	  i_exception_reflect = 0;
	  Shape_pointer->get_surfacenormal(i_cell_x, 
					   i_cell_z) ;
	reflection_SiClx:
	  if( i_exception_reflect > 50 )
	    break ;
	  
	  random_reflection //(FALSE, 0.0, 0.0 ) ;
	    (true,  // ��ȿ��
	     Shape_pointer->surfacenormal_x[i_cell_x][i_cell_z],
	     Shape_pointer->surfacenormal_z[i_cell_x][i_cell_z]) ;
	  
	  move_trans_togo(FREE_FLIGHT_PATH_Si) ;
	  
	  if( Shape_pointer->put_shape_matrix(pos_v_togo.x, 
					      pos_v_togo.z ,
					      &i_x, &i_z ) != SHAPE_SPACE )
	    {
	      i_exception_reflect++ ;
	      goto reflection_SiClx ;
	    }
	} 
    } //=== while loop ��λ�ʣ�γ�Ҥ�trajectory��
}
