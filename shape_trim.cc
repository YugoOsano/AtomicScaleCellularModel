// shape_trim.cc

#include <iostream>
#include <fstream.h>
#include <stddef.h>
#include <math.h>


#include "shape_trim.h"
#include "random_number.h"

//-- constructor

Shape_trim_class::Shape_trim_class()
{}

Shape_trim_class::Shape_trim_class(bool flag_open_space)
{
  //-- �������
  p_shape_matrix = new int[N_CELL_X * N_CELL_Z] ;
  p_n_Clbond     = new int[N_CELL_X * N_CELL_Z] ;
  p_n_oxygen     = new int[N_CELL_X * N_CELL_Z] ;

  //-- �����󥿥ꥻ�å�
  cntr_desorbed_Si   = 0 ;
  cntr_deposited_Si  = 0 ;
  cntr_desorbed_mask = 0 ;
  flag_get_bottom    = false ;

  for(int i_x = 0; i_x < N_CELL_X ; i_x++)
    {
      shape_matrix[i_x] = p_shape_matrix + ( i_x * N_CELL_Z ) ;
      n_Clbond[i_x]     = p_n_Clbond     + ( i_x * N_CELL_Z ) ;
		n_oxygen[i_x]     = p_n_oxygen + ( i_x * N_CELL_Z ) ;
    }


  //--- �ƥ���ξ��֤����Ƥ���
  for(int i_z = 0; i_z < N_CELL_Z ; i_z++)
    {
      if(i_z < N_TOP_Si ) // marginal space beyond the surface 
	{
	  for(int i_x = 0; i_x < N_CELL_X ; i_x++)
            {
	      // -- �ޥ�����ʬ
	      if( i_z >= N_TOP_MASK &&
		  (i_x <  N_LEFT_MASK ||
		   i_x >= N_RIGHT_MASK ) )
		{
		  shape_matrix[i_x][i_z] = HARD_MASK ;
		  //== �ޥ���¦�ɤ˷��гѤ������� ==
		  if( i_x >= N_LEFT_MASK - ((N_TOP_Si - i_z) * tan(SLOPE_ANGLE_SIDEWALL)) &&
		      i_x <  N_RIGHT_MASK+ ((N_TOP_Si - i_z) * tan(SLOPE_ANGLE_SIDEWALL)) )
		    shape_matrix[i_x][i_z] = SHAPE_SPACE ;
		  
		  //== ��¦ open space �ξ�硧��¦�Υޥ����Ϥʤ�
		  if( flag_open_space == true &&
		      i_x >= N_RIGHT_MASK  )
		    shape_matrix[i_x][i_z] = SHAPE_SPACE ;
		}
	      else
		shape_matrix[i_x][i_z] = SHAPE_SPACE ;

#ifdef       _NO_MASK_
	      shape_matrix[i_x][i_z] = SHAPE_SPACE ;
#endif
            }
        }
      else if(i_z < N_BOTTOM_Si)// �ѥ��������̤�겼��
	{ 
	  for(int i_x = 0; i_x < N_CELL_X ; i_x++)
	    shape_matrix[i_x][i_z] = SHAPE_Si ;
	}
      else
	{
	  for(int i_x = 0; i_x < N_CELL_X ; i_x++)
	    shape_matrix[i_x][i_z] = HARD_MASK ;
	}
    }
  for(int i_z = 0; i_z < N_CELL_Z ; i_z++)
    {
      for(int i_x = 0; i_x < N_CELL_X ; i_x++)
	{
	  n_Clbond[i_x][i_z] = 0 ;
	  n_oxygen[i_x][i_z] = 0 ;
	}
    }
  // -- ������� �����ͳѰ�¸��

  p_surfacenormal_x = new double[N_CELL_X * N_CELL_Z] ;
  p_surfacenormal_z = new double[N_CELL_X * N_CELL_Z] ;
  for(int i_x = 0; i_x < N_CELL_X ; i_x++)
    {
      surfacenormal_x[i_x] = p_surfacenormal_x + ( i_x * N_CELL_Z ) ;
      surfacenormal_z[i_x] = p_surfacenormal_z + ( i_x * N_CELL_Z ) ;
    }
  get_surfacenormal_all() ;



  //______________________________
  first_ptr    = NULL ;
  n_sticking   = 0 ;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
//---------------------------------
Shape_trim_class::~Shape_trim_class()
{
  delete [] p_shape_matrix ;
  delete [] p_n_Clbond   ;
  delete [] p_n_oxygen  ;


  delete [] p_surfacenormal_x ;
  delete [] p_surfacenormal_z ;

}

//----------------------------
int Shape_trim_class::put_shape_matrix(double x, double z , 
				       int *i_x, int *i_z )
{
  if(x >= 0.0   &&    x < SYSTEM_WIDTH_X  &&
     z >= 0.0   &&    z < SYSTEM_HEIGHT_Z  )
    {
      //cout << shape_matrix[ int( x / CELL_SIZE ) ][ int( z / CELL_SIZE ) ]
      //   << endl ;

      *i_x = int( x / CELL_SIZE ) ;
      *i_z = int( z / CELL_SIZE ) ;

      return shape_matrix[*i_x][*i_z] ;
    }
  else
    {
      return 0 ;
    }
}
 
//===============================================
int  Shape_trim_class::
find_solid_nearest_neighbor(double x, double z,  
			    double v_x, double v_z,
			    int *i_x, int *i_z  ,
			    int *i_x_particle, int *i_z_particle) 
{
  //---- �ɲá�γ�Ҥΰ���(x, z)�Υ����index���֤�
  put_shape_matrix(x, z , i_x_particle, i_z_particle );

  // (1)�ޤ�������
  if(put_shape_matrix(x, z , i_x, i_z ) != SHAPE_SPACE)
    return shape_matrix[*i_x][*i_z];

  else if(put_shape_matrix(x, z , i_x, i_z ) == SHAPE_SPACE &&
	  z >= 0.0 + CELL_SIZE  &&  z < SYSTEM_HEIGHT_Z - CELL_SIZE)
    {
      //==(2)�������Τ�����®�٥٥��ȥ�������˶ᤤ��
      //   ®��(v_x, v_z) ������˹�碌�ơ�-1, 1 �Ȥʤ�褦��
      //   �ѿ����Ѱ�
      int  flag_x = 0, flag_z = 0;
      if(v_x != 0.0)
	flag_x = int(v_x / fabs(v_x)) ;
      if(v_z != 0.0)
	flag_z = int(v_z / fabs(v_z)) ;

      //== x�����μ������������θ
      double x_left, x_right ;
      if(x < 0.0 + CELL_SIZE )
	{
	  x_left  = x - flag_x * CELL_SIZE + SYSTEM_WIDTH_X;
	  x_right = x + flag_x * CELL_SIZE ;
	}
      else if(x >= SYSTEM_WIDTH_X  - CELL_SIZE)
	{
	  x_left  = x - flag_x * CELL_SIZE ;
	  x_right = x + flag_x * CELL_SIZE - SYSTEM_WIDTH_X;
	}
      else 
	{
	  x_left  = x - flag_x * CELL_SIZE ;
	  x_right = x + flag_x * CELL_SIZE ;
	}

      // x������z�����Τɤ������ʬ���礭����Ƚ�ꤷ������˱�����
      // ��˽����򤹤�
      if(fabs(v_x) >= fabs(v_z))
	{
	  if(put_shape_matrix(x_right, z, 
			      i_x, i_z ) != SHAPE_SPACE)
	    return shape_matrix[*i_x][*i_z];
	  if(put_shape_matrix(x , z + flag_z * CELL_SIZE, 
			      i_x, i_z ) != SHAPE_SPACE)
	    return shape_matrix[*i_x][*i_z];
	  if(put_shape_matrix(x , z - flag_z * CELL_SIZE, 
			      i_x, i_z ) != SHAPE_SPACE)
	    return shape_matrix[*i_x][*i_z];
	  if(put_shape_matrix(x_left , z, 
			      i_x, i_z ) != SHAPE_SPACE)
	    return shape_matrix[*i_x][*i_z];
	}
      else 
	{
	  if(put_shape_matrix(x , z + flag_z * CELL_SIZE, 
			      i_x, i_z ) != SHAPE_SPACE)
	    return shape_matrix[*i_x][*i_z];
	  if(put_shape_matrix(x_right , z, 
			      i_x, i_z ) != SHAPE_SPACE)
	    return shape_matrix[*i_x][*i_z];
	  if(put_shape_matrix(x_left  , z, 
			      i_x, i_z ) != SHAPE_SPACE)
	    return shape_matrix[*i_x][*i_z];
	  if(put_shape_matrix(x , z - flag_z * CELL_SIZE, 
			      i_x, i_z ) != SHAPE_SPACE)
	    return shape_matrix[*i_x][*i_z];
	}
    }
  return SHAPE_SPACE ;
}

// --- ���塧 (bare) Si ��Cl bond �������䤹
bool Shape_trim_class::
settle_Cl_into_bareSi(int i_x, int i_z, 
		      double adsorption_probability) 
{
  // -- ����ϣ��Ĥޤ� / ��hard mask �Ǥʤ���Si �Ǥ��뤫��Ƚ��
  //if( n_Clbond[i_x][i_z] < NMAX_ADATOM &&
  //  adsorption_probability > RAN0()  ) 
  // oxidation ���θ���ƽ񤭴�����ȿ���������Ǥ򲡤��ऱ�Ƶ��夹�뤳�ȤϤʤ���
  if( n_Clbond[i_x][i_z] + n_oxygen[i_x][i_z] * 2  
      < NMAX_ADATOM                        && 
      shape_matrix[i_x][i_z] == SHAPE_Si   &&
      adsorption_probability > RAN0()  )
    {
      n_Clbond[i_x][i_z]++ ;
      return true ;
    }
  //== ������γ��Cl��flux count �򤹤����TRUE ===
  return false;//TRUE;//
}

// --- ���塧 Si ��Cl bond �������䤷��NMAX_ADATOM + 1 �ˤʤä���æΥ
//void Shape_trim_class::settle_Cl_and_desorption(int i_x, int i_z) 
//{
//  n_Clbond[i_x][i_z]++ ;
//  if(n_Clbond[i_x][i_z] >= NMAX_ADATOM + 1 )
//    {
//      shape_matrix[i_x][i_z] = SHAPE_SPACE ;
//      n_Clbond[i_x][i_z]     = 0 ;
//    }     
//}
//******** 02 Dec 2003 ******************
// --- have a Cl atom adsorbed
// --- ��Ͽ�Ѥ˰���(x,z)�������ͤˤ���褦�ʴؿ������
// ---> (shape_stickingCl.cc)
void Shape_trim_class::adsorb_Cl(int i_x, int i_z) 
{
  if( shape_matrix[i_x][i_z] == SHAPE_Si)
    n_Clbond[i_x][i_z]++ ;
}
// --- desorb Si atom if adatoms are full and 
//                 when a Cl ion impacts thereon.
void Shape_trim_class::desorb_Si(int i_x, int i_z) 
{
  if(shape_matrix[i_x][i_z] == SHAPE_Si &&
     n_Clbond[i_x][i_z] + 2 * n_oxygen[i_x][i_z] >= NMAX_ADATOM)
    {
      shape_matrix[i_x][i_z] = SHAPE_SPACE ;
      n_Clbond[i_x][i_z]     = 0 ;
      n_oxygen[i_x][i_z]     = 0 ;

      // æΥSi �Υ�����
      cntr_desorbed_Si++ ;
      
      if(FLAG_INCIDENT_ANGLE == true)
	get_surfacenormal_around(i_x, i_z) ;

    }
  // -- domain �β����飲���ܤޤ���ã�����顢flag_get_bottom ��
  //    TRUE �ˤ���ʰ��ֲ��� hard mask �ˤ����
  if( i_z >= N_CELL_Z - 1)
    flag_get_bottom  =  true  ;
}

//== �ϡ��ɥޥ����ο��� ==
void Shape_trim_class::desorb_mask(int i_x, int i_z) 
{
  if(shape_matrix[i_x][i_z] == HARD_MASK)
    {
      shape_matrix[i_x][i_z] = SHAPE_SPACE ;
      cntr_desorbed_mask++ ;
      
      if(FLAG_INCIDENT_ANGLE == true)
	get_surfacenormal_around(i_x, i_z) ;

    }
}

//***********************************
// ***  05 Aug 2004 : deposition 
bool Shape_trim_class::
deposit_Si(int i_deposit_x, int i_deposit_z, 
	   int  i_neighbor_x, int i_neighbor_z, 
	   int n_Cl,  int n_O,
	   double deposition_probability)
{
  // -- ����ϣ��Ĥޤ� / ���ܥ��뤬Si or �ϡ��ɥޥ����Ǥ��뤫��
  // --  i_z = 0,1 �ˤ��������Ѥ�̵����ΤȤ���
  if( deposition_probability > RAN0() &&
      i_deposit_z > 1  &&
      shape_matrix[i_deposit_x][i_deposit_z]    == SHAPE_SPACE &&
      (shape_matrix[i_neighbor_x][i_neighbor_z] == SHAPE_Si  ||
       shape_matrix[i_neighbor_x][i_neighbor_z] == HARD_MASK ))
    {
      shape_matrix[i_deposit_x][i_deposit_z] = SHAPE_Si ;
      n_Clbond[i_deposit_x][i_deposit_z]     = n_Cl ;
      n_oxygen[i_deposit_x][i_deposit_z]     = n_O ;

      cntr_deposited_Si++ ;

      return true ;
    }
    return false ;
}

//***********************************
// ***  24 Aug 2004 : oxidation 
bool Shape_trim_class::oxidation_Si(int i_x, int i_z, 
				   double reaction_probability ) 
{
	// -- ���뤬Si �Ǥ��뤫��
  if(  reaction_probability  > RAN0() && //(1)
       n_oxygen[i_x][i_z] < 2 && // (2)
		 shape_matrix[i_x][i_z] == SHAPE_Si  )
    { // n_oxygen < 2 �Ǥ���С�ȿ������(n_oxygen++)
      
      n_oxygen[i_x][i_z]++ ;
      
      if(n_oxygen[i_x][i_z] == 1 &&  // header file (3)
	 n_Clbond[i_x][i_z] > 2  )
	n_Clbond[i_x][i_z] = 2 ;
      
      if(n_oxygen[i_x][i_z] == 2)  // header file (4)
	n_Clbond[i_x][i_z] = 0 ;
      
      return true ;
    }
  return false ;
}





// Four point calculation �ˤ�ä�ɽ�̤�ˡ����������ꤹ�롣
void  Shape_trim_class::get_surfacenormal(int i_x, int i_z)
{
  double tmp_nx = 0.0 ;//--- ˡ���٥��ȥ�
  double tmp_nz = 0.0 ;
 
  //-- ������ֹ��Ƚ����ΰ�ζ��ˤ��륻��Ǥ���� 0,0 �Ȥ����
  if( i_z > 0  &&  i_z < N_CELL_Z - 1   )
    {
      //--�岼�Υ����Ƚ��
      if(shape_matrix[i_x    ][i_z - 1] != SHAPE_SPACE) //-- ���֤Ǥʤ����
	tmp_nz  +=  - 1.0 ;
      if(shape_matrix[i_x    ][i_z + 1] != SHAPE_SPACE)
	tmp_nz  +=  + 1.0 ;
      
      //--�����Υ����Ƚ��ʼ������������θ��
      if( i_x > 0  &&  i_x < N_CELL_X - 1)
	{
	  if(shape_matrix[i_x - 1][i_z] != SHAPE_SPACE)
	    tmp_nx  +=  - 1.0 ;
	  if(shape_matrix[i_x + 1][i_z] != SHAPE_SPACE)
	    tmp_nx  +=  + 1.0 ;
	}
      else if( i_x == 0 )
	{
	  if(shape_matrix[N_CELL_X - 1][i_z] != SHAPE_SPACE)
	    tmp_nx  +=  - 1.0 ;
	  if(shape_matrix[i_x + 1][i_z] != SHAPE_SPACE)
	    tmp_nx  +=  + 1.0 ;
	}
      else if( i_x == N_CELL_X - 1 )
	{
	  if(shape_matrix[i_x - 1][i_z] != SHAPE_SPACE)
	    tmp_nx  +=  - 1.0 ;
	  if(shape_matrix[0][i_z] != SHAPE_SPACE)
	    tmp_nx  +=  + 1.0 ;
	}
    }
  surfacenormal_x[i_x][i_z] = tmp_nx ;
  surfacenormal_z[i_x][i_z] = tmp_nz ;
}
//----
void  Shape_trim_class::get_surfacenormal_around(int i_x, int i_z )
{
  if(i_z > 0  &&  i_z < N_CELL_Z - 1   )
    {
      get_surfacenormal(i_x    , i_z - 1) ;
      get_surfacenormal(i_x    , i_z + 1) ;

      if( i_x > 0  &&  i_x < N_CELL_X - 1  )
	{
	  get_surfacenormal(i_x - 1, i_z    ) ;
	  get_surfacenormal(i_x + 1, i_z    ) ;
	}
      else if(i_x == 0)
	{
	  get_surfacenormal(N_CELL_X - 1, i_z) ;
	  get_surfacenormal(i_x + 1, i_z) ;
	}
      else if(i_x == N_CELL_X - 1)
	{
	  get_surfacenormal(i_x - 1, i_z) ;
	  get_surfacenormal(      0, i_z) ;
	}
    }
}

void Shape_trim_class::get_surfacenormal_all()
{
  for(int i_z = 1 ; i_z < N_CELL_Z - 1 ; i_z++)
    {
      for(int i_x = 0 ; i_x <= N_CELL_X - 1 ; i_x++)
	{
	  get_surfacenormal(i_x , i_z) ;
	}
    }
}

// === �������®����ʬ���顢���ͳѤ����(rad)
double Shape_trim_class::get_incident_angle(int i_x, int i_z,
					    double v_x, double v_z) 
{
  //== ˡ���٥��ȥ�ʤ⤷����®�١ˤ� 0,0 �ξ�硢0 rad ���֤�
  if( surfacenormal_x[i_x][i_z] * surfacenormal_x[i_x][i_z] +
      surfacenormal_z[i_x][i_z] * surfacenormal_z[i_x][i_z] <= 0.0 ||
      v_x * v_x + v_z * v_z <= 0.0   )
    {
      return 0.0 ;
    }
  //== ���Ѥμ� |n||v| cos�� = n_x v_x + n_z v_z
  //  ���Ȥ����
  double tmp_cos = ( surfacenormal_x[i_x][i_z] * v_x + 
		     surfacenormal_z[i_x][i_z] * v_z ) 
    / ( sqrt(surfacenormal_x[i_x][i_z] * surfacenormal_x[i_x][i_z] +
	     surfacenormal_z[i_x][i_z] * surfacenormal_z[i_x][i_z]) *
	sqrt( v_x * v_x + v_z * v_z ) ) ; 

  return acos(tmp_cos) ;

}
  
//=== Y(��) : etching yield�γ��ٰ�¸��
double etch_yield_angle(double theta, int n_oxygen ) 
{
  double tmp_coef = 0.0; 

  //===������ξ�硧Mahorowala, Sawin JVST 20, 1064 (2002)
  //   paper �ǤϺ����ͤ򸵤�����������Ƥ��뤬��
  //   �����Ǥ� ��=0 ����ˤ��뤿�ᡢ0.24 ����
  if(n_oxygen >= 1)
    {
      tmp_coef = 
	0.4 * (18.7  * cos(theta)          -  64.7 * pow(cos(theta),2.0) +
	       145.2 * pow(cos(theta),3.0) - 206.0 * pow(cos(theta),4.0) + 
	       147.3 * pow(cos(theta),5.0) -  39.9 * pow(cos(theta),6.0))/0.24;

      if(tmp_coef < 0.0)
	tmp_coef = 0.0 ;
      //== ���Ǥο��ˤ����� ==
      if( n_oxygen == 1 )
	tmp_coef *= SELECTIVITY_SiO ;
      else 
	tmp_coef *= SELECTIVITY_SiO2 ;

      //std::cout << theta << "\tCoef: " << tmp_coef << "\n" ;
    }
  else //== poly-Si�ξ�� ==
    {
      if(fabs(theta) < PI / 4.0 )
	{
	  tmp_coef = 1.0 ;
	}
      else if (fabs(theta) < PI / 2.0 )
	{
	  tmp_coef =  cos(theta) /  0.707106781186548 ;
	}
      else
	{
	  tmp_coef = 0.0 ;
	}
    }
  return  tmp_coef ;
}


