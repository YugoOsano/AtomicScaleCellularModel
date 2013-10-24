// shape_trim2.cc
//  
//   shape_trim.cc ���礭���ʤä��Τ�ʬ��

#include <iostream.h>
#include <fstream.h>
#include <stddef.h>
#include <math.h>


#include "shape_trim.h"
#include "random_number.h"
#include "fileio_particle.h"
#include "common_utils.h"

//=====================================================
//== chemical etching
void Shape_trim_class::desorb_Si_chemical_etch(int i_x, int i_z)
{
  if(shape_matrix[i_x][i_z] == SHAPE_Si &&
     n_Clbond[i_x][i_z]     >= 4        &&//== ����Cl��4�ʾ�
     n_oxygen[i_x][i_z]     == 0  )  //== ���Ǥε��夬�ʤ����
    {
      shape_matrix[i_x][i_z] = SHAPE_SPACE ;
      n_Clbond[i_x][i_z]     = 0 ;
      n_oxygen[i_x][i_z]     = 0 ;

      // æΥSi �Υ�����
      cntr_desorbed_Si++ ;
      //  cout << i_x << "\t" << i_z << "\tchemical_etch!\n" ;
    }
}

//=====================================================
//==   SiClx(x = 1��4)��æΥ�δؿ�

void Shape_trim_class::
desorb_SiClx(int i_x, int i_z) 
{
  //==  SiClx ��æΥ�γ�Ψ��SiCl4�� x/4 �Ȥ���
  if(shape_matrix[i_x][i_z] == SHAPE_Si &&
     RAN0() < double(n_Clbond[i_x][i_z]) / double(NMAX_ADATOM))
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

double etch_yield_disperse_oxidation(double theta, int n_oxygen )
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
      if(n_oxygen > 6)
	n_oxygen = 6 ;

      tmp_coef = ( n_oxygen / 6.0 ) * tmp_coef * SELECTIVITY_SiO2 
	+        ( 1.0 - n_oxygen / 6.0 ) ;

      //std::cout << theta << "\tCoef: " << tmp_coef << "\n" ;
    }
  else //== poly-Si�ξ�� ==
    {
      tmp_coef = 1.0 ;
    }
  return  tmp_coef ;
}



//***************************************
// Cl density distribution - depth�����פ�Ȥ롣

void Shape_trim_class::output_Cl_density_depth(char file_name[],
					       int  time_step)
{
  //-- ���פ�Ȥ뤿��Υ����󥿡�������� 0 - CL_DEPTH_RANGE��
  int *counter_Cl ;
  counter_Cl = new int[CL_DEPTH_RANGE] ;
 
  for (int i = 0 ; i < CL_DEPTH_RANGE ; i++)
    counter_Cl[i] = 0 ; //-- �����

  bool flag_below_surface ; // ��ɽ�̤��Ⲽ�Ǥ��뤫�ɤ�����flag
  int  tmp_depth ;          // ��ɽ�̤���ο���

  for(int i_x = 0; i_x < N_CELL_X ; i_x++)
    {
      flag_below_surface = false ;
      tmp_depth          = 0     ;

      for(int i_z = 0; i_z < N_CELL_Z ; i_z++)
	{
	  if(shape_matrix[i_x][i_z] == SHAPE_Si)
	    flag_below_surface = true ;

	  if(flag_below_surface == true &&
	     tmp_depth < CL_DEPTH_RANGE )
	    {
	      counter_Cl[tmp_depth] += n_Clbond[i_x][i_z] ;
	      tmp_depth++ ;
	    }
	}
    }
  // -- �ե��������
  char OUT1[50] ; 
  sprintf( OUT1, "%s%d.dat", file_name, time_step);
  ofstream  cl_density_depth_file(OUT1) ;
  for(int i = 0 ; i < CL_DEPTH_RANGE ; i++)
    cl_density_depth_file << i << "\t" 
			  << counter_Cl[i] << endl ;

  delete [] counter_Cl ;
}

//--- �����ե�������ɤ߹���
//    ���ϡ�file name(Si, Cl)
void Shape_trim_class::
input_profile(char file_Si[],char file_Cl[]) 
{
  input_array_splot(file_Si, shape_matrix ,
		     N_CELL_X, N_CELL_Z ) ;

  input_cellinfo(n_oxygen, n_Clbond , file_Cl, 
		 N_CELL_X, N_CELL_Z ) ;
}
	  
//----------------------------
void  Shape_trim_class::remove_isolated_Si() 
{
  // -- 
  int i_x, i_z ;

  for( i_z = 2; i_z < N_CELL_Z - 2 ; i_z++ )
    {
      for( i_x = 1; i_x < N_CELL_X - 1; i_x++ )
	{
	  //-- ���Υ����Si ���Ҥ����뤫��
	  if(shape_matrix[i_x][i_z] == SHAPE_Si)
	    {
	      //-- ���ϤΣ����뤬���٤ƶ���Ǥ���С�
	      //    Si���Ҥ������
	      if(shape_matrix[i_x - 1][i_z - 1] == SHAPE_SPACE &&
		 shape_matrix[i_x - 1][i_z    ] == SHAPE_SPACE &&
		 shape_matrix[i_x - 1][i_z + 1] == SHAPE_SPACE &&
		 shape_matrix[i_x    ][i_z - 1] == SHAPE_SPACE &&
		 shape_matrix[i_x    ][i_z + 1] == SHAPE_SPACE &&
		 shape_matrix[i_x + 1][i_z - 1] == SHAPE_SPACE &&
		 shape_matrix[i_x + 1][i_z    ] == SHAPE_SPACE &&
		 shape_matrix[i_x + 1][i_z + 1] == SHAPE_SPACE )
		{
		  shape_matrix[i_x][i_z] = SHAPE_SPACE ;
		  n_Clbond[i_x][i_z]     = 0 ;
		  cout << i_x << "\t" << i_z << "\n" ;
		}
	    }
	}
    }
}
