
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "atom_struct.h"

#include "common_utils.h"
#include "fileio_particle.h"

double rotating_angle(double x, double y )
{
  if(x > 0.0 && y >= 0.0 ) //���ݸ�
    {
      return atan(y/x) ;
    }
  else if(x < 0.0 ) //���󡢻��ݸ�
    {
      return atan(y/x) + PI ;
    }
  else if(x > 0.0 && y < 0.0) //��;ݸ�
    {
      return atan(y/x) + 2.0 * PI ;
    }
  else if( x == 0.0 && y >= 0.0 )
    {
      return 0.5 * PI ;
    }
  else if( x == 0.0 && y < 0.0 )
    {
      return 1.5 * PI ;
    }

  // --- ���Υ����ɤǤ� (3/2)�Ф���ľ����������ˤʤ롣
  return 1.5 * PI  ;
}

//==  -> ��ĤΥ���ˤ�����Cl�ο��Ȼ��Ǥο���Ʊ��Υե�����ǵ�Ͽ

void  output_cellinfo(int  **n_oxygen, int  **n_Cl,
		      int i_ion, char filename[], int  n_x,  int  n_z,
		      double  cell_size ) 
{
  int *tmp_matrix[N_CELL_X]; int *p_tmp_matrix; // <- tmp �ѿ��γ���
  p_tmp_matrix = new int[N_CELL_X * N_CELL_Z] ;
  
  for(int i_x = 0; i_x < N_CELL_X ; i_x++) //== ������� ==
    tmp_matrix[i_x] = p_tmp_matrix + ( i_x * N_CELL_Z ) ;
 
  for(int i_z = 0; i_z < N_CELL_Z ; i_z++)
    {
      for(int i_x = 0; i_x < N_CELL_X; i_x++)
	{
	  tmp_matrix[i_x][i_z] = 10 * n_oxygen[i_x][i_z] + n_Cl[i_x][i_z] ;
	}
    }
  
  output_array_splot(filename ,
		     tmp_matrix ,
		     n_x,  n_z ) ;
 
  delete [] p_tmp_matrix ;
}

//== �ե��������� ==
void  input_cellinfo(int  **n_oxygen, int  **n_Cl,
		     char filename[] , int  n_x,  int  n_z ) 
{
  int *tmp_matrix[N_CELL_X]; int *p_tmp_matrix; // <- tmp �ѿ��γ���
  p_tmp_matrix = new int[N_CELL_X * N_CELL_Z] ;
  
  for(int i_x = 0; i_x < N_CELL_X ; i_x++) //== ������� ==
    tmp_matrix[i_x] = p_tmp_matrix + ( i_x * N_CELL_Z ) ;

  input_array_splot(filename ,
		    tmp_matrix ,
		    n_x,  n_z ) ;

  for(int i_z = 0; i_z < N_CELL_Z ; i_z++)
    {
      for(int i_x = 0; i_x < N_CELL_X; i_x++)
	{
	  n_Cl[i_x][i_z]     = tmp_matrix[i_x][i_z] % 10 ;
	  n_oxygen[i_x][i_z] = 
	    (tmp_matrix[i_x][i_z] - (tmp_matrix[i_x][i_z] % 10)) / 10 ;
	}
    }

  delete [] p_tmp_matrix ;
}
