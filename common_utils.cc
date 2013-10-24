
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "atom_struct.h"

#include "common_utils.h"
#include "fileio_particle.h"

double rotating_angle(double x, double y )
{
  if(x > 0.0 && y >= 0.0 ) //第一象現
    {
      return atan(y/x) ;
    }
  else if(x < 0.0 ) //第二、三象現
    {
      return atan(y/x) + PI ;
    }
  else if(x > 0.0 && y < 0.0) //第四象現
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

  // --- このコードでは (3/2)πが鉛直上向き方向になる。
  return 1.5 * PI  ;
}

//==  -> 一つのセルにおけるClの数と酸素の数を同一のファイルで記録

void  output_cellinfo(int  **n_oxygen, int  **n_Cl,
		      int i_ion, char filename[], int  n_x,  int  n_z,
		      double  cell_size ) 
{
  int *tmp_matrix[N_CELL_X]; int *p_tmp_matrix; // <- tmp 変数の確保
  p_tmp_matrix = new int[N_CELL_X * N_CELL_Z] ;
  
  for(int i_x = 0; i_x < N_CELL_X ; i_x++) //== メモリ確保 ==
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

//== ファイル入力 ==
void  input_cellinfo(int  **n_oxygen, int  **n_Cl,
		     char filename[] , int  n_x,  int  n_z ) 
{
  int *tmp_matrix[N_CELL_X]; int *p_tmp_matrix; // <- tmp 変数の確保
  p_tmp_matrix = new int[N_CELL_X * N_CELL_Z] ;
  
  for(int i_x = 0; i_x < N_CELL_X ; i_x++) //== メモリ確保 ==
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
