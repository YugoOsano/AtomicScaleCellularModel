#ifndef _FILEIO_PARTICLE_H_DEFINED_

//#include "header_flux_model.h" 

void input_array(char file_name[], 
		 double  array1[] ,  double array2[],
		 long int  array_size );
//---------------------------------------------------
// 出力
void output_array(char file_name[], double  array1[] , 
		  long int  array_size );
void output_array(char file_name[], double  array1[] , double  array2[] ,
		  long int  array_size );
void output_array(char file_name[], double  array1[] , double  array2[] ,
		  double  array3[] , double  array4[] ,
		  long int  array_size );
// ****************************************************
// n_per_grid のような２次元配列をsplot するための
// 関数も作成
// array1[][] と書けない。下記のようにしないとエラー。


void input_array_splot(char file_name[], 
                       double **array1, 
                       int  array_size_x, int  array_size_y  );
void input_array_splot(char file_name[], 
                       int **array1, 
                       int  array_size_x, int  array_size_y  );
// ---------------------------------------------------

void output_array_splot(char file_name[], 
			double **array1 , 
			int  array_size_x, int  array_size_y,
                        double cell_size );
void output_array_splot(char file_name[], 
			int **array1 , 
			int  array_size_x, int  array_size_y,
                        double cell_size );

//--- 09/Sep/2004 HDDの節約のため x,y 座標をプロットしない関数も用意
//     （データサイズは約 1/10 になる）
void  output_array_splot(char file_name[], 
			 int **array1 , 
			 int  array_size_x, int  array_size_y ) ;


#define _FILEIO_PARTICLE_H_DEFINED_
#endif
