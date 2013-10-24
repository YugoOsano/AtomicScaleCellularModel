#ifndef _FILEIO_PARTICLE_H_DEFINED_

//#include "header_flux_model.h" 

void input_array(char file_name[], 
		 double  array1[] ,  double array2[],
		 long int  array_size );
//---------------------------------------------------
// ����
void output_array(char file_name[], double  array1[] , 
		  long int  array_size );
void output_array(char file_name[], double  array1[] , double  array2[] ,
		  long int  array_size );
void output_array(char file_name[], double  array1[] , double  array2[] ,
		  double  array3[] , double  array4[] ,
		  long int  array_size );
// ****************************************************
// n_per_grid �Τ褦�ʣ����������splot ���뤿���
// �ؿ������
// array1[][] �Ƚ񤱤ʤ��������Τ褦�ˤ��ʤ��ȥ��顼��


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

//--- 09/Sep/2004 HDD������Τ��� x,y ��ɸ��ץ�åȤ��ʤ��ؿ����Ѱ�
//     �ʥǡ������������� 1/10 �ˤʤ��
void  output_array_splot(char file_name[], 
			 int **array1 , 
			 int  array_size_x, int  array_size_y ) ;


#define _FILEIO_PARTICLE_H_DEFINED_
#endif
