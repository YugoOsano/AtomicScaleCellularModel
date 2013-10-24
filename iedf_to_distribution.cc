// iedf_to_distribution.cc

// Nomura ��IEDF���ϤΥե����ޥå�
/* ���ͥ륮��(eV)    ��Ψ(a.u.) 
0.000000        0.000000e+000
0.100000        3.894278e-006
0.200000        1.559246e-005
0.300000        2.078770e-005
0.400000        7.351976e-006
0.500000        0.000000e+000
0.600000        9.140265e-006
0.700000        2.766339e-006

IADF���ϤΥե����ޥå�(iadf_angle_to_energy.dat)

����(deg)  ��Ψ(a.u.)
-90     2
-89     0
-88     0
 .....

-2      2117
-1      4844
0       8318
1       10585
2       8869

......

87      1
88      3
89      1
90      0




���Ȥˡ�ʬ�۴ؿ� ���(x)dx ���ᡢ����Ȥ���Ϳ���롣
��������򡢵���줿ʬ�۴ؿ��εմؿ���Ϳ���뤳�Ȥ�
IEDF �˽��ä�����������롣
http://hi.sakura.ne.jp/~nmaeda/web/fran.shtml
�򻲾�
 
����Ϸ����֤���������ɬ�פ�����Τǡ�
1  IEDF�Υե������ɤ߹��� -> distribution �κ���

2  distribution ���б��������ʬ�ۤκ���

���̤δؿ��ˤ���

1  iedf_to_distribution(char[] �ե�����̾ , int ����Υ�����,
                       -> ���� double[] ���ͥ륮��  ,double[] ʬ��)
 
2  randomized_distribution   (long int ����μ�, int ����������ο�,
                              double[] ���ͥ륮��  ,double[] ʬ��,
                     -> ����  double[] (�������), double[] ����줿���) 
*/
#ifndef _IEDF_TO_DISTRIBUTION_
 
#include<iostream.h>
#include<fstream.h>
#include<string.h>
//#include "common_const.h"
#include "fileio_particle.h"
#include "random_number.h" 

#include "iedf_to_distribution.h"

//void input_iedf_array(char[],double[], double[],long int ) ;

//void iedf_to_distribution(char[] , int ,
//			  double[]   ,double[] ) ;

//void randomized_distribution (long int , int ,
//			      double[]   ,double[],
//			      double[], double[] ) ; 

// test main
// ��������ϴؿ���������˰����ˤ���
/*
const long int ARRAY_SIZE = 2502 ;
long int RANDOM_SEED = 2 ;
const int N_CREATING_RANDOM  = 1000 ;
*/

void iedf_to_distribution(char file_name[], int array_size,
			  double energy_in_eV[] , 
			  double distribution[]  	)

{
  // ---- ���� distribution ��index �� 0 ���� array_size - 1 
  // ---- �Ǥ��뤳�Ȥ���� 14 Feb 2003

  // -- ���󥵥�����variable �Ǥ�����ˤϡ�ưŪ�˳��ݤ���
  
  double *amount_function;
  amount_function = new double[array_size + 1] ;


  input_array(file_name , energy_in_eV, amount_function, array_size) ;

  distribution[0] = amount_function[0] ;

  for(int i = 0; i < array_size - 1 ; i++) // <- ���[i+1] ��index��ȤäƤ��롣
    {
      //cout << energy_in_eV[i] << "\t" << amount_function[i] << 
      //	"\t" << distribution[i] << endl;

      distribution[i+1] = distribution[i] + amount_function[i+1] ;
    }

  // ��Ψʬ�ۤˤʤ�褦����������Ԥ�
  // (distribution �κǸ���ͤ� 1 �Ȥʤ�褦�ˤ���)

  double normalizer = 1.00 / distribution[array_size - 1] ;

  for(int i = 0; i < array_size - 1 ; i++)
    {
      distribution[i] *= normalizer ;
    
      //cout << energy_in_eV[i] << "\t" << amount_function[i] << 
      //"\t" << distribution[i] << endl;
    }
  distribution[array_size - 1] = 1.00 ; // �Ǹ��1�ˤ��Ƥ���
  
  delete [] amount_function ;
      
} 

//-------------------------------
// ���ȯ��(0-1) -> distribution ����������� -> 
// ������б����� energy_in_eV �����

void randomized_distribution(long int random_seed, int n_creating_random,
			     double energy_in_eV[]   ,double energy_threshold,
                             double distribution[] ,
			     double random_sample[] , 
			     double randomized_distribution_array[] ) 
{
  for (int i_create = 0 ; i_create < n_creating_random ; i_create++)
    {
      // == introduce the threshold ==
      randomized_distribution_array[i_create] = 0.0 ;
      while (randomized_distribution_array[i_create] <= energy_threshold)
        {
          random_sample[i_create] = RAN0() ;
          
          int i = 0;
          while( distribution[i+1] < random_sample[i_create])
            {
              i += 1 ;
            }
          // distribution[i+1] > random_sample[i_create]
          // �Ȥʤä��顢 i �����Ǥ������
          
          // cout << random_sample << "\t" << energy_in_eV[i] << endl ;
          randomized_distribution_array[i_create] = energy_in_eV[i] ;
        }
    }
}
//===========================================
double distribution_from_energy(double input_energy_eV,
                                double energy_in_eV[] , 
                                double distribution[]  ,
				int array_size )
{
  int i = 0;
      
  while( energy_in_eV[i+1] < input_energy_eV &&
	 i < array_size - 1 )
	{
	  i++ ;
	}
  cout << input_energy_eV << "\t" << energy_in_eV[i] << "\t"
       << distribution[i] << endl ;
  return distribution[i] ;

}
//----------------------------


#define _IEDF_TO_DISTRIBUTION_ 
#endif

/*

void input_iedf_array(char file_name[], 
		      double  array1[] ,  double array2[],
		      long int  array_size )
{
  ifstream  data_file(file_name) ;
  for (int i = 0 ; i < array_size ; i++ )
    {
      data_file >> array1[i] >> array2[i];
    }
 
} 
*/
