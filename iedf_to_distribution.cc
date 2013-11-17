// iedf_to_distribution.cc

// Nomura のIEDF出力のフォーマット
/* エネルギー(eV)    確率(a.u.) 
0.000000        0.000000e+000
0.100000        3.894278e-006
0.200000        1.559246e-005
0.300000        2.078770e-005
0.400000        7.351976e-006
0.500000        0.000000e+000
0.600000        9.140265e-006
0.700000        2.766339e-006

IADF出力のフォーマット(iadf_angle_to_energy.dat)

角度(deg)  確率(a.u.)
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




をもとに、分布関数 ∫ρ(x)dx を求め、配列として与える。
一様乱数を、求められた分布関数の逆関数に与えることで
IEDF に従った乱数が得られる。
http://hi.sakura.ne.jp/~nmaeda/web/fran.shtml
を参照
 
乱数は繰り返し作成する必要があるので、
1  IEDFのファイル読み込み -> distribution の作成

2  distribution に対応した乱数分布の作成

は別の関数にする

1  iedf_to_distribution(char[] ファイル名 , int 配列のサイズ,
                       -> 出力 double[] エネルギー  ,double[] 分布)
 
2  randomized_distribution   (long int 乱数の種, int 作りだす乱数の数,
                              double[] エネルギー  ,double[] 分布,
                     -> 出力  double[] (元の乱数), double[] 得られた乱数) 
*/
#ifndef _IEDF_TO_DISTRIBUTION_
 
#include<iostream>
#include<fstream>
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
// 下の定数は関数化する時に引数にする
/*
const long int ARRAY_SIZE = 2502 ;
long int RANDOM_SEED = 2 ;
const int N_CREATING_RANDOM  = 1000 ;
*/

void iedf_to_distribution(char file_name[], int array_size,
			  double energy_in_eV[] , 
			  double distribution[]  	)

{
  // ---- 配列 distribution のindex は 0 から array_size - 1 
  // ---- であることに注意 14 Feb 2003

  // -- 配列サイズがvariable である時には、動的に確保する
  
  double *amount_function;
  amount_function = new double[array_size + 1] ;


  input_array(file_name , energy_in_eV, amount_function, array_size) ;

  distribution[0] = amount_function[0] ;

  for(int i = 0; i < array_size - 1 ; i++) // <- 中で[i+1] のindexを使っている。
    {
      //cout << energy_in_eV[i] << "\t" << amount_function[i] << 
      //	"\t" << distribution[i] << endl;

      distribution[i+1] = distribution[i] + amount_function[i+1] ;
    }

  // 確率分布になるように正規化を行う
  // (distribution の最後の値が 1 となるようにする)

  double normalizer = 1.00 / distribution[array_size - 1] ;

  for(int i = 0; i < array_size - 1 ; i++)
    {
      distribution[i] *= normalizer ;
    
      //cout << energy_in_eV[i] << "\t" << amount_function[i] << 
      //"\t" << distribution[i] << endl;
    }
  distribution[array_size - 1] = 1.00 ; // 最後は1にしておく
  
  delete [] amount_function ;
      
} 

//-------------------------------
// 乱数発生(0-1) -> distribution を前から比較 -> 
// それに対応する energy_in_eV を出力

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
          // となったら、 i の要素を取りだす
          
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
  std::cout << input_energy_eV << "\t" << energy_in_eV[i] << "\t"
	    << distribution[i] << std::endl ;
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
