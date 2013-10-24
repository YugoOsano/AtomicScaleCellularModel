//== common_utils.h ===


#ifndef _COMMON_UTILS_DEFINED_

//== 2次元ベクトルが与えられた時に、その回転角を返す。
// x 軸を 0 rad として、0〜2π
//   arctan を場合分けして使う

//const double PI = atan(1.0) * 4.0 ;

double rotating_angle(double x, double y ) ;



//==  -> 一つのセルにおけるClの数と酸素の数を同一のファイルで
//     記録したい。-> 10進数で十の位を酸素の数、一の位をClの数とする
//     -> プロットする際にmodをとってそれぞれの数を抽出する
//    入力：Shape.n_oxygen , Shape.n_Clbond ,
//			     N_CELL_X, N_CELL_Z , 
//			     CELL_SIZE

void  output_cellinfo(int  **n_oxygen, int  **n_Cl,
		      int i_ion, char filename[] , int  n_x,  int  n_z,
		      double  cell_size ) ;

//== ファイル入力 ==
void  input_cellinfo(int  **n_oxygen, int  **n_Cl,
		     char filename[] , int  n_x,  int  n_z ) ;


//==  入射エネルギーに対して etch yield を返す ==
//    yield の係数β = C(√E - √Eth): C = 0.77, Eth = 20.0 
//    *: 0.77 * sqrt(20.0) = 3.44354468534968

inline double yield_beta(double incident_energy)
{
  double  tmp ;
  tmp = 0.77 * sqrt(incident_energy) - 3.44354468534968 ;

  if(tmp > 0.0)
    return tmp ;
  else 
    return 0.0 ;
}


//*********************


#define  _COMMON_UTILS_DEFINED_
#endif
