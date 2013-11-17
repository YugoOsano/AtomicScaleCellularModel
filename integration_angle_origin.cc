#include <iostream.h>
#include <math.h>
#include "atom_struct.h"

// 入力値：m 出力：eV
inline double V_sw(double r)
{
  // --- SWポテンシャルにおける規格化
  // r: 2.0951Å, V: 2.17eV を１単位と見なす
  double r_sw  = r / 2.0951e-10 ;

  double tmp_v = 2.17 * 28.0 * (0.67 * pow(r_sw, -2.2) 
                                - pow(r_sw, -0.9)) 
    *    exp(  1.3 / (r_sw - 1.8) ) ;
  
  if(tmp_v >= 0.0)
    return tmp_v ;
  
  return 0.0 ;
}

const double Dr = 0.005e-10 ; // m
const double p  = 1.5e-10 ; //   m
const double Ec = 200.0  ; // eV

int main ()
{
  double r, tmp_f ;
  double tmp_integral = 0.0 ;
  
  for(int i = 1 ; i < 1000 ; i++)
    {
      r = Dr * i ;
      
      if(1.0 - V_sw(r) / Ec - (p/r)*(p/r) >= 0.0 )
	{
	  tmp_f = p * Dr / (r * r * sqrt( 1.0 - V_sw(r) / Ec 
					  - (p/r)*(p/r) ) ) ;
	  
	  tmp_integral += tmp_f ;

	  cout << r     << "\t" << V_sw(r) << "\t" 
	       << tmp_f << "\t" << tmp_integral << endl ;
	}
    }

  cout << PI - 2.0 * tmp_integral << endl ;


  return 0 ;
}
