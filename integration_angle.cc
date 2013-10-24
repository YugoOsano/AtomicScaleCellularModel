#include <iostream>
#include <math.h>

#include "integration_angle.h"
#include "atom_struct.h"

// 入力値：m 出力：eV
inline double V_sw(double r)
{
  // --- Stillinger-Weberポテンシャルにおける規格化
  // r: 2.0951Å, V: 2.17eV を１単位と見なす

  double tmp_v ;
  double r_sw  = fabs( r / 2.0951e-10 );

  // == r_cutoff = 1.8 (Feil et al.) == 
  if( r_sw < 1.8 )
    {
      tmp_v = 2.17 * 28.0 
	* (0.67 * pow(r_sw, -2.2) - pow(r_sw, -0.9)) 
	*  exp(  1.3 / (r_sw - 1.8) ) ;
    }
  else
    {
      tmp_v = 0.0 ;
    }

  return tmp_v ;
  
  //return 0.0 ; <- 以前、値が負の場合はゼロを返り値としていたが
  //                隣接格子点での衝突も考慮すると、負のポテンシャル
  //                領域も含める必要がある。
}

const double Dr     = 0.001e-10 ; // m
const double Lambda = 0.0001e-10 ;
//const double p  = 1.5e-10 ; //   m
//const double Ec = 200.0  ; // eV

//---- 入力：インパクトパラメータ p (m)
//           入射エネルギー(eV)
double theta_integrate(double p, double Ec)
{
  double r, tmp_f ;
  double tmp_integral = 0.0 ;
  
  for(int i = 1 ; i < 100000 ; i++)
    {
      r = Dr * i ;
      
      //√の中の式の値が非常に小さい場合、計算がおかしくなるので、
      // 1.0e-5 を超える場合のみカウントする
      
      if(1.0 - V_sw(r) / Ec - (p/r)*(p/r) > 1.0e-5 )
	{         
          tmp_f = p * Dr / (r * r * sqrt( 1.0 - V_sw(r) / Ec 
                                          - (p/r)*(p/r) ) ) ;
	  
	  tmp_integral += tmp_f ;

	  //std::cout << r     << "\t" << V_sw(r) << "\t" 
	  //    << 1.0 - V_sw(r) / Ec - (p/r)*(p/r) << "\t" 
	  //    << tmp_f << "\t" << tmp_integral << "\n" ;
	}
    }

  //  std::cout << PI - 2.0 * tmp_integral << "\n" ;


  return PI - 2.0 * tmp_integral ;
}

//------------------------
double theta_integrate_NC(double p, double Ec) 
{
  //  被積分関数が r^(-2)のorder なので
  //  積分間隔を λ((i+1)^2 - i^2) = λ(2i + 1) 
  //   とする。（λは定数）
  double r1, r2, tmp_f ;
  double tmp_integral = 0.0 ;
  
  for(int i = 1 ; i < 10000 ; i++)
    {
      r1 = Lambda * i * i ;
      r2 = Lambda * (i+1) * (i+1) ;
      if(1.0 - V_sw(r1) / Ec - (p / r1)*(p / r1) > 1.0e-5 )
	{
	  tmp_f = 0.5 * 
	    (p * Lambda * (2*i + 1) 
	     / (r1 * r1 * sqrt( 1.0 - V_sw(r1) / Ec 
				- (p/r1)*(p/r1))) +
	     p * Lambda * (2*i + 1) 
	     / (r2 * r2 * sqrt( 1.0 - V_sw(r2) / Ec 
				- (p/r2)*(p/r2))) ) ;
	  
	  tmp_integral += tmp_f ;
	}
    }
  return PI - 2.0 * tmp_integral ;

}

// test of function "theta_integrate"
/*
int main()
{
  double tmp_theta , tmp_psi;
  double Ec = 50.0 ;
  for(double tmp_p = 0.1e-10; tmp_p < 5.0e-10; tmp_p += 0.1e-10)
  //for(double tmp_p = 2.5e-10; tmp_p < 2.6e-10; tmp_p += 0.1e-10)
    {
      tmp_theta = theta_integrate_NC(tmp_p, Ec) ;

      tmp_psi   = atan( sin(tmp_theta) / 
			(cos(tmp_theta) 
			 + Cl.atomic_weight / Si.atomic_weight_solid )) ;

      std::cout << "p: "   << tmp_p * 1.0e+10 
		<< "\tE: " << Ec
		<< "\ttheta: " << tmp_theta 
		<< "\tpsi(deg): "   << tmp_psi * 180.0 / PI << "\n" ;

      //std::cout << tmp_p * 1.0e+10 
      //	<< "\t"  << tmp_psi * 180.0 / PI << "\n" ;
      
    }
  return 0 ;
}
*/
//== SW potential のプロット ==
/*
  int main()
  {
  for(double r = 0.1e-10; r < 5.0e-10; r += 0.1e-10)
  {
    std::cout << r << "\t" << V_sw(r) << "\n" ;
    }
    return 0 ;
    }
*/
