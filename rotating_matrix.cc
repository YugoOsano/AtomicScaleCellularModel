
#include <math.h>
#include <iostream.h>
#include "atom_struct.h"
#include "rotating_matrix.h"

void rotating_matrix(double theta , double psi,
                     bool  sign_flag , 
                     double *a11, double *a12, double *a13,
                     double *a21, double *a22, double *a23,
                     double *a31, double *a32, double *a33 )
{
  double f , a21_min ;
  double f_minimum = 5.0 ;
  double tmp_a11, tmp_a12, tmp_a21, tmp_a22 ;
  a21_min = 0.0 ;

  for(tmp_a21 = - 1.0; tmp_a21 <= 1.0; tmp_a21 += 0.005)
    {
      tmp_a11 = cos(psi) - tmp_a21 * tan(theta) ;

      tmp_a12 = (tmp_a11 - 1.0) * tan(theta) ;
      tmp_a22 = 1.0 + tmp_a21 * tan(theta) ;
      
      f = (tmp_a11 * tmp_a21 + tmp_a12 * tmp_a22)
        * (tmp_a11 * tmp_a21 + tmp_a12 * tmp_a22)
        - (1.0 - tmp_a11 * tmp_a11 - tmp_a12 * tmp_a12)
          * (1.0 - tmp_a21 * tmp_a21 - tmp_a22 * tmp_a22) ;

      if(f < f_minimum)
        {
          f_minimum = f ;
          a21_min = tmp_a21 ;
        }
    }
  //cout << a21_min << "\t" << f_minimum << endl ;

  //---- s—ñŒˆ’è
  *a21 = a21_min ;
  
  *a11 = cos(psi) - (*a21) * tan(theta) ;

  *a12 = ((*a11) - 1.0) * tan(theta) ;
  *a22 = 1.0 + (*a21) * tan(theta) ;

  //-- ˆÈ‰ºsqrt‚ÉŠÖ‚µ‚Ä³•‰‚Ì”CˆÓ«—L‚è
  if( 1.0 - (*a11) * (*a11) -  (*a12) * (*a12) <= 0.0 )
    *a13 = 0.0 ;
  else if(sign_flag == true)
    *a13 = sqrt(1.0 - (*a11) * (*a11) -  (*a12) * (*a12) ) ;
  else
    *a13 = - sqrt(1.0 - (*a11) * (*a11) -  (*a12) * (*a12) ) ;

  
  if( fabs(*a13) > 0.0001 )
    {
      *a23 = - ((*a11) * (*a21) + (*a12) * (*a22)) / (*a13) ;

      if(sign_flag == true)
        *a31 = - 1.0 
          / sqrt(1.0 + tan(theta) * tan(theta)
                 + pow( ((*a11) + (*a12) * tan(theta))/ (*a13), 2.0 ) ) ;
      else
        *a31 = 1.0 
          / sqrt(1.0 + tan(theta) * tan(theta)
                 + pow( ((*a11) + (*a12) * tan(theta))/ (*a13), 2.0 ) ) ;
      
      *a32 = (*a31) * tan(theta) ;
      *a33 = - ((*a11) * (*a31) + (*a12) * (*a32)) / (*a13) ;
      //sqrt(1.0 - a31 * a31 -  a32 * a32) ;
    }
  else
    {
      *a23 = 0.0 ;
      *a31 = 0.0 ;
 
      *a32 = 0.0 ;
      *a33 = 1.0 ;
    }
}
//--- test main
/*
int main()
{
  double a11, a12, a13 ;
  double a21, a22, a23 ;
  double a31, a32, a33 ;
  
  rotating_matrix( - 0.75 * PI , 0.3 * PI ,
                   &a11, &a12, &a13 ,
                   &a21, &a22, &a23 ,
                   &a31, &a32, &a33 ) ;
        
  cout << a11 << "\t" << a12 << "\t" << a13 << endl ;
  cout << a21 << "\t" << a22 << "\t" << a23 << endl ;
  cout << a31 << "\t" << a32 << "\t" << a33 << "\n\n" ;

  // -- ŒŸŽZ
  cout << a11 * a11 + a12 * a12 + a13 * a13 << "\t"
  << a11 * a21 + a12 * a22 + a13 * a23 << "\t"
  << a11 * a31 + a12 * a32 + a13 * a33 << "\n"
  
  << a11 * a21 + a12 * a22 + a13 * a23 << "\t"
  << a21 * a21 + a22 * a22 + a23 * a23 << "\t"
  << a21 * a31 + a22 * a32 + a23 * a33 << "\n"

  << a11 * a31 + a12 * a32 + a13 * a33 << "\t"
  << a21 * a31 + a22 * a32 + a23 * a33 << "\t"
  << a31 * a31 + a32 * a32 + a33 * a33 << "\n" ;
  return 0 ;


}
  */
