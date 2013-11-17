
#include <iostream.h>
#include <math.h>
#include "atom_struct.h"

// max, min 関数
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))




//----------------------------------------------------
// proton electronic stopping powers (RPSTOP関数/p-230)

// 入力：固体原子の原子番号、入射エネルギー、struct
// 出力：stopping power

double range_proton_stop(int z2, double energy_i,
			 Atom_struct Atom )
{
  if(energy_i < 25.0)
    energy_i = 25.0 ;
  
  double sl, sh ;
  sl = Atom.atomic_n * pow(energy_i, Atom.atomic_mass) +
    Atom.atomic_weight * pow(energy_i, Atom.atomic_weight_solid) ;
  
  sh = Atom.density_mass / pow(energy_i, Atom.density_atom) 
    *  log( Atom.fermi_v / energy_i + Atom.factor_ion_screening * energy_i) ;
  
  if(energy_i > 25.0)
    return sl * sh / (sl + sh) ;

  if(z2 <= 6)
    return sl * sh / (sl + sh) * pow(energy_i / 25.0 ,0.25)  ;

  return sl * sh / (sl + sh) * pow(energy_i / 25.0 ,0.45)  ;
}

//==========================================
// -- Ion stopping cross-section function used in TRIM85
//  (RSTOP/p-228)
//
// 入力： ion (struct), target (struct), 入射エネルギー、
//      lambda screening factor for ions, Fermi velocity of solid/v0 (?)

// 出力： electronic stopping (get 1000 values)

void range_stop(Atom_struct  Ion, Atom_struct Target_solid,
		double energy_i,  double lambda_screening,
		double v_fermi ,  
		double elec_stop[] )
{
  double V_EFFECTIVE_MIN = 0.13 ;
  double vr_min ;
  double tmp_energy ;
  double tmp_v , tmp_vr, v_effective ;

  double tmp_a , tmp_b  ; // -- ionization level を表す時の指数に付く係数
  double tmp_l0, tmp_l1, tmp_l ;
  double tmp_zeta , tmp_power ;

  double q_ionization_level ;


  for(int i = 1; i <= 1000; i++)
    {
      vr_min = 1.0  ;
      // 0 〜 energy_i までを1000段階に分ける
      tmp_energy = i * 0.001 * energy_i / Ion.atomic_weight ;

      //--- heavy ion のための換算
      tmp_v = sqrt( tmp_energy / 25.0) / v_fermi ;
      
      if (tmp_v < 1.0)
	tmp_vr = (3.0 * v_fermi / 4.0) 
	  * (1.0 + (2.0 * tmp_v * tmp_v / 3.0 ) - pow(tmp_v , 4.0)/15.0 ) ;
      else
	tmp_vr = tmp_v * v_fermi * ( 1.0 + 1.0 / (5.0 * tmp_v * tmp_v)) ;

      v_effective = max( V_EFFECTIVE_MIN, tmp_vr / pow(Ion.atomic_n, 0.6667) ) ;
      v_effective = max( v_effective , vr_min / pow(Ion.atomic_n, 0.6667) ) ;
		    
      tmp_a = - 0.803 * pow(v_effective, 0.3) + 1.3167 * pow(v_effective, 0.6)
	+    0.38157 * v_effective  + 0.008983 * v_effective * v_effective  ;

      q_ionization_level = min( 1.0, max (0.0, - exp( - min (tmp_a, 50.0)))) ;

      // -- 'effective charge' に変換
      tmp_b = (min (0.43, max(0.32, 0.12 + 0.025 * Ion.atomic_n ))) 
	/   pow( Ion.atomic_n, 0.3333 ) ;

      tmp_l0 = (0.8 - q_ionization_level * (min(1.2, 0.6 + Ion.atomic_n / 30.0)))
	/  pow( Ion.atomic_n, 0.3333 ) ;

      if(q_ionization_level < 0.2)
	tmp_l1 = 0.0 ;

      else if (q_ionization_level < ( max(0.0, 0.9 - 0.025 * Ion.atomic_n)) )
	tmp_l1 = tmp_b * (1.0 - q_ionization_level) 
	  /   fabs(max(0.0, 0.9 - 0.025 * Ion.atomic_n) - 0.2000001 ) ;

      else if (q_ionization_level < 
	       (max(0.0, 1.0 - 0.025 * min(16.0, 1.0 * Ion.atomic_n ))) )
	tmp_l1 = tmp_b ;

      else 
	tmp_l1 = tmp_b * (1.0 - q_ionization_level)
	  /    (0.025 * min(16.0, 1.0 * Ion.atomic_n )) ;

      tmp_l = max(tmp_l1, tmp_l0 * lambda_screening ) ;
      tmp_zeta = q_ionization_level 
	+ (1.0 / (2.0 * pow(v_fermi, 2.0))) * (1.0 - q_ionization_level) 
	* log(1.0 + pow( (4.0 * tmp_l * v_fermi / 1.919), 2.0) )  ;
      
      tmp_a = - pow( (7.6 - max(0.0, log(tmp_energy) )) , 2.0) ;
      tmp_zeta = tmp_zeta * (1.0 + (1.0/ pow(Ion.atomic_n, 2.0)) * 
			     (0.18 + 0.0015 * Ion.atomic_n ) * exp(tmp_a) ) ;

      if(v_effective > max(V_EFFECTIVE_MIN, vr_min / pow(Ion.atomic_n, 0.6667)) )
	{
	  elec_stop[i] = range_proton_stop(Target_solid.atomic_n, tmp_energy,
					   Target_solid )
	    * pow( (tmp_zeta * Ion.atomic_n), 2.0);
	}
      else
	{
	  vr_min = max(vr_min, v_effective * pow(Ion.atomic_n, 0.6667)) ;
	  
	  double v_min = 0.5 * (vr_min + sqrt( max(0.0, vr_min * vr_min 
						   - 0.8 * v_fermi * v_fermi )) ) ;

	  if( (Target_solid.atomic_n == 6) ||
	      ( ( Target_solid.atomic_n == 14 ) ||
		( Target_solid.atomic_n == 32 ) ) &&
	      ( Ion.atomic_n == 19 ) )
	    tmp_power = 0.375 ;
	  else
	    tmp_power = 0.5 ;
     

	  elec_stop[i] = (range_proton_stop(Target_solid.atomic_n, 25.0 * v_min * v_min,
					    Target_solid )
			  * pow( (tmp_zeta * Ion.atomic_n) , 2.0)) 
	    * pow( tmp_energy / (25.0 * v_min * v_min) , tmp_power ) ;
	}
    }
}

//------------------------------------------
int main()
{
  
  double *tmp; 
  tmp = new double[1100] ;

  range_stop(Cl, Si,
	     10.0, Si.factor_ion_screening ,
	     Si.fermi_v,
	     tmp) ;
	    
  for(int i = 1 ; i <= 1000 ; i++)
    cout << tmp[i] << endl ;


  
  delete [] tmp ;
  return 0 ;
}




