
#ifndef _IEDF_TO_DISTRIBUTION_H_DEFINED_


void iedf_to_distribution(char file_name[], int array_size,
			  double energy_in_eV[] , 
			  double distribution[]  	);

//-------------------------------
// 乱数発生(0-1) -> distribution を前から比較 -> 
// それに対応する energy_in_eV を出力
//
// append the (energy) threshold to promote the simulation efficiency

void randomized_distribution(long int random_seed, int n_creating_random,
			     double energy_in_eV[], double energy_threshold, 
                             double distribution[] ,
			     double random_sample[] , 
			     double randomized_distribution_array[] ) ; 

// =======================
// iedf_to_distribution で得られた energy_in_eV / distribution
// のセットにおいて、
// エネルギーを入力 -> 対応するdistribution を返す関数

double distribution_from_energy(double input_energy_eV,
                                double energy_in_eV[] , 
                                double distribution[]  ,
				int array_size );

#define _IEDF_TO_DISTRIBUTION_H_DEFINED_
#endif
