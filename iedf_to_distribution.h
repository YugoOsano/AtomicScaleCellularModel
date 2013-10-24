
#ifndef _IEDF_TO_DISTRIBUTION_H_DEFINED_


void iedf_to_distribution(char file_name[], int array_size,
			  double energy_in_eV[] , 
			  double distribution[]  	);

//-------------------------------
// ���ȯ��(0-1) -> distribution ����������� -> 
// ������б����� energy_in_eV �����
//
// append the (energy) threshold to promote the simulation efficiency

void randomized_distribution(long int random_seed, int n_creating_random,
			     double energy_in_eV[], double energy_threshold, 
                             double distribution[] ,
			     double random_sample[] , 
			     double randomized_distribution_array[] ) ; 

// =======================
// iedf_to_distribution ������줿 energy_in_eV / distribution
// �Υ��åȤˤ����ơ�
// ���ͥ륮�������� -> �б�����distribution ���֤��ؿ�

double distribution_from_energy(double input_energy_eV,
                                double energy_in_eV[] , 
                                double distribution[]  ,
				int array_size );

#define _IEDF_TO_DISTRIBUTION_H_DEFINED_
#endif
