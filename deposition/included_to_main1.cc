
// included_to_main1.cc 
//=== monte_carlo.cc ��Ǥ� main �ˤ����륯�饹���ѿ������

// == �׻����ε�Ͽ


if(FLAG_SiCl4_INCLUDED == true)
{
  cout               << "sticking probability of SiCl4: " 
		     <<  SiCl4_DEPOSITION_PROBABILITY << endl ;
  condition_file     << "sticking probability of SiCl4: " 
		     <<  SiCl4_DEPOSITION_PROBABILITY << endl ;
}

if(FLAG_SiCl2_INCLUDED == true)
{
  cout               << "SiCl2 reincidence probability: " 
		     << SiCl2_REINCIDENCE_PROBABILITY << endl ;
  condition_file << "\nSiCl2 reincidence probability: " 
		 << SiCl2_REINCIDENCE_PROBABILITY << endl ;
}
ofstream file_SiCl2("file_SiCl2.dat");

if(FLAG_SiOCl2_INCLUDED == true )
{
  std::cout << "SiOCl2 reincidence probability: " 
	    << SiOCl2_REINCIDENCE_PROBABILITY 
	    << "\nsticking probability of SiOCl2: "
	    << SiOCl2_DEPOSITION_PROBABILITY << "\n" ;
  condition_file << "\nSiOCl2 reincidence probability: "
		 << SiOCl2_REINCIDENCE_PROBABILITY 
		 << "\nsticking probability of SiO: "
		 << SiOCl2_DEPOSITION_PROBABILITY << "\n" ;
}
ofstream file_SiOCl2("file_SiOCl2.dat");

//***************************************
//--- Etch product�Υ��󥹥��󥹤����

class Etch_product_class  SiCl4_neutral(CL_MASS, 0, FLAG_OPEN_SPACE, 4, 0,
					0.0, 
					SiCl4_DEPOSITION_PROBABILITY	) ;
class Etch_product_class  SiCl2_neutral(CL_MASS, 0, FLAG_OPEN_SPACE, 2, 0,
					SiCl2_REINCIDENCE_PROBABILITY, 
					SiCl2_DEPOSITION_PROBABILITY	) ;
class Etch_product_class  SiOCl2_neutral(CL_MASS, 0, FLAG_OPEN_SPACE, 2, 1,
					 SiOCl2_REINCIDENCE_PROBABILITY, 
					 SiOCl2_DEPOSITION_PROBABILITY	) ;
//int interval_step_SiO = 0 ;

//***************************************
//--- ���ǤΥ��󥹥��󥹤����
// -- Si�λ����ν����� Shape_trim ���饹�Υ��ФǹԤ��Τ�
// ���Τ�����ä˿��������饹�Ϻ��ʤ�

class Neutral_class   Oxygen_atom(CL_MASS, 0) ;

if(FLAG_OPEN_SPACE == true)//== ��¦ open space�ξ��
{
  SiCl4_neutral.set_flag_boundary(SPECULAR_REFLECT) ;
  SiCl2_neutral.set_flag_boundary(SPECULAR_REFLECT) ;
  SiOCl2_neutral.set_flag_boundary(SPECULAR_REFLECT) ;
  Oxygen_atom.set_flag_boundary(SPECULAR_REFLECT) ;
}

//======
// �롼����� f��˰�������¹�
// ��infiniti��C++ practice����򻲾ȡ�
int interval_step_oxygen = 0 ;



bool    flag_deposition  ;
flag_deposition = false  ;
//bool    flag_SiCl4_desorption = false ;

int i_exception_move ;
int i_exception_reflect ;

//cout << "included_to_main1.cc is included.\n" ;

if(FLAG_OXYGEN_INCLUDED == true)
{
  //***************************************
  cout           << "oxygen to ion flux ratio: " << OXYGEN_ION_RATIO 
		 << "\netch rate ratio SiO/Si: " << SELECTIVITY_SiO 
		 << "\netch rate ratio SiO2/Si: " << SELECTIVITY_SiO2 << endl ;
  condition_file << "\noxygen to ion flux ratio: " << OXYGEN_ION_RATIO 
		 << "\netch rate ratio SiO/Si: " << SELECTIVITY_SiO 
		 << "\netch rate ratio SiO2/Si: " << SELECTIVITY_SiO2 << endl ;
}



cout << "-------------------------" << endl ;
