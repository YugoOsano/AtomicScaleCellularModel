
// included_to_main1.cc 
//=== monte_carlo.cc 内での main におけるクラス、変数の宣言

// == 計算条件の記録


if(FLAG_SiCl4_INCLUDED == true)
{
  std::cout               << "sticking probability of SiCl4: " 
		     <<  SiCl4_DEPOSITION_PROBABILITY << std::endl ;
  condition_file     << "sticking probability of SiCl4: " 
		     <<  SiCl4_DEPOSITION_PROBABILITY << std::endl ;
}

if(FLAG_SiCl2_INCLUDED == true)
{
  std::cout               << "SiCl2 reincidence probability: " 
		     << SiCl2_REINCIDENCE_PROBABILITY << std::endl ;
  condition_file << "\nSiCl2 reincidence probability: " 
		 << SiCl2_REINCIDENCE_PROBABILITY << std::endl ;
}
std::ofstream file_SiCl2("file_SiCl2.dat");

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
std::ofstream file_SiOCl2("file_SiOCl2.dat");

//***************************************
//--- Etch productのインスタンスを作成

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
//--- 酸素のインスタンスを作成
// -- Siの酸化の処理は Shape_trim クラスのメンバで行うので
// そのために特に新しいクラスは作らない

class Neutral_class   Oxygen_atom(CL_MASS, 0) ;

if(FLAG_OPEN_SPACE == true)//== 片側 open spaceの場合
{
  SiCl4_neutral.set_flag_boundary(SPECULAR_REFLECT) ;
  SiCl2_neutral.set_flag_boundary(SPECULAR_REFLECT) ;
  SiOCl2_neutral.set_flag_boundary(SPECULAR_REFLECT) ;
  Oxygen_atom.set_flag_boundary(SPECULAR_REFLECT) ;
}

//======
// ループ中で f回に一回処理を実行
// （infinitiのC++ practiceの例を参照）
int interval_step_oxygen = 0 ;



bool    flag_deposition  ;
flag_deposition = false  ;
//bool    flag_SiCl4_desorption = false ;

int i_exception_move ;
int i_exception_reflect ;

//std::cout << "included_to_main1.cc is included.\n" ;

if(FLAG_OXYGEN_INCLUDED == true)
{
  //***************************************
  std::cout      << "oxygen to ion flux ratio: " << OXYGEN_ION_RATIO 
		 << "\netch rate ratio SiO/Si: " << SELECTIVITY_SiO 
		 << "\netch rate ratio SiO2/Si: " << SELECTIVITY_SiO2 
		 << std::endl ;
  condition_file << "\noxygen to ion flux ratio: " << OXYGEN_ION_RATIO 
		 << "\netch rate ratio SiO/Si: " << SELECTIVITY_SiO 
		 << "\netch rate ratio SiO2/Si: " << SELECTIVITY_SiO2 
		 << std::endl ;
}



std::cout << "-------------------------" << std::endl ;
