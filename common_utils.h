//== common_utils.h ===


#ifndef _COMMON_UTILS_DEFINED_

//== 2�����٥��ȥ뤬Ϳ����줿���ˡ����β�ž�Ѥ��֤���
// x ���� 0 rad �Ȥ��ơ�0��2��
//   arctan ����ʬ�����ƻȤ�

//const double PI = atan(1.0) * 4.0 ;

double rotating_angle(double x, double y ) ;



//==  -> ��ĤΥ���ˤ�����Cl�ο��Ȼ��Ǥο���Ʊ��Υե������
//     ��Ͽ��������-> 10�ʿ��ǽ��ΰ̤���Ǥο�����ΰ̤�Cl�ο��Ȥ���
//     -> �ץ�åȤ���ݤ�mod��ȤäƤ��줾��ο�����Ф���
//    ���ϡ�Shape.n_oxygen , Shape.n_Clbond ,
//			     N_CELL_X, N_CELL_Z , 
//			     CELL_SIZE

void  output_cellinfo(int  **n_oxygen, int  **n_Cl,
		      int i_ion, char filename[] , int  n_x,  int  n_z,
		      double  cell_size ) ;

//== �ե��������� ==
void  input_cellinfo(int  **n_oxygen, int  **n_Cl,
		     char filename[] , int  n_x,  int  n_z ) ;


//==  ���ͥ��ͥ륮�����Ф��� etch yield ���֤� ==
//    yield �η����� = C(��E - ��Eth): C = 0.77, Eth = 20.0 
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
