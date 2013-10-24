// neutral_particle.h

// ����γ�ҤΥ��饹��Particle class ��Ѿ����롣

#ifndef  _NEUTRAL_PARTICLE_H_DEFINED_

#include "particle.h"

class Neutral_class : public Particle_class
{
protected:
  bool flag_adsorption ; //== ����򼨤������flag ;
  int  i_exception_reflect; //== ��ȿ�ͤˤ������㳰�����Τ���Υ�����

public:
  //-- constructor 
  Neutral_class() ;
  Neutral_class(double mass_input,  int flag_boundary_input);
  ~Neutral_class();

  //-- x��ɸ�������®�٤�x,z������Ϳ����
  // ��y ������®�٤�Ϳ����ɬ�פϤʤ���

  // -- ������cosine distribution�Ƿ���
  //  -> ���Υ����ɤǤ� x��������0��ʤΤǡ�sin��
  // distribution �����ȡ�[-cos��]�ȡ�0��/ [-cos��]180��� 0��
  // 0��1�����R ���Ф���
  // acos(1-2R)

  // ���ϡ��¿�®��
  void inject_from_top( double v ) ;

  //-- ��¦���������
  void inject_from_right_side(double v) ;
  
  //==================================
  // 06 Aug 2004
  // ***  substrate surface �ˤ�������ȿ��
  // monte_carlo.cc ��main ��ǹԤäƤ����������Υ��饹���
  // �Ԥ��褦�� refactoring 
 

  //=================================
  // ����γ�Ҥ����ͤ������С�����ޤ�
  // ���ϡ������������ΰ����(false)����¦(true)
  void all_process(class Shape_trim_class    *Shape,
		   class Shape_counter_class *Neutral_counter,
		   bool   flag_inject_from_side ,
		   bool   flag_flux_count ,
		   double incident_energy ,
		   bool   flag_chemical_etch,
		   double yield_chemical) ;
  

};


#define  _NEUTRAL_PARTICLE_H_DEFINED_
#endif
