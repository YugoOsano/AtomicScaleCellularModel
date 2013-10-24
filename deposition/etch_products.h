// etch_product.h 
// 
// 2004/07/29 -
//  ���å�������ʪ�Υ��饹��Neutral_particle 
//  ���饹��Ѿ����롣deposition ��᥽�åɤȤ��Ƽ�������

#ifndef _ETCH_PRODUCT_H_DEFINED_

#include "neutral_particle.h"

//== �����Υ���ξ��֤�Ƚ�ꤹ�뤿�ᡢshape���饹��include���롣
#include "shape_trim.h"

class Etch_product_class : public Neutral_class
{
protected:
  //=== SiClx �ˤ�������ǲ��ο� x ====
  int    n_chlorination ;

  //=== SiOy �ˤ�����oxidation�ο� y ==== 
  int    n_oxidation ;

  //=== ���巸�� Sp ===================
  double sticking_probability ;
  double sticking_probability_1st_contact ;

  //=== ���ߤΥ�������ͤ���������Υ����ֹ�
  int  i_previous_x, i_previous_z ;

public:
  //===�¸�Ū�����ͳ�Ψ
  //     �ѥ��������ˤ�ä����ͳ�Ψ���Ѳ����롧
  // �������Τ�Si��Ϫ�Ф�x�ѡ�����ȤȤ��ơ�simulation domain��
  // Si��Ϫ�Ф�y�ѡ�����ȤȤ���ȡ������ΰ�Ǥ����Τ�y/x�ܤ����Ф����롣
  //  �������äơ�Pr �� (x/y)���¸�Ū�����ͳ�Ψ�Ȥʤ� �� x = 0.5 �Ȥ����
  
  double reincidence_probability ;

public:
  //-- constructor 
  Etch_product_class() ;
  Etch_product_class(double mass_input, int flag_boundary_input,
		     bool   flag_open_space, 
		     int    n_chlorination_input ,
		     int    n_oxidation_input ,
		     double reincident_probability_whole,
		     double sticking_probability_input   );
  ~Etch_product_class();


  //-- deposition �Υ᥽�å�
  // ñ��� ballistic deposition model ���������
  // --> ���ꡢ����(2002, 2003) �����ز��򻲾�

  // 1. ������ʬ(SHAPE_Si)�����ͤ�����Ƚ��

  // =============================================
  // ���ߤΰ��֡�®�٥٥��ȥ��ꡢ�ɤΥ��뤫��
  // ���ߤΥ�������ͤ������׻����롣

  //(1) v_x ������ˤ�äƾ��ʬ����Ԥ�
  //(2) ®�٥٥��ȥ�ε��������Ĺ�������Ҥΰ��֤���Ƚ���Ԥ�

  // v_x > 0 �ʤ�С�0���Ҥ���
  //  -> ���Ҥε������(��x, ��z)���̤ꡢ�����٥��ȥ뤬
  //   ( - v_x , - v_z) ��ľ����0���Ҥϡ�
  //    ��z - (v_z/v_x)��x 

  // v_x < 0 �ʤ�С�x = CELL_SIZE �����Ҥ���
  // ���Ʊ���������Ҥ�
  //    ��z + (v_z/v_x)( CELL_SIZE - ��x)

  void get_previous_cell() ;


  //=========================================
  // ��SiCl4��æΥ������Ȥ��ơ�Siɽ�̤�γ�Ҥ����֤�
  //��®��Ϳ���롧���� whole_flight �ؿ��γȻ�ȿ�ͤ���ʬ�����

  // ���ϡ�Shape_trim_class�Υݥ��󥿡�
  //       ���֡�®�٤ι�¤�Ρʰ��־���Τߡˡ�®��
  // �֤��͡����֤��Ԥ�줿��� TRUE �����Ǥʤ���� FALSE

  bool allocate_on_surface(Shape_trim_class *Shape_pointer,
			   Particle_location_velocity_struct  pos_v_input,
			   double v ) ;

  // ===================================
  // ���饹��Shape_trim_class�ˤΥݥ��󥿤����Ѥ��ơ�
  // γ�Ҥ����Ͱʹߤν�����ؿ�������

  // �ʵ��夻���ˡ��ΰ�γ��˽Ф���硢true ���֤�
  //  :  æΥ���� SiCl4 �ο��˱��������ͤ���SiCl2�ο�����ꤹ�뤿��
  void whole_flight(Shape_trim_class *Shape_pointer,
		    bool flag_inject_from_side ) ;

};




#define _ETCH_PRODUCT_H_DEFINED_
#endif

