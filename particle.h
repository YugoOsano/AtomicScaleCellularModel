
// particle.h

// γ�Ҥ⥯�饹�Ȥ��ư�����

#ifndef _PARTICLE_H_DEFINED_

#include "atom_struct.h"

//== �������Τ����flag == 
const int NO_PERIODIC      = 1;//== ���������ʤ� ==
const int SPECULAR_REFLECT = 2;//== �����Ƕ���ȿ�� ==

//
class Particle_class
{
public:
  struct Particle_location_velocity_struct pos_v ;
  //--- �ʤޤ������β��ΰ���
  struct Particle_location_velocity_struct pos_v_togo ;

  //===== 2004/08/04 �ɲ�
  //  ����γ�Ҥ����֤��륻���ֹ桢�����������Ū�ʰ���
  // Shape_trim_class::put_shape_matrix �ǹԤäƤ����׻�
  // �ΰ���ʬ�Ǥ��롣
  /*     dx
       +--------+
       |   |    |
     dz|   |    |
       |---*    |
       |        |
       +--------+
  */
  int    i_cell_x, i_cell_z ;
  double dx_cell , dz_cell  ;


  //--- ����
  double mass ;
protected:
  //---��ư���ͥ륮��
  double energy_k ;
public:


  //=========================================
  //===  ���ͳ���ʬ�� =======================
  int n_adf_array ;  // ����ʬ�۴ؿ��γ������ 1 ����ꤹ����������ͤȤ���
  double *angular_df ;
  double *angle_array_degree ;
  
  //--- ɽ�̤˿������Ƥ��뤫�ɤ����� flag
  // ������Ȥ�æΥ��Ԥ����ɤ���Ƚ�Ǥ��롣
  // -- ����ͤ� FALSE, inject ��FALSE ���᤹��
  bool   flag_contact ; 

protected: 
  //=============================
  //===  ����γ�ҿ��Υ����� ===
  int n_injection ;

  //== �㳰�����Τ���Υ�����
  int i_exception_move ;

  //== �ΰ������ˤ��뤫Ƚ�ꤹ�뤿��Υե饰��������true��
  bool flag_inside ;
  bool flag_inside_togo ;

  //== x�����ζ�������ɽ���ե饰 === 
  // ��0������������1�����������ʤ���2�������Ƕ���ȿ�͡�
  int  flag_boundary ;

public:
  //-- constructor
  Particle_class() ;
  Particle_class(double mass_input, int flag_boundary_input ) ;
  ~Particle_class() ;

  int  put_n_injection(){ return n_injection ;}
  void add_n_injection(){ n_injection++ ;}
  
  bool put_flag_inside(){ return flag_inside ;}

  //== flag_boundary ���ѹ� ==
  void set_flag_boundary(int input_flag) 
  { flag_boundary = input_flag ; }

  //=======================================
  //== ����ʬ�ۤ�ե����뤫��������� =====
  // from "flux_model" directory
  // input : file name
  void read_angular_df(char filename[]);

  //-- γ�����͡���ɸ��x �����濴��®�٤�z�����Τ�Ϳ����
  void inject_from_center( double v_z ) ;

  //-- γ�����͡�x��ɸ�������®�٤�z�����Τ�Ϳ����
  void inject_from_top(double v_z) ;

  //-- γ�����͡��Ф������˽��®�٤�Ϳ����
  // �Ȥ��ϰϤ� -��/2 �� ��/2 
  //  ���ϡ����֡����ͳѡ�®��
  void inject_oblique(double x, double psi, double v) ;


  //===== 2004/08/04 �ɲ�
  //  ����γ�Ҥ����֤��륻���ֹ桢�����������Ū�ʰ���
  // ( i_cell_x, i_cell_z , dx_cell, dz_cell ) ��׻�����
  void get_position() ;

  //=================================
  //-- Ĺ�� l �����¿ʱ�ư������
  //     �ΰ褫��Ф���� flag_inside ��false�ˤ���
  //
  // -- ����������������
  //================================
  void  move_trans(double  l) ;
    
  // -- pos_v_togo �����
  void  move_trans_togo(double  l ) ;

  //=== 06 Aug 2004 modified 
  // ����ɽ�̤ˤ����� thermal reflection
  //-- 
  // flag FALSE -> ñ���®�٤��ѹ� 
  //      TRUE  -> cosine distribution ��Ƴ��
  // ��ˡ�����������ϡ�
  // ˡ�������򣰡�Ȥ���distribution �����ȡ�
  // 0��1 �����R ���Ф���
  // asin(2R-1)

  void random_reflection(bool   flag_cosine_dist ,
			 double normal_x, double normal_z ) ;

  //-------------------------
  //-- nuclear energy loss (T) ����Ϥ���
  // ���ϡ����ͥ����󡢥������åȸ��Ρ����ͥ��ͥ륮��(eV)��
  //                   impact parameter (p)
  // ���ϡ�scattering angle

  double put_energy_loss(Atom_struct  Ion,   Atom_struct  Solid,
			 double  energy_i ,
			 double  b_impact_para,
			 double  *theta_c ,    double  *psi_scattering ) ;

  //--------------------------
  // ���Υ����󡿸��Ҥξ���
  // �֤��͡�ķ���֤ä����ҤΥ��ͥ륮��: (tmp_energy)
  void  collision_with_solid_atom(double p_impact_parameter ) ;

  //==========================
  // �����󡿸��Ҥξ���
  // ����������������ɽ�����ᡢimpact parameter��
  // γ�Ҥΰ��֡�®�٤���׻����Ƶ��롣
  // ��ά���Τ��ᡢ�������η׻��ΤߤȤ��롣
  // �ʱ��Ԥ�������®�٤ϰ���Ȥ����

  // �����͡����θ��Ҥΰ���
  void  collision_accurate(double x_solid, double z_solid,
			   double *tmp_energy ) ;

  // �����󡿸��Ҥξ��͡�3�����׻�

  /* ������ΰ��֤�(x, y, z), ®�٤�(vx, vy, vz), 
     Si���Ҥΰ��֤�(X, Y, Z) �Ȥ��ơ�������-���Ҥκ��ܶ�ΰ��֤�
     ���롣�ѥ�᡼���ˡ��٥��ȥ�r = (rx, ry, rz) �����Ѥ��ơ�
     
     (x, y, z)+��(vx, vy, vz) = (X, Y, Z)+(rx, ry, rz) 
     ��򤯡�
     �� = [(X-x)vx + (Y-y)vy + (Z-z)vz]/(vx^2 + vy^2 + vz^2)

     rx = x - X + �� vx
     ry = y - Y + �� vy
     rz = z - Z + �� vz 

     impact parameter �� |r|

     ������®�٥٥��ȥ�� v = (vx, vy, vz) �� r
     �κ��ʿ�̾�ˤ���Τǡ������ΰ켡��硣
     ���줾��ñ�̥٥��ȥ�[v],[r]���ᡢ����Ѥ��פǤ���Ȥ���ȡ�
     v' = |v'|([v] cos�� + [r] sin��)
     �Ȥʤ롣

     ���֥٥��ȥ��Si���Ҥ�ꡢv'�ȿ�ľ�������� |r|����Υ�줿
     ���� (X, Y, Z) + |r|([v](-sin��) + [r]cos��) �Ȥ��롣
  */
  void  collision_accurate3D(double x_solid, double y_solid,
			     double z_solid, double *tmp_energy ) ;

} ;



#define _PARTICLE_H_DEFINED_
#endif

