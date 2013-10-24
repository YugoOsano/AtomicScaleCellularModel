// shape_trim.h
//    -> ������Ÿ�ξ��֤�ɽ������Υ��饹��

// ���ĤΥ���ǣ��Ĥ�Si ���Ҥ�ɽ����ΤȤ��롣�ʤ�������������
// shate_state ���饹���顢trim ��ɬ�פ����ǤΤߤ�ȴ���Ф����Ѥ��롣



#ifndef _SHAPE_TRIM_DEFINED_

const int SHAPE_SPACE = 0 ;
const int SHAPE_Si    = 10 ;
const int HARD_MASK   = 5 ;

const int CL_DEPTH_RANGE = 20 ;
#include "atom_struct.h"

class Shape_trim_class
{
public:
  // �ᥤ���matrix
  int *shape_matrix[N_CELL_X]; int *p_shape_matrix; 
  
  // --- ���줾���Si ���Ҥ˷�礷�Ƥ��� Cl �ο�
  int *n_Clbond[N_CELL_X];     int *p_n_Clbond  ;

  // --- Ʊ���� Si���Ҥ˷�礷�Ƥ��� O�饸����ο�
  int *n_oxygen[N_CELL_X];     int *p_n_oxygen  ;
  
  // *** ɽ�̤ˤ�����ˡ��������

  // * ���٤ƤΥ���ˤ����Ʒ׻����롣��ɽ�̤�̵�����ϣ��٥��ȥ��
  // * ���� -> ��������������
  // * Four-point calculation �ˤ�äƷ׻�

  double *surfacenormal_x[N_CELL_X] ; double *p_surfacenormal_x ;
  double *surfacenormal_z[N_CELL_X] ; double *p_surfacenormal_z ;

  //---- æΥ����Si���ҤΥ����󥿡��ϡ��ɥޥ����Υ�����
  //     (desorb_Si/desorb_mask�ؿ�)
  int cntr_desorbed_Si ; 
  int cntr_desorbed_mask ;
  
  //---- ���Ѥ���Si���ҤΥ�����
  //     (deposit_Si�ؿ�)
  int cntr_deposited_Si ;
  
  // -- domain �β��ޤ���ã�����顢 TRUE�Ȥ��롣
  bool flag_get_bottom ;

  //========================================
  // ���夷�����֤�Ͽ���Ƥ�������Υ��С�
  // md2d ��Ʊ�ͤˡ����饹��Υ��饹�ˤ���
  class Sticking_Cl_class
  {
  public:
    double x ; // ����
    double z ; 
    
  private:
    Sticking_Cl_class *next_ptr   ;
    friend class Shape_trim_class ;
  };
public:
  //______________________________
  Sticking_Cl_class *first_ptr ;
  int   n_sticking ;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //**********************************************
public:
  
  // ���󥹥ȥ饯�� 
  Shape_trim_class() ;
  Shape_trim_class(bool flag_open_space) ;
  ~Shape_trim_class() ;

  // -- ��ɸ�����Ϥ����� shape_matrix ���֤��ؿ�
  // -- ����͡��⤷�����ϰϳ����ͤ����Ϥ������ϣ����֤�

  //===============================================
  //===  ����Υ���ǥå�����ݥ��󥿤��֤�
  //===============================================
  int put_shape_matrix(double x, double z, 
		       int *i_x, int *i_z ) ;

  //===============================================
  //    ����ʻ����ˤ����ơ����������⤷������������
  //     ������ʬ(SPACE�Ǥʤ���) ������С�����index��ݥ��󥿤��֤�

  //    �ɲá�γ�Ҥΰ���(x, z)�Υ����index���֤�
  //     put_shape_matrix�ؿ���Ʊ�ͤλȤ����������������®�٤�����
  //  ͥ���١�(1)�ޤ�������
  //          (2)�������Τ�����®�٥٥��ȥ�������˶ᤤ��

  //  �֤��͡���ĤǤ������ʬ�������index�ˤ�������ʬ��shape_matrix���֤�
  //===============================================
  int find_solid_nearest_neighbor(double x,   double z, 
				  double v_x, double v_z,
				  int *i_x, int *i_z ,
				  int *i_x_particle, int *i_z_particle) ;

  // --- ���塧 (bare) Si ��Cl bond �������䤹��Cl neutral�ε����
  // ���ϡ�����Υ���ǥå����ʾ�δؿ��������
  //       �����Ψ
  // �֤��͡����夷�� -> TRUE
  bool settle_Cl_into_bareSi(int i_x, int i_z, 
			     double adsorption_probability) ;
  
  // --- ���塧 Si ��Cl bond �������䤷�����ˤʤä���æΥ
  // ���ϡ�����Υ���ǥå����ʾ�δؿ��������
  //  void settle_Cl_and_desorption(int i_x, int i_z) ;

  //***************************************************
  //***  02 Dec 2003  : modify the desorption algorithm
  //       (The above two functions are to be replaced)
  //***************************************************
  // --- have a Cl atom adsorbed
  void adsorb_Cl(int i_x, int i_z) ;
  void adsorb_Cl(int i_x, int i_z, 
		 double x, double z); // ��Ͽ��
  // --- just desorb Si atom thereat
  void desorb_Si(int i_x, int i_z) ;
  void desorb_Si(int i_x, int i_z, int dummy) ;  // ��Ͽ��

  void desorb_mask(int i_x, int i_z) ;//== �ϡ��ɥޥ����ο���

  //=====================================================
  //==   2005 / 01 / 18 chemical etching �δؿ����ɲ�
  //    ȿ���� SiCl4(s) -> SiCl4(g) �˽�������Cl neutral�����ͤ���������
  //    æΥ��������ʪ��flight�Ϲ�θ���ʤ���
  //    
  void desorb_Si_chemical_etch(int i_x, int i_z);

  //=====================================================
  //==   2005 / 08 / 02 SiClx(x = 1��4)��æΥ�δؿ����ɲ�
  //     SiClx ��æΥ�γ�Ψ��SiCl4�� x/4 �Ȥ���
  void desorb_SiClx(int i_x, int i_z) ;
  
  //***********************************
  // ***  05 Aug 2004 : deposition ���ɲ�
  // desorb_Si �εա�SPACE -> SOLID�ˤ�Ԥ�

  // ���Ѥ���Τ�Si���ҤΤߤ˸��ꤹ�뤿�ᡢ
  // ���ܤ������Υ���Υ���ǥå���������

  // SiClx ��������θ���뤿�ᡢCl�θ��ҿ�������
  // settle_Cl_into_bareSi �ؿ���Ʊ�ͤˡ������Ψ������Ȥ���
  // �֤��͡����夷�� -> TRUE
  bool deposit_Si(int i_deposit_x, int i_deposit_z, 
		  int  i_neighbor_x, int i_neighbor_z, 
		  int  n_Cl ,  int n_O,
		  double deposition_probability) ;
  
  //***********************************
  // ***  24 Aug 2004 : oxidation ���ɲ�
  // ���Ǥ�ȿ�����Ƥ�����ǿ���chlorination �ο��˱����ƽ������Ѥ���
  // 
  // 1: �ޤ���reaction probability �򥯥ꥢ�����Tuda ��
  //   ���ȡ����ǥ饸����ϳ�Ψ100���ȿ����
  // 2:  n_oxygen < 2 �Ǥ���С�ȿ������(n_oxygen++)
  // 3:  n_oxygen = 1 �ˤʤä��� -> n_Clbond > 2 �Ǥ����
  //     n_Clbond = 2 �Ȥ���
  // 4:  n_oxygen = 2 �ˤʤä��� -> n_Clbond = 0 �Ȥ���
  // �ʤĤޤꡢCl + O��2 ������Ķ���ʤ��褦�ˤ����
  
  // ���ϡ������ȿ���˳�Ψ ��
  // �֤��͡�ȿ������  -> true

  bool oxidation_Si(int i_x, int i_z, 
		    double reaction_probability ) ;

  //***************************************************
  //--- �����ϴ�Ϣ
  // Cl density distribution - depth
  // �����פ�Ȥ롣
  // ��Z���ξ夫��shape_matrix��ߤƤ�����������ʬ�ˤʤä��Ȥ�����
  // Cl ��ǻ�٤�ȤäƤ�����
  // -> �����ƥե�������ϡ��ե�����̾�˥����ॹ�ƥåפ�ޤ�뤿�ᡢ
  // �����ॹ�ƥåפ������ͤȤ���
  // 
  void output_Cl_density_depth(char file_name[],
			       int  time_step );

  //--- �����ե�������ɤ߹���
  //    ���ϡ�file name(Si, Cl)
  void input_profile(char file_Si[],char file_Cl[]) ;

  //***************************************************
  // --- "­��"��̵���ʤä�Si ���Ҥ������
  void remove_isolated_Si() ;

  //***************************************************
  // --- ���룱�ĤΥ���μ��ϤΣ������ʾ岼�����ˤ�Ĵ�١�
  // Four point calculation �ˤ�ä�ɽ�̤�ˡ����������ꤹ�롣
  // ˡ���٥��ȥ�϶��� -> ���������������Ȥ���Τǡ�
  //�岼�����Υ����Si���Ҥ�¸�ߤ����顢���줾��
  // (x,z) = (0,-1)(0,+1)(-1,0)(1,0) �νŤߤ�ä��Ƥ��Ф褤��

  // �ǽ�˥�����ֹ��Ƚ��ʾ岼�����˥��뤬¸�ߤ��뤫����
  // ���ϡ������ֹ�
  void get_surfacenormal(int i_x, int i_z ) ;

  // === Si ���Ҥ�æΥ�������ˡ����ϤΣ��ĤΥ�����Ф���
  // get_surfacenormal ��¹Ԥ���
  void get_surfacenormal_around(int i_x, int i_z ) ;

  // === ��������Ф��ƹԤ���constructor��Ǽ¹ԡ�
  void get_surfacenormal_all() ;

  // === �������®����ʬ���顢���ͳѤ����(rad)
  // ���ϡ������ֹ桿�������®��
  double get_incident_angle(int i_x, int i_z,
			    double v_x, double v_z) ;

protected:
  // -- ���ϡ��������� (X, Y) / �ɤ�Ĥ֤����κ�����
  void paint_fill(int start_point_x, int start_point_y,
		  int n_maxpoint ) ;


};

//==================================
//=== Y(��) : etching yield�γ��ٰ�¸�� (Hwang, Giapis)
//    |��| <  ��/4 �� 1,
//    |��| >= ��/4 �� cos��/cos(��/4)

//=== ɽ�̤λ����ˤ��Yield ���Ѳ���Ƴ������
// -> �����ͤϤ��γʻҤ˴ޤޤ����Ǥο� (n_oxygen)
//   �Ȥ��ơ�������εտ���ݤ��롣

double etch_yield_angle(double theta, int n_oxygen ) ;

//== (3rd paper)���Ǥ����夹�뤳�Ȥˤ������԰�����ɤ����ᡢ
//   ���å��󥰼�Ψ��ʬ�������롣nearest-neighbor�Υ���˴ޤޤ����Ǥο����פ���
//   6 �ǳ�ä�coverage���Ȥ��롣
//   ��Ψ�� Y = �� Y(SiO2) + (1 - ��)Y(Si)
double etch_yield_disperse_oxidation(double theta, int n_oxygen ) ;

//==================================
//=== flux �򥫥���Ȥ��뤿��Υ��饹

class Shape_counter_class
{
public:
  int *ctr[N_CELL_X]; int *p_ctr; 

public:
  //== constructor�ʥ���γ��ݤ��̤δؿ��ǹԤ���
  Shape_counter_class()
  {}
  ~Shape_counter_class()
  {
    delete [] p_ctr ;
  }

  // -- ������� & �����󥿽���� --
  void allocate_memory()
  {
    p_ctr = new int[N_CELL_X * N_CELL_Z] ;
    for(int i_x = 0; i_x < N_CELL_X ; i_x++)
      ctr[i_x] = p_ctr + ( i_x * N_CELL_Z ) ;
	
    for(int i_z = 0; i_z < N_CELL_Z ; i_z++)
      {
	for(int i_x = 0; i_x < N_CELL_X ; i_x++)
	  ctr[i_x][i_z] = 0 ;
      }
  }   
  // -- counting �ʥ����󥿤ˣ���ä����
  // ���ϡ�cell index (i_x, i_z)
  void count(int i_x, int i_z )
  {
    ctr[i_x][i_z]++ ;
  }
};

#define _SHAPE_TRIM_DEFINED_
#endif
