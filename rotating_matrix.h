//
// �ՓˑO�̑��x���������Ƃ������ΓI�ȋ��ʍ��W
// -> laboratory ���W�ւ̕ϊ����s��

/*----------------------------------
   ���W�̎����i���x�ɂ͋��ʍ��W��p����j
      y
   \  |  /
    \ | /
     \|/) ��
      \--------------------->z
      |\ 
      | \     z �����:��
          x
*/

//�ՓˑO�̑��x������(v, ��, ��) �Ƃ���B
// �x�N�g�� (-sin��, cos��, 0) ��ۂϊ��Ƃ������Ƃ���A

// a12 = (a11 - 1)tan��
// a22 = 1 + a21 tan��
// a31 = a31 tan��

// �܂��A(cos��, sin��, 0)�ƁA�����ϊ������x�N�g��
// A (cos��, sin��, 0)�@�̂Ȃ��p���Ղł���̂ŁA
//
// (a11 cos��+ a12 sin��)cos��
// + (a21 cos��+ a22 sin��)sin��= cos��

// �܂��AA�͐��K�����s��ł��邱�Ƃ𗘗p���āA�e�W�������߂�B

// --------��ӏ�sqrt �̐����ɔC�Ӑ�������̂ŁA
// TRUE/FALSE �ŏ������򂷂�B
void rotating_matrix(double theta , double psi,
                     bool   sign_flag , 
                     double *a11, double *a12, double *a13,
                     double *a21, double *a22, double *a23,
                     double *a31, double *a32, double *a33 ) ;


//---- dust box
/*

 a11 = ( cos(Cl_incident.v_psi) +
                  sin(Cl_incident.v_theta) * sin(Cl_incident.v_theta) +
                  sin(Cl_incident.v_theta) * sin(Cl_incident.v_theta) *
                  tan(Cl_incident.v_theta) * tan(Cl_incident.v_theta) )
            /   ( 1.0 +
                  sin(Cl_incident.v_theta) * sin(Cl_incident.v_theta) +
                  sin(Cl_incident.v_theta) * sin(Cl_incident.v_theta) *
                  tan(Cl_incident.v_theta) * tan(Cl_incident.v_theta) ) ;

          a21 = (a11 - 1.0) * tan(Cl_incident.v_theta) ;

          a31 = - sqrt(1.0 - a11 * a11 - a21 * a21) ;
          //--
          a12 = (a11 - 1.0) * tan(Cl_incident.v_theta) ;

          a22 = 1.0 + a21 * tan(Cl_incident.v_theta) ;

          a32 = a31 * tan(Cl_incident.v_theta) ;
          //--
          a13 = sqrt(1.0 - a11 * a11 - a12 * a12)  ;
          a23 = - sqrt(1.0 - a21 * a21 - a22 * a22)  ;
          a33 = sqrt(1.0 - a31 * a31 - a32 * a32)  ;
*/
