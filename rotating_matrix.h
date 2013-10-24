//
// 衝突前の速度方向を軸とした相対的な球面座標
// -> laboratory 座標への変換を行う

/*----------------------------------
   座標の取り方（速度には球面座標を用いる）
      y
   \  |  /
    \ | /
     \|/) ψ
      \--------------------->z
      |\ 
      | \     z 軸回り:θ
          x
*/

//衝突前の速度方向を(v, θ, ψ) とする。
// ベクトル (-sinθ, cosθ, 0) を保つ変換ということから、

// a12 = (a11 - 1)tanθ
// a22 = 1 + a21 tanθ
// a31 = a31 tanθ

// また、(cosθ, sinθ, 0)と、これを変換したベクトル
// A (cosθ, sinθ, 0)　のなす角がψであるので、
//
// (a11 cosθ+ a12 sinθ)cosθ
// + (a21 cosθ+ a22 sinθ)sinθ= cosψ

// また、Aは正規直交行列であることを利用して、各係数を求める。

// --------一箇所sqrt の正負に任意性があるので、
// TRUE/FALSE で条件分岐する。
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
