
#ifndef _INTEGRATION_ANGLE_H_

// --- Stillinger-Weberポテンシャルが 
//     minimumとなる原子間距離 
const double LENGTH_SW_MINIMUM = 2.1e-10 ;

//=== 衝突による散乱角の計算 ===
//    重心座標における跳ね返り角を求める
//    Ion stopping .. に記載している
//    積分式を解く

double theta_integrate(double p, double Ec) ;


//=== 上と同じく散乱角の計算 ===
//    積分の間隔が ∞なので、
//    non-equally spaced のΔr を用いる。
//    また、台形公式も用いる(Newton-Cotes)
double theta_integrate_NC(double p, double Ec) ;


#define _INTEGRATION_ANGLE_H_
#endif
