
#ifndef _INTEGRATION_ANGLE_H_

// --- Stillinger-Weber$B%]%F%s%7%c%k$,(B 
//     minimum$B$H$J$k86;R4V5wN%(B 
const double LENGTH_SW_MINIMUM = 2.1e-10 ;

//=== $B>WFM$K$h$k;6Mp3Q$N7W;;(B ===
//    $B=E?4:BI8$K$*$1$kD7$MJV$j3Q$r5a$a$k(B
//    Ion stopping .. $B$K5-:\$7$F$$$k(B
//    $B@QJ,<0$r2r$/(B

double theta_integrate(double p, double Ec) ;


//=== $B>e$HF1$8$/;6Mp3Q$N7W;;(B ===
//    $B@QJ,$N4V3V$,(B $B!g$J$N$G!"(B
//    non-equally spaced $B$N&$(Br $B$rMQ$$$k!#(B
//    $B$^$?!"Bf7A8x<0$bMQ$$$k(B(Newton-Cotes)
double theta_integrate_NC(double p, double Ec) ;


#define _INTEGRATION_ANGLE_H_
#endif
