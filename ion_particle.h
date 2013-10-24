// ion_particle.h

// $B%$%*%s$N%/%i%9!'(BParticle class $B$r7Q>5(B

// 19Oct2004
// SiCl4 $B$NC&N%$r9MN8$9$k$?$a$N%U%i%0$r$3$3$GF3F~(B


#ifndef  _ION_PARTICLE_H_DEFINED_

#include <iostream>
#include <fstream>
#include "particle.h"
#include "shape_trim.h"

class Ion_class : public Particle_class
{
  
public:
  //== $B%$%*%s$N>WFM$K$h$C$FI=LL$N(BSi$B86;R$,C&N%$7$?$+$r<($9%+%&%s%?(B ==
  // $B!J%$%*%s$NF~<M$N$?$S$K(Bfalse$B$KLa$9!K(B
  bool  flag_desorb_Si ;
  //== $BC&N%$7$?@8@.J*Cf$N;@AG$N?t(B (SiClxOy$B$N(By)
  int   n_oxy_desorb_Si ; 
  
  //== $BC&N%$,5/$3$C$?;~$N%$%*%s$N:BI8(B: $B9=B$BN$G5-O?$7$F$*$/!JB.EY$OIT;HMQ!K(B
  // $B!J(BSiCl4$B$N(Bflight$B$,$3$3$+$i;O$^$k$H$9$k!K(B
  struct Particle_location_velocity_struct position_at_desorption ;

protected:
  //========================================================
  //    $B%$%*%s$NA0J};6Mp$K4X$9$k5-O?$r9T$&(B
  bool flag_reflection ; // $BA0J};6Mp$,5/$3$C$?$+$I$&$+$N(Bflag
  int  ctr_reflection  ; // $B%+%&%s%?(B
  std::ofstream velocity_file ; //== $B;6Mp8e$NB.EY$r5-O?$9$k$?$a$N%U%!%$%k(B

protected:
  //========================================================
  //    $B%$%*%s$,:GI=LL$K>WFM$7$?8e!"H?<M$9$k$+!?4pHD>e$K?/F~$7$F(Bstop$B$9$k$+(B
  //    $BH=CG$7!"(Byield$B$r6hJL$9$k>l9g$KMQ$$$kJQ?t(B
  //    $B$"$i$+$8$a>WFM;~$N0LCVEy$r5-O?$7!";v8eE*$K%(%C%A%s%0$N=hM}$r9T$&(B
  //==  $B>WFM!J%(%C%A%s%0!K0LCV!"%(%M%k%.!<$r5-O?$9$kJQ?t(B  ==
  int    i_x_etch, i_z_etch ;
  double energy_etch ;
  struct Particle_location_velocity_struct pos_v_etch ;

public:
  void  record_etch_position(int i_x, int i_z, double input_energy)
  {
    i_x_etch    = i_x ;          i_z_etch = i_z ;
    energy_etch = input_energy ; pos_v_etch = pos_v ; 
  }

public:
  //-- constructor 
  Ion_class() ;
  Ion_class(double mass_input, int flag_boundary_input,
	    int n_adf_input );
  ~Ion_class();

  //== $BC&N%%+%&%s%?$N%j%;%C%H(B ==
  //void reset_ctr_desorb_Si(){ ctr_desorb_Si = 0 ;}
  
  //== $BA0J};6Mp4XO"(B ==
  void set_flag_reflection(bool input_flag) //== flag$B$N@_Dj(B ==
  { flag_reflection = input_flag ; }
  
  void add_ctr_reflection()//== $B%+%&%s%?$N2C;;(B ==
  { if(flag_reflection == true) ctr_reflection++ ;  }

  int  put_ctr_reflection()//== $B%+%&%s%?$N=PNO(B ==
  { return  ctr_reflection ; }

  //=============================================
  //*********************************************************
  // -- microstructure/charged_particle::inject_from_top $B4X?t(B
  // $B$r!"2~NI$7$F;HMQ$9$k!#$3$3$G$O!"$^$:(Bangular distribution $B$N$_(B
  // $BF3F~$9$k!#$^$?!"EjF~$9$kN3;R$O#1$D$H$9$k!#(B
  void inject_iadf( //int n_injection, 
		   //	 double energy_array[],double energy_threshold,
		   //	 double energy_df[],
		   double v_z) ;  
  // double angle_array[],  double angle_df[] );


  //==============================================
  //== $BC&N%$,5/$3$C$?;~$K%U%i%0$rN)$F$F0LCV$r5-O?(B
  //   $BF~NO!'!J%(%C%A%s%0;~E@$N!K0LCV!?B.EY(B, $B;@AG$N?t(B
  void record_desorption(Particle_location_velocity_struct pos_v_recorded,
			 int input_n_oxygen) ;


  //==============================================
  //== $B%$%*%s$N>WFM$K$h$k%(%C%A%s%0$N=hM}(B
  //   $B%$%*%s$,H?<M$7$?>l9g$N%(%C%A%s%0<}N($O(B Er < Eth $B$G$"$l$P(B C($B"e(BEi - $B"e(BEth),
  //                                          Er > Eth $B$G$"$l$P(B C($B"e(BEi - $B"e(BEr)$B$H$9$k(B
  //   $BH?<M$GL5$$>l9g$O(B Er = 0 $B$rF~NO(B
  //   $BF~NO!'(BShape_trim class 
  //   $B!J%(%C%A%s%0;~E@$N!KF~<M%(%M%k%.!<(B Ei, $BH?<M%(%M%k%.!<(B Er,
  //   $B0LCV!?B.EY!"(B $B3J;RE@$N(Bindex (i_x, i_z)
  void ion_enhanced_etch(class  Shape_trim_class *Shape_trim,
			 double incident_energy, double reflected_energy,
			 Particle_location_velocity_struct pos_v_recorded,
			 int i_x, int i_z ) ;

  //== $B%O!<%I%^%9%/$N%9%Q%C%?%j%s%0(B
  //   $B<}N($O;@2=$5$l$?(BSi$B4pHD$N$b$N$K$5$i$K78?t$r$+$1$k(B
  void hardmask_sputter(class  Shape_trim_class *Shape_trim,
			double incident_energy, double reflected_energy,
			Particle_location_velocity_struct pos_v_recorded,
			int i_x, int i_z ) ;

  //==============================================
  //==   $B%Q%?!<%sI=LL$KE~C#$7$?$+$I$&$+$NH=Dj(B 
  //   ->$B%(%C%A%s%0!">WFM$N=hM}(B
  //  $BJV$jCM!'(Btrue: $BNN0h$r=P$k!"$b$7$/$O%$%*%s$,(Bstop$B$9$k(B
  //   $B!J(Bwhile loop $B$r=P$k!K(B
  bool impact_on_surface(class  Shape_trim_class *Shape_trim) ;

  //==============================================
  //==   $B>e$HF1$8$/>WFM$N=hM}!'A0J};6Mp$N$b4^$a$k(B

  //     1 $B%$%*%s$N$$$k3J;RE@$,(BSi$B$G$"$k$+!)(B 
  //        YES -> $B=>MhDL$j$N>WFM$N=hM}(B 
  //  NO -> 2  $BNY@\3J;RE@$,(BSi/Hard mask $B$G$"$k$+!)(B
  //        YES -> impact parameter$B$r@53N$K7W;;$9$k(B
  //               (collision_accutate)
   bool impact_scattering(class  Shape_trim_class *Shape_trim,
			  bool  flag_mask_erosion) ;

};


#define  _ION_PARTICLE_H_DEFINED_
#endif
