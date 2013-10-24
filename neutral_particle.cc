
#include <iostream>
#include <math.h>
#include "neutral_particle.h"
#include "random_number.h"
#include "atom_struct.h"
#include "shape_trim.h"

Neutral_class::Neutral_class() 
{
}

Neutral_class::Neutral_class(double mass_input , int flag_boundary_input) 
  : Particle_class(mass_input, flag_boundary_input) 
{
  //mass = mass_input ; 
}
Neutral_class::~Neutral_class() 
{}


void Neutral_class::inject_from_top( double v )
{
  pos_v.x = RAN0() * SYSTEM_WIDTH_X ; 
  pos_v.y = 0.0 ;
  pos_v.z = 0.0 ;

  // acos(1-2R)
  //double theta = RAN0() * PI ;// obsolete: 
  double theta = acos(1.0 - 2.0 * RAN0()) ;

  pos_v.v_x = v * cos(theta) ;
  pos_v.v_y = 0.0 ;
  pos_v.v_z = v * sin(theta) ;

  pos_v.v_r     = v ;
  pos_v.v_theta = 0.0 ;
  pos_v.v_psi   = sin(theta) ;

	//	cout << pos_v.v_x << "\t" << pos_v.v_y << "\t"
     	//	<< pos_v.v_z << "\n" ;
  i_exception_move = 0 ;//== $BNc30=hM}MQ%+%&%s%?$N=i4|2=(B ==
}
void Neutral_class::
inject_from_right_side(double v)
{
  pos_v.x = SYSTEM_WIDTH_X ; 
  pos_v.y = 0.0 ;
  pos_v.z = RAN0() * SYSTEM_HEIGHT_Z ;

  // acos(1-2R)
  //double theta = RAN0() * PI ;// obsolete: 
  double theta = acos(1.0 - 2.0 * RAN0()) ;

  pos_v.v_x = - v * sin(theta) ;
  pos_v.v_y =   0.0 ;
  pos_v.v_z =   v * cos(theta) ;

  pos_v.v_r     = v ;
  pos_v.v_theta = 0.0 ;
  pos_v.v_psi   = sin(theta) ;
}
//=================================
// $BCf@-N3;R$NF~<M$+$iJ|=P!"5[Ce$^$G(B
void Neutral_class::
all_process(class Shape_trim_class *Shape,
	    class Shape_counter_class *Neutral_counter,
	    bool   flag_inject_from_side ,
	    bool   flag_flux_count ,
	    double incident_energy,
	    bool   flag_chemical_etch,
	    double yield_chemical )
{
  //=== $B%;%k%$%s%G%C%/%9(B ===
  int i_x, i_z, i_x_particle, i_z_particle ;
  //--- $BF~<M(B ---
  inject_from_top
    ( sqrt(2.0 * incident_energy * 
	   Q_ELEMENTAL / mass) ) ; 

  i_exception_move = 0 ;
  
  //== $B2#J}8~$+$i$NF~<M(B ==
  if(flag_inject_from_side == true)
    {
      inject_from_right_side
	( sqrt(2.0 * incident_energy * 
	       Q_ELEMENTAL / mass) ) ; 
      //== $BF~<M$N0LCV$,6u4V$GL5$1$l$P(B i_exception_move $B$KBg$-$$(B
      //   $BCM$rBeF~$7$F(B while loop $B$r(Bskip $B$9$k(B
      //   put_shape.. $B4X?t$ONN0h1&C<$G(B0$B$rJV$9$N$G!"(BL/2$B$r0z$/(B
      if(Shape->put_shape_matrix(pos_v.x - 0.5 * CELL_SIZE, pos_v.z ,
				     &i_x, &i_z ) != SHAPE_SPACE )
	i_exception_move = N_CELL_Z * 200 + 1 ;
    }
 
  while(1) 
    {
      //-- $BNN0hFbIt$K$"$^$jD9$/$H$I$^$C$F$$$k>l9g!"(B
      //   $BNc30=hM}!JN3;R$r>C5n!K(B
      i_exception_move++ ;
      if( i_exception_move > N_CELL_Z * 200 )
	break ;
		
      // -- $B$^$:JB?J1?F0(B ----
      move_trans(FREE_FLIGHT_PATH_Si) ;

      if (flag_inside == false) // -- $B30$K=P$?!)(B
	break ;
	     
      //--- $B%;%k$,6uGr$G$J$$!)(B ---
      //if(Shape.put_shape_matrix(pos_v.x, pos_v.z ,
      //			&i_x, &i_z ) != SHAPE_SPACE)
      if( Shape->find_solid_nearest_neighbor
	  (pos_v.x,    pos_v.z ,
	   pos_v.v_x,  pos_v.v_z ,
	   &i_x, &i_z, &i_x_particle, &i_z_particle ) != SHAPE_SPACE )
	{
	  //== flux counting 
	  if(flag_flux_count == true )
	    Neutral_counter->count(i_x, i_z ) ;
	  
	  //---- $B2=3XE*%(%C%A%s%0(B ----
	  // YIELD_CHEMICAL$B$N3NN($GC&N%$,5/$3$k(B 
	  if( RAN0() < yield_chemical &&
	      flag_chemical_etch == true  )
	    {
	      Shape->desorb_Si_chemical_etch(i_x, i_z);
	    }
	  else //---- $B5[Ce(B -----
	    {
	      flag_adsorption = 
		Shape->settle_Cl_into_bareSi(i_x, i_z,  // $B5[Ce$7$?!)(B
					    ADSORPTION_PROBABILITY );
	      if(flag_adsorption == true )
		break;   // $B5[Ce$7$?$i%k!<%W$r=P$k(B
	    }
	  // $BMpH?<M(B : 
	  // $BL58B%k!<%W$r$5$1$k$?$a$K!"%+%&%s%?$rMQ0U$9$k!#(B
	  // $BN3;R$r2>$K?J$^$;$F$_$F!"(B(pos_v_togo)
	  // $B$b$7?J$s$@@h$,6uGrItJ,$G$J$1$l$P!"85$KLa$C$F(B
	  //     $B$d$jD>$7(B(50$B2s$K$J$C$?$i$b$&$d$a$k(B)
	  i_exception_reflect = 0;
	  Shape->get_surfacenormal(i_x, i_z) ;
	reflection:
	  if( i_exception_reflect > 50 )
	    break ;
	  
	  random_reflection //(FALSE,0.0,0.0);//$BMpH?<M(B
	    (true,//<-cosine diffusion $B$rF3F~$9$k$H:81&HsBP>N$K$J$k(B<- this bug was fixed (29Aug2004)
	     Shape->surfacenormal_x[i_x][i_z],
	     Shape->surfacenormal_z[i_x][i_z]) ;
		
	  move_trans_togo(FREE_FLIGHT_PATH_Si) ;
	  
	  if(Shape->put_shape_matrix(pos_v_togo.x, pos_v_togo.z ,
				     &i_x, &i_z ) != SHAPE_SPACE )
	    {
	      i_exception_reflect++ ;
	      goto reflection ;
	    }
	  // cout << Shape.surfacenormal_x[i_x][i_z] * 10000.0 
	  //	<< "\t" << Shape.surfacenormal_z[i_x][i_z] * 10000.0
	  //	<<  "\t" << Cl_neutral.pos_v.v_x << "\t" << Cl_neutral.pos_v.v_z << endl ;
	}
    }//=== while loop $B=*N;(B
  
}
