// random_number.cc

#include <stdlib.h>                           
#include "random_number.h"




static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializing the array with a NONZERO seed */
void sgenrand(unsigned long seed)
{
  /* setting initial seeds to mt[N] using         */
  /* the generator Line 25 of Table 1 in          */
  /* [KNUTH 1981, The Art of Computer Programming */
  /*    Vol. 2 (2nd Ed.), pp102]                  */
  mt[0]= seed & 0xffffffff;
  for (mti=1; mti<N; mti++)
    mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

/* generating reals */
/* unsigned long */ /* for integer generation */
double RAN0()
{
  /*
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    // mag01[x] = x * MATRIX_A  for x=0,1 
    
    if (mti >= N)
    { // generate N words at one time 
    int kk;
      
    if (mti == N+1)   // if sgenrand() has not been called, 
    sgenrand(4357); // a default initial seed is used   
    
    for (kk=0;kk<N-M;kk++) {
    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
    mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    for (;kk<N-1;kk++) {
    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
    mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
    }
    y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
    mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];
    
    mti = 0;
    }
    
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
  */

  return (double(rand()) / RAND_MAX ) ;

  //return ( ((double)y++) * 2.3283064370807974e-10 ); /* translation from reals: [0,1)-interval to reals: (0,1]-interval*/
  
  /* return ( ((double)y++) * 2.3283064365386963e-10 ); */ /* reals: [0,1)-interval */ 
  /* return y; */ /* for integer generation */

}

// ----------------------------
// Žg‚¢•û

/* this main() outputs first 1000 generated numbers  */
/*
main()
{ 
  int j;
  
  sgenrand(4357); // any nonzero integer can be used as a seed 
  for (j=0; j<1000; j++) 
    {
      // printf("%10.8f ", genrand());
      printf("%10.8f \n", RAN0());
      //if (j%8==7) printf("\n");
    }
  printf("\n");
}
*/

/* —””­¶ƒ‹[ƒ`ƒ“ (From Nomura's random.cpp) */	
// 0 - 1 ‚Ì—”‚ð”­¶‚³‚¹‚é
/*   
double RAN0(long *idum)

{  // idum==2 ‚ª‚¢‚¢
  const int IA=16807;
  const long IM=2147483647;
  const long IQ=127773;
  const int IR=2836;
  const long MASK=123459876;
  const double AM=1.0/IM;
  
  long k;
  double ans;
  
  *idum ^= MASK;
  k = (*idum)/IQ;
  *idum = IA*(*idum-k*IQ)-IR*k;
  
  if (*idum < 0) 
    *idum +=IM;
    
  ans=AM*(*idum);
  *idum ^= MASK;
    
  return ans;

}
*/
