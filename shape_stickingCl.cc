// shape_stickingCl.cc

#include <stddef.h>
#include "shape_trim.h"


// --- 記録用に位置(x,z)も入力値にするような関数も作成
// ---> (shape_stickingCl.cc)
void Shape_trim_class::adsorb_Cl(int i_x, int i_z, 
				 double x, double z) 
{ 
  adsorb_Cl(i_x, i_z) ;

  //================================
  //______________________________
  Sticking_Cl_class *new_ptr ;
  new_ptr           = new Sticking_Cl_class ;
  new_ptr->next_ptr = first_ptr ;
  first_ptr         = new_ptr ;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  new_ptr->x = x ;
  new_ptr->z = z ;

  n_sticking++ ;
}

void Shape_trim_class::desorb_Si(int i_x, int i_z, int dummy) 
{
  desorb_Si(i_x, i_z) ;

  //--- 全てのsticking_Cl に対して、i_x, i_z のセル内にいるかチェック
  //______________________________
  Sticking_Cl_class *tmp_ptr ;  
  Sticking_Cl_class *tmp_ptr2 ; Sticking_Cl_class *tmp_first_ptr ;

  tmp_ptr = first_ptr ;
  tmp_first_ptr = NULL ; //-- initialize

  if(i_x == int( tmp_ptr->x / CELL_SIZE ) &&
     i_z == int( tmp_ptr->z / CELL_SIZE ) &&
     tmp_ptr->next_ptr != NULL )
    {
      tmp_first_ptr  = tmp_ptr ;
      first_ptr      = tmp_ptr->next_ptr ;
    }
  
  while(1){
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(i_x == int( tmp_ptr->next_ptr->x / CELL_SIZE ) &&
       i_z == int( tmp_ptr->next_ptr->z / CELL_SIZE ) &&
       tmp_ptr->next_ptr != NULL )
      {
	//-- 消去する
	tmp_ptr2 = tmp_ptr->next_ptr ;
	tmp_ptr->next_ptr = tmp_ptr->next_ptr->next_ptr ;

	delete tmp_ptr2 ;
      }
    
    //______________________________
    tmp_ptr = tmp_ptr->next_ptr ;
    if(tmp_ptr == NULL){  break ; }
  }//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  delete tmp_first_ptr ;
}
