
#ifndef _HEADER_DEPOSITION_INCLUDED_

#include "atom_struct.h"

//== �ޥ��������etch products/���Ǥ��ǥ������뤫��

const bool FLAG_SiCl4_INCLUDED = true;//false;// <- basically true

const bool FLAG_SiCl2_INCLUDED = false;//true;//

const bool FLAG_OXYGEN_INCLUDED = false;//true;//

const bool FLAG_SiOCl2_INCLUDED = false;//true;//


//=== SiCl2 ��flux : Tuda [1996] �ˤ��� 
// ��p/��i = 0, 0.1, 0.5
// �Ȥ��Ƥ��롣

// -> ����������Υ����󥿤ˤ���ΤϤ��Τޤޤǡ�
// i_ion % [flux ratio] == 0 �λ��� SiCl2 �����ͤ���
// ��flux��Ͼ�εտ� ��i/��p �������

//const int ION_SiCl2_RATIO = 2;//1;//10 ;


// == ������˽Ф����å�������ʪ��SiCl2�˲�Υ����
//    ���Ĥ���äƤ����Ψ

const  double  SiCl2_REINCIDENCE_PROBABILITY = 0.0;//0.5 ;//0.2;//0.05;//

// @@@@@@  SiCl4 �ε����Ψ  @@@@@@@
// Tuda, 1996 JVST B �Ǥϡ��ºݤε����Ψ��
//   Sp <= 0.002 �ȵ��ܤ���Ƥ��뤬��
// ����礭���ͤⲾ�ꤷ�����β��Ϸ�̤��󼨤���Ƥ���
const  double  SiCl4_DEPOSITION_PROBABILITY = 0.02;//0.05;//0.1;//0.5;//0.2;//0.002;//0.01;//0.005;//0.001;

// @@@@@@  SiCl2 �ε����Ψ  @@@@@@@
// Tuda, 1996 JVST B �Ǥϡ������ΨSp�� 0��1 ��
// ������Ѳ������롢�ȵ��Ҥ���Ƥ��롣
//  Fig. 7 �Ǥ� Sp = 0.1 �Ȥ��ơ� ��p/��i = 0, 0.1, 0.5
// �Ȥ��Ƥ��롣 
const  double  SiCl2_DEPOSITION_PROBABILITY = 0.1;//0.0;//0.5;//0.2;//


//=== O �ʻ��ǡˤ�flux: Tuda [1996] �ˤ��� 
// ��o/��i = 0, 0.02, 0.05

const  double OXYGEN_ION_RATIO = 0.0;//1.0;//2.0;//0.5;//5.0;//0.1;//0.05;
//20.0;//50.0;//100.0;//30.0;//

// ��SiO/��i = 0, SiOCl2 �ε����Ψ
//const  double SiO_ION_RATIO = 0.0;//0.05 ;

const  double SiOCl2_REINCIDENCE_PROBABILITY = 0.0;//0.5;//0.2;//0.05;//
const  double SiOCl2_DEPOSITION_PROBABILITY  = 0.1;//0.5;//0.2; 

#define _HEADER_DEPOSITION_INCLUDED_
#endif
