// shape_trim.h
//    -> 形状進展の状態を表すためのクラス。

// １つのセルで１つのSi 原子を表すものとする。（ただし２次元）
// shate_state クラスから、trim で必要な要素のみを抜き出して用いる。



#ifndef _SHAPE_TRIM_DEFINED_

const int SHAPE_SPACE = 0 ;
const int SHAPE_Si    = 10 ;
const int HARD_MASK   = 5 ;

const int CL_DEPTH_RANGE = 20 ;
#include "atom_struct.h"

class Shape_trim_class
{
public:
  // メインのmatrix
  int *shape_matrix[N_CELL_X]; int *p_shape_matrix; 
  
  // --- それぞれのSi 原子に結合している Cl の数
  int *n_Clbond[N_CELL_X];     int *p_n_Clbond  ;

  // --- 同じく Si原子に結合している Oラジカルの数
  int *n_oxygen[N_CELL_X];     int *p_n_oxygen  ;
  
  // *** 表面における法線方向：

  // * すべてのセルにおいて計算する。（表面で無い場合は０ベクトル）
  // * 空間 -> 固体内部の方向
  // * Four-point calculation によって計算

  double *surfacenormal_x[N_CELL_X] ; double *p_surfacenormal_x ;
  double *surfacenormal_z[N_CELL_X] ; double *p_surfacenormal_z ;

  //---- 脱離したSi原子のカウンタ／ハードマスクのカウンタ
  //     (desorb_Si/desorb_mask関数)
  int cntr_desorbed_Si ; 
  int cntr_desorbed_mask ;
  
  //---- 堆積したSi原子のカウンタ
  //     (deposit_Si関数)
  int cntr_deposited_Si ;
  
  // -- domain の下まで到達したら、 TRUEとする。
  bool flag_get_bottom ;

  //========================================
  // 吸着した位置を記録しておくためのメンバ。
  // md2d と同様に、クラス内のクラスにする
  class Sticking_Cl_class
  {
  public:
    double x ; // 位置
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
  
  // コンストラクタ 
  Shape_trim_class() ;
  Shape_trim_class(bool flag_open_space) ;
  ~Shape_trim_class() ;

  // -- 座標を入力したら shape_matrix を返す関数
  // -- 負の値、もしくは範囲外の値を入力した場合は０を返す

  //===============================================
  //===  セルのインデックスをポインタで返す
  //===============================================
  int put_shape_matrix(double x, double z, 
		       int *i_x, int *i_z ) ;

  //===============================================
  //    ある格子点において、その点かもしくは隣接点に
  //     固体部分(SPACEでない点) があれば、そのindexをポインタで返す

  //    追加：粒子の位置(x, z)のセルのindexも返す
  //     put_shape_matrix関数と同様の使い方だが、イオンの速度も入力
  //  優先度：(1)まずその点
  //          (2)隣接点のうち、速度ベクトルの方向に近い点

  //  返り値：一つでも固体部分があればindexにあたる部分のshape_matrixを返す
  //===============================================
  int find_solid_nearest_neighbor(double x,   double z, 
				  double v_x, double v_z,
				  int *i_x, int *i_z ,
				  int *i_x_particle, int *i_z_particle) ;

  // --- 定着： (bare) Si のCl bond を一つ増やす（Cl neutralの吸着）
  // 入力：セルのインデックス（上の関数で得る）
  //       吸着確率
  // 返り値：吸着した -> TRUE
  bool settle_Cl_into_bareSi(int i_x, int i_z, 
			     double adsorption_probability) ;
  
  // --- 定着： Si のCl bond を一つ増やし、４になったら脱離
  // 入力：セルのインデックス（上の関数で得る）
  //  void settle_Cl_and_desorption(int i_x, int i_z) ;

  //***************************************************
  //***  02 Dec 2003  : modify the desorption algorithm
  //       (The above two functions are to be replaced)
  //***************************************************
  // --- have a Cl atom adsorbed
  void adsorb_Cl(int i_x, int i_z) ;
  void adsorb_Cl(int i_x, int i_z, 
		 double x, double z); // 記録用
  // --- just desorb Si atom thereat
  void desorb_Si(int i_x, int i_z) ;
  void desorb_Si(int i_x, int i_z, int dummy) ;  // 記録用

  void desorb_mask(int i_x, int i_z) ;//== ハードマスクの浸食

  //=====================================================
  //==   2005 / 01 / 18 chemical etching の関数を追加
  //    反応式 SiCl4(s) -> SiCl4(g) に従う。（Cl neutralが入射した時点）
  //    脱離した生成物のflightは考慮しない。
  //    
  void desorb_Si_chemical_etch(int i_x, int i_z);

  //=====================================================
  //==   2005 / 08 / 02 SiClx(x = 1〜4)の脱離の関数を追加
  //     SiClx の脱離の確率はSiCl4の x/4 とする
  void desorb_SiClx(int i_x, int i_z) ;
  
  //***********************************
  // ***  05 Aug 2004 : deposition を追加
  // desorb_Si の逆（SPACE -> SOLID）を行う

  // 堆積するのをSi原子のみに限定するため、
  // 隣接する固相のセルのインデックスも入力

  // SiClx の付着を考慮するため、Clの原子数も入力
  // settle_Cl_into_bareSi 関数と同様に、吸着確率も引数とする
  // 返り値：吸着した -> TRUE
  bool deposit_Si(int i_deposit_x, int i_deposit_z, 
		  int  i_neighbor_x, int i_neighbor_z, 
		  int  n_Cl ,  int n_O,
		  double deposition_probability) ;
  
  //***********************************
  // ***  24 Aug 2004 : oxidation を追加
  // すでに反応している酸素数、chlorination の数に応じて処理を変える
  // 
  // 1: まず、reaction probability をクリアする（Tuda に
  //   よると、酸素ラジカルは確率100％で反応）
  // 2:  n_oxygen < 2 であれば、反応する(n_oxygen++)
  // 3:  n_oxygen = 1 になった？ -> n_Clbond > 2 であれば
  //     n_Clbond = 2 とする
  // 4:  n_oxygen = 2 になった？ -> n_Clbond = 0 とする
  // （つまり、Cl + O×2 が４を超えないようにする）
  
  // 入力：吸着（反応）確率 ／
  // 返り値：反応した  -> true

  bool oxidation_Si(int i_x, int i_z, 
		    double reaction_probability ) ;

  //***************************************************
  //--- 入出力関連
  // Cl density distribution - depth
  // の統計をとる。
  // 各Z軸の上からshape_matrixをみていき、固体部分になったところから
  // Cl の濃度をとっていく。
  // -> そしてファイル出力。ファイル名にタイムステップを含めるため、
  // タイムステップを入力値とする
  // 
  void output_Cl_density_depth(char file_name[],
			       int  time_step );

  //--- 形状ファイルの読み込み
  //    入力：file name(Si, Cl)
  void input_profile(char file_Si[],char file_Cl[]) ;

  //***************************************************
  // --- "足場"が無くなったSi 原子を取り除く
  void remove_isolated_Si() ;

  //***************************************************
  // --- ある１つのセルの周囲の４方向（上下左右）を調べ、
  // Four point calculation によって表面の法線方向を決定する。
  // 法線ベクトルは空間 -> 固体内部の方向とするので、
  //上下左右のセルにSi原子が存在したら、それぞれ
  // (x,z) = (0,-1)(0,+1)(-1,0)(1,0) の重みを加えてやればよい。

  // 最初にセルの番号を判定（上下左右にセルが存在するか？）
  // 入力：セル番号
  void get_surfacenormal(int i_x, int i_z ) ;

  // === Si 原子が脱離した時に、周囲の４つのセルに対して
  // get_surfacenormal を実行する
  void get_surfacenormal_around(int i_x, int i_z ) ;

  // === 全セルに対して行う（constructor内で実行）
  void get_surfacenormal_all() ;

  // === イオンの速度成分から、入射角を求める(rad)
  // 入力：セル番号／イオンの速度
  double get_incident_angle(int i_x, int i_z,
			    double v_x, double v_z) ;

protected:
  // -- 入力：開始点の (X, Y) / 塗りつぶし点の最大値
  void paint_fill(int start_point_x, int start_point_y,
		  int n_maxpoint ) ;


};

//==================================
//=== Y(θ) : etching yieldの角度依存性 (Hwang, Giapis)
//    |θ| <  π/4 で 1,
//    |θ| >= π/4 で cosθ/cos(π/4)

//=== 表面の酸化によるYield の変化も導入する
// -> 入力値はその格子に含まれる酸素の数 (n_oxygen)
//   として、選択比の逆数を掛ける。

double etch_yield_angle(double theta, int n_oxygen ) ;

//== (3rd paper)酸素が吸着することによる形状不安定を防ぐため、
//   エッチング収率を分散させる。nearest-neighborのセルに含まれる酸素の数を合計し、
//   6 で割ってcoverageΘとする。
//   収率は Y = Θ Y(SiO2) + (1 - Θ)Y(Si)
double etch_yield_disperse_oxidation(double theta, int n_oxygen ) ;

//==================================
//=== flux をカウントするためのクラス

class Shape_counter_class
{
public:
  int *ctr[N_CELL_X]; int *p_ctr; 

public:
  //== constructor（メモリの確保は別の関数で行う）
  Shape_counter_class()
  {}
  ~Shape_counter_class()
  {
    delete [] p_ctr ;
  }

  // -- メモリ確保 & カウンタ初期化 --
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
  // -- counting （カウンタに１を加える）
  // 入力：cell index (i_x, i_z)
  void count(int i_x, int i_z )
  {
    ctr[i_x][i_z]++ ;
  }
};

#define _SHAPE_TRIM_DEFINED_
#endif
