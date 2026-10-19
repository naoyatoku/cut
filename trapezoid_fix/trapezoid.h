#ifndef	__TRAPEZOID_H__
#define	__TRAPEZOID_H__
//台形駆動関連のサブルーチンをまとめてみます。
//
#include <cmath>
#include <type_traits>

#include "CommonModule.h"
#include "accurized.h"
#include "tuple.h"


//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
//	加減速度に関する構造体です。（pana製のコードに合わせた形なので、ここに置いておきます。）
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
struct acc
{
	accurized<>	a1u;		//速度精度倍率が考慮された加速度とします。
	accurized<>	a1d;		//		"
	accurized<>	a3d;		//		"
	bool empty(void)const{
		return (a1u.raw() == 0)&&(a1d.raw()==0)&& (a3d.raw()==0);
	}
	//全部に共通する精度倍率を一度に参照できるようにしておきます。
	long accuracy(void) const {
		_Assert( (a1u.accuracy() == a1d.accuracy()) && (a1u.accuracy() == a3d.accuracy()) , "acc::accuracy() accuracy differnt");
		return a1u.accuracy();
	}
	//なんかいっぺんにセットできるように
	void set(const accurized<>&maxspd , long a1u_ms , long a3d_ms , long a1d_ms=0){
		if(!a1d_ms){	a1d_ms = a3d_ms;	}		//省略されている場合に
		_Assert(a1d_ms>0 && a1d_ms >0 && a3d_ms>0 , "acc_set() : accel time constant is 0 (%d,%d,%d,)" , a1u_ms , a1d_ms ,a3d_ms );
		a1u.set( (long long)maxspd.raw()/a1u_ms	,	maxspd.accuracy());
		a1d.set( (long long)maxspd.raw()/a1d_ms	,	maxspd.accuracy());
		a3d.set( (long long)maxspd.raw()/a3d_ms	,	maxspd.accuracy());
	}
	acc(){;}	//カラが帰るはずですが
	acc&operator=(const acc& a){	//代入
		a1u = a.a1u;
		a1d = a.a1d;
		a3d = a.a3d;
		return *this;
	}
};


//==========================================================================================================
//	三角関数
//==========================================================================================================

/*
    long N;      // 全時間 [tick] (整数)
	F    n_a;    // 加速時間 (= 減速時間) [tick]
    F    n_c;    // 等速時間 [tick]
    F    v_pk;   // 実際の最高速度 [pulse/tick]
    F    a_t;    // 実際の加速度   [pulse/tick^2]
*/
//----------------------------------------------------------------------------------------------------------------------------------
//  台形速度プロファイル  (C++11 対応 / 1ms 割り込み向け)
//
//  typename Num : 距離 L・速度 v・加速度 a の型 (整数でも可)
//  F            : 小数が必要な量の型 (Num が浮動小数点なら Num、整数なら double)
//
//  【使い方】
//    construct(L, v, a)            L[pulse], v[pulse/s], a[pulse/s^2] で指定。(dt は 1tick の秒数、既定 1ms)
//    vel(k)                        k 番目の tick の速度 [pulse/tick]
//    step(k)                       k 番目の tick で進む距離 [pulse]。k=1..T() の総和がちょうど L
//    s(k)                          k tick 時点の移動距離 (パラメータ s : 0〜L)
//    T()                           全時間 [tick]。必ず整数値
//
//  ・時間 t は「tick 番号」です。整数の k をそのまま渡せます。(割り込みごとに ++k するだけ)
//  ・全時間は整数 tick に切り上げ、その分だけ時間方向に引き伸ばします。
//    形は変わらず、速度・加速度は指定値より小さくなる方向にのみ変わります。(超えません)
//  ・メンバ関数 T() と衝突しないよう、テンプレート引数名は Num にしています。
//----------------------------------------------------------------------------------------------------------------------------------
//パナ方式
class trapezoid
{
public:	
	//各フェーズインデックスの定義です。
	enum{	ACC,	CONST,	DEC	,	N_PHASES};		//動作フェーズ
	enum{	ACC_UP,	ACC_DN,		};					//

	//=============================================================================================================
	//	各フェーズの動作シナリオを意味します。
	//=============================================================================================================
	struct	prof
	{
		int			n;					//	パルス出力回数
		//-----------------------------------------------------------------------------------
		//	移動距離です。これは計算で速度精度倍率を割った数値が入っています。
		//-----------------------------------------------------------------------------------
		long		l;												//	総移動距離(動作開始からの移動距離です。途中開始からも０スタートでよいのでは)
		long		l_amari;										//	総あまりパルス。これはどこかで吐き出さないとなりません。（これは距離の単位です。速度精度倍率は考慮されない）
		long		l_all(void) const {	return l + l_amari; 	}	//	上記 l + l_amariが吐き出すべきパルスです。
		//-----------------------------------------------------------------------------------
		//	速度
		//-----------------------------------------------------------------------------------
//		_tuple<accurized<>> v;				//	始点と終点の速度です。
		accurized<>		vc;
		accurized<>		vt;
		accurized<>		ac;					//	加速度（減速度符号付です)
		accurized<>		dv_amari;			//	あまりパルス分散のために、均等に付加する速度です。（一番初めの速度に加えます）
											//	※これは速度成分です。（速度精度倍率を含みます。
		void reset(){n=0;l=0;l_amari=0;}

		//※あまりパルスは基本的に計算出します。で
		prof():n(0),l(0),l_amari(0){}
		prof(int _n,long _l, const accurized<>&_vc , const accurized<>_vt , accurized<>_ac ):n(_n),l(_l),vc(_vc),vt(_vt), ac(_ac), dv_amari(), l_amari(0) {

		}

		prof & operator=(const prof &a){
			n = a.n;
			l = a.l;
			l_amari = a.l_amari;
			dv_amari = a.dv_amari;
			vc = a.vc;
			vt = a.vt;
			ac = a.ac;
			dv_amari = a.dv_amari;
			return *this;
		}
	};	//prof define
private:
protected:
	//これらは、
	const long				_L;						//	トータル距離 [nm または pulse?]
	const accurized<>		_vc;					//	初速度
	const accurized<>		_vt_order;				//	目標速度（三角動作の場合には、この通りにはならない場合があるので注意）
	accurized<>				_vt;					//	実際の到達速度
	const acc				_a;						//	加速
	//計算された動作プロファイルです。
	prof					_prof[N_PHASES];		//各フェーズの動作プロファイルです。

public:
	trapezoid() : _L(0) , _vc()  , _vt() , _a() {}			//すべて0です。
	trapezoid(long L , long vc , long vt , const acc& a ,int accuracy =1000) :
		_L(L) , _vc(accurized<>((long long)vc,accuracy)) , _vt_order(accurized<>((long long)vt,accuracy)),_a(a)
	{
		build();
	}
	//-----------------------------------------------------------------------------------
	//	あらかじめ望みの動作条件でオブジェクト構築した前提でそれを利用したメンバ関数としてみます。
	//-----------------------------------------------------------------------------------
	void build();		//これを呼べば、指定条件による動作プロファイルを作るとする。
	//以下は上記関数のサブルーチンでもはやよいのでは
//	bool	make_trianglar_profiles();		//	三角駆動のプロファイルを作る。


	//=================================================================================================
	//	一般的なサブルーチンとして定義します。
	//
	//	このクラスメンバの構築のため以外に一般的に台形計算でつかるようにstatic定義しておきます。
	//=================================================================================================
	static prof	_accel(const accurized<> &vc ,const accurized<> & vt  , const accurized<> &acc);	//	
	static prof	_constant(const accurized<>& v, long l);											//	等速フェーズ
	static prof	decel( const accurized<>& vc , const acc &a){	return	_accel( vc, accurized<>( 0LL , vc.accuracy() ) , a.a3d );	}
	static prof	accel(const accurized<> &vc ,const accurized<>&vt , const acc& a ){		return	_accel(vc , vt , (vt.raw()>vc.raw()) ? a.a1u : a.a1d );		}

	//
	static prof	dec_stop_at_target_distance(const accurized<>& _vc, const accurized<> &_acc, long target_distance);
	static accurized<>	correct_calc_miss(const accurized<> &vc, const accurized<>& vt, long l_target, const acc &ptp, prof &l1p,prof &l2, prof &l3);
	static int max_n_with_constant_accel(const accurized<>&vt, const accurized<>& a, long L);		//初速vc, int aで Lに一番近い距離まで進むnです。

	//三角駆動のための関数定義
	//これらは、台形計算のためのサブルーチンとして機能します（クラスメンバを使わない）
	static bool	triangle_is_possible(int n, const accurized<>& vc, const acc& a,long _L);
	static int	triangle_max_acc_n(const accurized<>&vc, const acc & a , long _L );
};


#endif	// __TRAPEZOID_H__
