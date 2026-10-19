#include "CommonModule.h"
#include "trapezoid.h"

#include <cmath>
#include <cstdint>



//-----------------------------------------------------------------------------------------------
//	三角駆動用の関数
//-----------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
//	tri_fits()	:	
//	初速 vc ～ aで規定された加速度で n 回 加速し減速するときに _Lの距離に対して
//	acc aで定義された
//	0	:	m 回の加速で行き過ぎる
//	1	:	m 回の加速でLより小さい
//------------------------------------------------------------------------------------------------
static inline __int64 acc_dist(accurized<>vc, accurized<> a1u, int m)
{
    return vc.raw() * m + a1u.raw() * (((__int64)m * (m + 1)) >> 1);
}
static inline __int64 dec_dist(accurized<>vp, accurized<>a3d, int n3)
{
    return vp.raw() * n3 - a3d.raw() * (((__int64)n3 * (n3 + 1)) >> 1);
}
bool trapezoid::triangle_is_possible(int n, const accurized<>& vc, const acc& a,long _L)
{
	_Assert(n>0 , "tri_fits:: n = 0");
	const __int64	L = (__int64)_L*vc.accuracy();								//速度精度倍率で計算します。

	const accurized<>new_vt(vc.raw() + a.a1u.raw() * n , vc.accuracy());		//
    const	int			n3	=	(new_vt.raw() - 1) / a.a3d.raw();            	// 減速回数 

    /* 倍精度で概算し、明らかに超過なら整数計算せず不成立にする。
     * （超過する候補で int64 が桁あふれするのを防ぐ。誤差は 10^-2 程度なので余裕 1e6 で十分） */
	//----------------------------------------------------------------------------------------------------------------------
	//	「行き過ぎ」判定
	//----------------------------------------------------------------------------------------------------------------------
	//上記で決めた回数 (m) で加速減速+(等速一回）+減速するときの距離(gd)を計算し、L+double誤差よりも大きくなってしまうとき、
	//　m 回加速してしまうと、距離Lに「収まらない」として 0 を返す。
	{
    	const double dn = (double)n , dn3 = (double)n3, d_new_vt = (double)new_vt.raw();
	    const double gd = (double)vc * dn + 0.5 * (double)a.a1u.raw() * dn * (dn + 1.0)			/*    加速域      */
                    + new_vt.raw() * dn3 - 0.5 * (double)a.a3d.raw() * dn3 * (dn3 + 1.0);		/*    減速域      */
	    if (gd > (double)L + 1.0e6) return false;
	}
	//----------------------------------------------------------------------------------------------------------------------
	//	ここまでくると「Ｌを超えない」です
	//----------------------------------------------------------------------------------------------------------------------
    /* ここまで来れば各項は 2^46 程度までなので int64 で厳密に計算できる */
    return acc_dist(vc, a.a1u, n ) + dec_dist(new_vt, a.a3d, n3) <= L;
}
int trapezoid::triangle_max_acc_n(const accurized<>&vc, const acc & a , long _L )
{
    /* 連続近似の解:
     *   m* = ( -vc(a1+a3) + sqrt( a3(a1+a3)(vc^2 + a1*vc + 2*a1*L) ) ) / ( a1(a1+a3) ) */
	const __int64 L = (__int64)_L*vc.accuracy();		//速度精度倍率で計算します。

	//等差数列の和をnについて二次関数で解きます。（解の公式）
    const double A = (double)a.a1u.raw(), D = (double)a.a3d.raw(), V = (double)vc.raw();
    const double s  = sqrt(D * (A + D) * (V * V + A * V + 2.0 * A * (double)L));
    const double nd = (s - V * (A + D)) / (A * (A + D));
	//
	_Assert( ((__int64)nd  / 0x7fffffff)==0 , "trapezoid::max_n_in_triangular n is too large");		//intで収まらない場合にエラーです。
    int n = (nd > 0.0) ? (int)nd : 0;    /* 32bit 変換で十分（m <= 約 3×10^6） */
	//
    while (triangle_is_possible(n + 1, vc , a , _L )) ++n;          /* 小さすぎたら増やす */
    while (n > 0 && !triangle_is_possible(n, vc, a , _L )) --n;    /* 大きすぎたら減らす */
    return n;
}
//これはpana製の計算ルーチンそのものです。
trapezoid::prof	trapezoid::_accel(const accurized<> &vc ,const accurized<> & vt  , const accurized<> &a)
{
	prof s;
	accurized<>acc(a*( vt<vc?-1LL:1) );												//絶対値をとるため
	if (vc == vt	)	{	s.reset();		return s;	}							//加速なし
	s.n		=	(int)( (abs((vt - vc).raw()) - 1) / abs(acc.raw()));					//加速減速回数です。
	s.l		=	accurized<>( (long long)((double)acc.raw() * (s.n + 1.0) * s.n	/ 2.0	+ vc.raw()*s.n) , vt.accuracy()).val();		//計算のあまりを考慮してみます。
	s.vc	=	vc;															//初速度
	s.vt.set( vc.raw() + (acc.raw() * s.n) , vc.accuracy()) ;						//目的速度
	s.ac			=	acc;														//加速度
	return s;
}
/*
trapezoid::prof trapezoid::accel(const accurized<> &vc , const accurized<>& vt , const acc &a)
{
	return	_accel(vc , vt , (vt.raw()>vc.raw()) ? a.a1u : a.a1d );
}
trapezoid::prof trapezoid::decel(const accurized<>& vc, const acc& a)				//減速フェーズ
{
	return	_accel( vc, accurized<>( 0 , vc.accuracy() ) , a.a3d );
}
*/
//等速運動
trapezoid::prof	trapezoid::_constant(const accurized<>& v, long l)					//等速フェーズ
{
	_Assert(v.raw() != 0, "axis::_constant() v is 0");
	prof s;
	s.n			= (int)((long long)l*v.accuracy() / v.raw());						//結局これは整数の計算でよいのでは・・
	s.l			= 	accurized<>( v.raw() * s.n , v.accuracy()).val();		//ここまでやる必要があるのか。
	s.vc = s.vt = v;
	s.l_amari	=	l - s.l;	//あまりパルスはここで発生します。
	return s;
}
//現在速度から減速するときの距離
trapezoid::prof	trapezoid::dec_stop_at_target_distance(const accurized<>& _vc, const accurized<> &_acc, long target_distance)
{
	_Assert(_vc.accuracy() == _acc.accuracy() , "_dec_stop_at_targetdistance accuracy is different(vc:%d , acc:%d)" , _vc.accuracy() , _acc.accuracy()); 
	//以下doubleで計算を進めていくので、
	//vc,accはdoubleで扱うため、作り直します。
	accurized<double>vc((double)_vc.raw(), _vc.accuracy());
	accurized<double>acc((double)_acc.raw(), _acc.accuracy());
	prof s;
	double dval		=	(double)( (vc *2.0)- acc);
	dval			=	(dval * dval)	-	((acc * 8.0)* accurized<double>( target_distance*vc.accuracy() , vc.accuracy() ) ).raw();
	s.n				=	(long)(( (2.0 * vc.raw()) - acc.raw() - sqrt( dval)) / (2.0 * acc.raw()));
	s.l				=	accurized<>(	(long long)( ((double)vc * s.n) - (((double)acc * ((double)s.n + 1.0) * (double)s.n) / 2.0) ) , vc.accuracy()).val();
	s.vt			=	accurized<>( 0LL , _vc.accuracy());			//
	s.ac			=	_acc;
	return	s;
}

//==========================================================================================================================
//	パナのコードで以下の処理があり、これまでの経路計算で精度の部分でミスがあるためその補正。
//		/*--- memo ---------------------------------------------------------------------------------------------*/
//		/* 以下の条件が全て成立した場合に、事前計算による　L1U,L2U,L3U　の結果が　すべて0となってしまう。		*/
//		/* 　条件①　余りパルス分散設定																			*/
//		/* 	 条件②　速度精度倍率が１以外																		*/
//		/*	 条件③　加速度が送り速度よりも大きい場合（加減速域、減速域がない場合）								*/
//		/*   条件④　目標速度vtの値が小数部の値を持つ場合														*/
//		/*   条件⑤　移動量が送り速度の整数部と一致している場合													*/
//		/*   【原因】目標速度vtを速度精度倍率で割り算した際に、余りを切り捨てた後、目標速度に到達できる移動量が	*/
//		/*　　　　　 残移動量に確保できるか判断しているため、本来確保できないケースにも関わらず					*/
//		/*           確保できると誤った判断となり、目標速度を補正して移動量を再計算する処理が実行されなかった。	*/
//		/*			 そのため、この処理の前段処理で計算した結果が最終結果となり、								*/
//		/*           加減速、減速のパルス出力回数=0　また、等速域のパルス出力回数=0となり、事前計算にて移動量が */
//		/*			 算出されない。																				*/
//		/*------------------------------------------------------------------------------------------------------*/
//	tokuとりあえず、パナコードをそのまま流用します。
//==========================================================================================================================
//accurized<>	axis::_correct_calc_miss(long long vc , long long vt , int l_target ,acc ptp, pulse_scenario_CP &l1 , pulse_scenario_CP &l2 , pulse_scenario_CP &l3)
accurized<>	trapezoid::correct_calc_miss(const accurized<> &vc, const accurized<>& vt, long l_target, const acc &ptp, prof &l1,prof &l2, prof &l3)
{
	accurized<>new_vt(vt);
	long long	lltemp  = (long long)l1.l + vt.val() + l3.l;						//等速距離が１クロックだけあった場合の距離
	if(		(	lltemp > l_target	)		|| 										//↑が、目標距離より大きい。
			//	^^^^^^^^^^^^^^^^^	いきすぎる場合
			(	(lltemp == (long long)l_target) && ( vt.dcm() != 0 )	)
			//		^^^^^^上記距離がピッタリ目標位置	＆	(vt%_v_accuracy	!= 0)	精度倍率で足がでている状況
	){
		/*** 現指令速度(vc) ≧ 目標速度(vt)の場合 ***/
		if( vc >= vt)	{		//	目標速度に現指令速度を設定して開始時点から減速域にする */
								//	減速域のみの出力パルスでは残移動パルスに足りないので、
								//	減速中滑らかになるタイミングで余りパルス出力で調整
			new_vt = vc;	/*	三角駆動防止補正速度=vcとする */
			l1.reset();			//L1（加速部分はない)
			l2.reset();
			l3		=	decel(vc,ptp);			//toku_?	パナコードの場合は、この計算ですが
											//toku_?	本当に l3.l == l_targetとなるか？が心配。
											//			バグの場合、l_targetと同じ移動量になるように計算する。
		} else {	/*** 現指令速度(vc) < 目標速度(vt)の場合 ***/
			/* 目標速度を補正する */
			/* 補正した目標速度で「加減速域＋等速域１回（補正後速度）＋減速域＝残移動量」の式を組み立て新目標速度を求める */		
#if 1
{
	static double __debug = (double)accurized<>((long long)l_target * new_vt.accuracy(), new_vt.accuracy() );
	__debug=0;
}	
#endif		
			double dval  = (	((double)ptp.a3d * (double)vc*(double)vc) 
							+	((double)ptp.a3d * (double)ptp.a1u * (double)vc)
							+	(2.0 * (double)ptp.a3d * (double)ptp.a1u * (double)accurized<>( (long long)l_target * new_vt.accuracy(), new_vt.accuracy() ) ) 
							)
							*	((double)ptp.a3d + (double)ptp.a1u);
			dval		+=	(((double)ptp.a1u * (double)ptp.a3d * (double)ptp.a3d * (double)ptp.a1u)
			         - (2.0 *(double)ptp.a3d * (double)ptp.a3d * (double)ptp.a1u)  + ((double)ptp.a1u * (double)ptp.a1u));
			/*** Doubleの計算結果の小数部を四捨五入してlongデータvtに反映させる ***/
			new_vt.set( (long long)((((double)sqrt(dval) - (double)ptp.a3d*((double)ptp.a1u -1.0))/((double)ptp.a1u + (double)ptp.a3d)) +0.5 ) , new_vt.accuracy() );
			//vc -> new_vt で再度計算します。
			l1		=	accel( vc, new_vt , ptp );
			l3		=	decel( new_vt ,ptp);
			l2		=	_constant( new_vt, l_target - (l1.l+l3.l) );
		}
	}
	return new_vt;
}
//================================================================================================================================
//		台形と三角駆動を作ります。
//		(可能な _vt と、 メンバー_vc,_vt_aを使ったプロファイル_prof[]を作ります。)
//================================================================================================================================
void trapezoid::build()		//これを呼べば、指定条件による動作プロファイルを作るとする。
{
	//減速のみで目標位置まで
	if(decel(_vc,_a).l > _L){		//減速距離が動作距離よりも長いケース
		_prof[DEC]	=	decel(_vc,_a);
		_prof[ACC].reset();
		_prof[CONST].reset();
	}else{
		//総距離よりも大きな場合は、三角駆動として新たにプロファイルを作ります。
		if( accel(_vc,_vt_order,_a).l_all() + decel(_vt_order,_a).l_all() > _L) {	//三角駆動させて_Lを超えてしまう
			//新たなvtを作ります。
			int n_acc = triangle_max_acc_n(_vc,_a,_L);		//_Lに収まる最大の加速回数nを求めます。
			_vt = accurized<>(_vc.raw() + _a.a1u.raw() * n_acc, _vc.accuracy());		//新たな目標速度です。
		}else{										//そうでない場合は台形がつくれます。
			_vt = _vt_order;						//目標速度にできます。
		}
		//目標速度が決まったので動作プロファイルを作ります。
		_prof[ACC] = accel(_vc,_vt,_a);
		_prof[DEC] = decel(_vt,_a);
		_prof[CONST] = _constant( _vt , _L - (_prof[ACC].l_all() + _prof[DEC].l_all()) );		
	}
}
