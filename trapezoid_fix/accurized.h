#ifndef	__ACCURIZED_H__
#define	__ACCURIZED_H__

//==============================================================================================================================
//	精度倍率を考慮した値を表現します。
//==============================================================================================================================
//64ビットの値と、それを累積していくしくみを一つ作ってみます。
template<class T = long long, class DST = long>			//一応テンプレートにしておきます
class accurized
{
private:
	T			__raw;					//	×  accuracy倍された生の値。
	//toku 以下、val,dcmは__rawから計算で出せるので、なくせるか検討する。問題になるのは速度だろう。
	DST			__val;					//	最後に精度倍率で割り算した値です。
	DST			__dcm;					//	累積の際のあまり。次回の累積のときに足されます。
	DST			__accuracy;				//	このモジュールで最後に計算に使われた精度倍率です。（このモジュールでは精度倍率は持ちません）
public:
	void		reset(bool reset_accuracy = false) {
		__raw = 0;	__val = __dcm = 0;
		if (reset_accuracy) {
			__accuracy = 1;
		}
    }
	//中身がないかどうかを判定します。
	bool		empty(void)const {
        if ((__raw == 0) && (__val == 0) && (__dcm == 0)) { return true; }
	    return false;
	}

	accurized() {
		reset(true);
	}	//default
	accurized(const accurized &a) {
		*this = a;
	}	//copy
/*	accurized(DST accuracy)						{
		if (accuracy == 0) {
			__val = 0;
			return;
		}
		set(0,accuracy);					}*/
	accurized(T   raw, DST accuracy) {
		set(raw, accuracy);
	}
		//toku これは危険なので廃止です。
	accurized(DST val, DST accuracy) {
		_Assert(0,"accurized<>accurized(DST,DST) is not avail");
		set((T)val*accuracy, accuracy);
	}		//ちょっと紛らわしいか

		//setter
	DST		set(T raw, DST accuracy) {
		_Assert(accuracy > 0, "accurized::set() accuracy can't set  0");	//ここのチェックは必要ないけど一応
	//		reset();		//これまでの累積はなしにしてしまう
	//		return accumlate(raw,accuracy);
		//この操作は、持ち越されたあまり(dcm)は0、rawとval、accuracyをセットするようにします。
		__raw = raw;		//
		__accuracy = accuracy;	//
		__val = raw / __accuracy;	//
		__dcm = 0;					//
		return	__val;
	}
	//必要性が出てきたのでvalでセットするものを用意します。
	virtual DST		set_by_val(DST val) {
		return set((T)val*__accuracy, __accuracy);
	}

	//※ここの計算速度が問題なら、1024にすればビットシフトで大丈夫になるのでそうします。
	DST		accumlate(T raw, DST accuracy) {		//値を累積していきます。(前回までのあまりを考慮)
		if (__accuracy != accuracy) {			//指示された値の精度倍率が異なる場合、自分自身はカラの状態でないとなりません。
			_Assert(empty(), "accurate_val::accumlate() : object is not empty but accuracy is deffernt(%d -> %d) .", __accuracy, accuracy);
			__accuracy = accuracy;							//まず精度倍率を確定します。
		}
		//rawの値は今回のものに更新です。
		__raw = raw;									//これまでの累積あまりを考慮します。
		//精度倍率の低い値は、前回の
		{
			long long v = __raw + __dcm;		//前回使った古いdcm
	//		__val			=	(DST)( (__raw + (T)__dcm) / (T)__accuracy);		//									※精度倍率をキリのいい値にすればビットシフトでOK
			__val = (DST)(v / (T)__accuracy);						//									※精度倍率をキリのいい値にすればビットシフトでOK
			//※ここで__vがDSTの値の範囲内におさまっているかが心配。
			//計算速度のため
			//		dcm		=	(long)(v % (T)accuracy);				//新たな累積あまりです。
			//		__dcm = (__raw+(T)__dcm)%__accuracy;
			//		__dcm = (v	%	__accuracy;
			__dcm = (DST)(v - ((T)__val*__accuracy));								//精度倍率で割り算したあまりです。	※精度倍率をキリのいい値にすればビットシフトでOK
		}
		//ここで割り算後はlongの範囲内に収まっているはずですのでチェックはしておきます。
		return	__val;
    }
	DST		accumlate(const accurized &v) { return accumlate(v.raw(), v.accuracy()); }

	//	set_accuracy(long _acc)	{	_Assert(acc>0,"val_accuracy : acc is zero");	accuracy=_acc;	}
	bool	operator==(const accurized &a)const {
		//テンプレートの型が同じかどうか
		return		((__raw == a.raw())
			&& (__dcm == a.dcm())
			&& (__val == a.val())
            && (__accuracy == a.accuracy()));
			//toku ここは全部考慮する
	}
	accurized<> operator-(const accurized<> &a) const {
		_Assert(accuracy() == a.accuracy(), "accurized::operator- accuracy different(%d,%d)", accuracy(), a.accuracy());
		accurized<T, DST> b(raw() - a.raw(), accuracy());
		//		_Assert(a == b, "accurized::operator-() calc wrong");		//一応検算です。
		return b;		//これでいいのか
	}
	accurized<> operator+(const accurized<>&a)const {
#if 0
			if (__raw < 0) {
				__cm_debug = 1;
			}

#endif
		_Assert(accuracy() == a.accuracy(), "accurrized<> operator+:: different accuracy ");
		accurized<>b(__raw + a.raw(), accuracy());
		//toku ※これまで累積した値（__dcm）への考慮がなくなる。それがいいかどうか？
		return b;
    }
	accurized<T, DST> & operator=(const accurized<T, DST> &a) {
		__raw = a.raw();					//	×  accuracy倍された生の値。
		__val = a.val();					//	最後に精度倍率で割り算した値です。
		__dcm = a.dcm();					//	累積の際のあまり。次回の累積のときに足されます。
		__accuracy = a.accuracy();				//	このモジュールで最後に計算に使われた精度倍率です。（このモジュールでは精度倍率は持ちません）
		return *this;
	}
	//operators
	accurized<T, DST>		operator*(T a)const { accurized<T, DST>b(__raw * a, accuracy());return b;	/*toku これまで累積した__dcmが消える.問題ないかどうか判別できなかった*/ }
	accurized<T, DST> &		operator*=(T a) { *this = (*this)*a;	return *this; }
	bool					operator<(const accurized<>&a)const { return __raw < a.raw(); }
	bool					operator<=(const accurized<>&a)const { return __raw <= a.raw(); }
	bool					operator>=(const accurized<>&a)const { return __raw >= a.raw(); }
	bool					operator>(const accurized<>&a)const { return __raw > a.raw(); }
	bool					operator>(long long a)const { return __raw > a; }
	accurized<>&			operator+=(const accurized<>&a) { *this = operator+(a);	return *this; }
	//キャスト演算子を定義します。なるべく意味不明なのはやめる
		//					operator  T() const		{	return raw();}		//
	//	template<class NEW_T>		operator NEW_T() const {	return(NEW_T)__raw;}		//toku これはコンパイルが通らない。理由はわからなかった。
	operator double() const { return (double)__raw; }
    //	operator T() const			{	return	raw();	}
	template<class NEW_T>	operator accurized<NEW_T, T>() const { return accurized<NEW_T, T>((NEW_T)raw()), accuracy(); }

	//精度倍率の変更を作ってみます。
	accurized<T, DST> & chg_accuracy(DST new_accuracy) {
		if (accuracy() == new_accuracy) { return *this; }	//変化なしの場合には何もしない（オーバーヘッド軽減のため）
		//現状の__rawは、古い精度倍率でかけ算されたものです。
		_Assert(new_accuracy > 0, "accurized::chg_accuracy() is 0");
        //今のところ、累積がある場合には精度倍率変更ができないようにします。
		//何か手法があればいいが
		_Assert(__dcm == 0, "accurized<>::chg_accuracy() : not allowed when __dcm exist(%d)", __dcm);

		//かけ算のオーバーフローを見ます。
		{
			//			「a < LONG_MAX / b が成り立つ場合、a * b はオーバーフローしない」
			T max = std::numeric_limits<T>::max();
			_Assert(__raw < (max / new_accuracy), "__raw or accuracy is too big to covert accuracy");
		}
		//リセットして値を入れ直します。
		set(__raw * new_accuracy / __accuracy, new_accuracy);		//これが成り立てば、精度倍率の変更がいつあっても大丈夫です。
																		//ですが、現在は精度倍率変更は動作中は行わない制約を設けます。
		//※この処理は、いままでの累積がなくなってしまいます。			//
		//累積の値も新しい制度倍率にいれるか？
		return *this;
	}

	//絶対値をとるようにしてみます。
	accurized<> abs(void) const {
		if (__raw < 0) {
			return accurized<>(*this)*-1LL;
		}
		return *this;
	}

	//getter:
	T		raw(void)			const { return		__raw; }
	DST		val(void)			const { return		__val; }
	DST		dcm(void)			const { return		__dcm; }
	DST		accuracy(void)		const { return		__accuracy; }

	//あまり部分をクリアする関数を用意します。multi_pulse_generator(axis)が、累積あまりがある状態で、
	// generatorが次に進むと、あまり0スタートになるので、それに合わせないと計算があまり分ずれてしまうため。
	//
	void	clr_dcm()					{	__dcm=0;	}
};
#endif//
