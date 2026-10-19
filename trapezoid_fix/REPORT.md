# trapezoid 検証レポート

対象: `trapezoid.h` / `trapezoid.cpp` / `accurized.h`(Panasonic方式の台形・三角駆動)
方法: 欠けていた `CommonModule.h` / `tuple.h` をスタブ化し g++ 13 でビルド。
ASan/UBSan で実行し、L・vc・vt・精度倍率・加速時間を振って各フェーズの距離・余り・到達速度を検査。

## 結論
元のコードはそのままでは動作しない(コンストラクタが必ず assert)。
下記3点を修正すれば、台形・三角駆動の計算は正しく動作する。

## 修正したバグ
| # | 箇所 | 内容 | 修正 |
|---|------|------|------|
| 1 | `trapezoid::trapezoid(long,long,long,...)` | `accurized<>(vc, accuracy)` が `(DST,DST)` オーバーロードに解決され、常に `_Assert(0)` で停止 | `(long long)` キャスト |
| 2 | `prof()` / `prof::reset()` | `l_amari` ほか未初期化。空フェーズの `l_all()` にゴミが混入し、CONST長が負になる(例 `n=-1,l=-3`) | コンストラクタで 0 初期化、`reset()` でも `l_amari=0` |
| 3 | `acc::operator=` | `return` 無し(未定義動作)、`a1d = a.a3d`(`a.a1d` の誤り) | 修正 |
| 4 | `acc_dist` / `dec_dist` | 戻り値 `int`、`m*(m+1)` も32bitで桁あふれの恐れ | `__int64` 化 |
| 5 | `#endif __TRAPEZOID_H__` | 警告 | コメント化 |

## 修正後の確認結果
- 台形 L=300000: 加速99 + 等速201 + 減速99、合計 = L
- 三角駆動 L=50000: 到達速度 700 に制限、合計 = L
- vc>0、vt<vc も正常
- 余り (`CONST.l_amari`) は常に 0 以上 1tick 分以内、到達速度は指定を超えない
- 総当たり 414 パターンで不一致 0(下記2点を除く)
- 未初期化修正前は 1008 パターン中 210 件で余りが負

## 未対応(仕様判断が必要)
1. **短距離で例外**: `vc=0` で L が 1回分の加速距離より短いと三角駆動の `n=0` → `_vt=0` → `_constant()` が assert。
   例: 精度倍率1・加速度10/tick² で L<10。最低速度で L パルス出す等のフォールバックが必要。
2. **減速距離 > L**: `vc>0` で止まりきれないとき DEC のみとなり合計が L を超える
   (例 vc=500, L=1000 → 12250 進む)。呼び出し側での扱いを決める必要あり。

## 備考
- ctor の `vc`/`vt` は精度倍率込みの raw 値(pulse/tick × accuracy)。単位に注意。
- `correct_calc_miss` は `build()` から呼ばれておらず未検証。`max_n_with_constant_accel` は宣言のみで未定義。
- `accurized::operator+/-` は `__dcm` を落とす(コード内コメントの通り)。
- ACC終了速度(`vc+ac*n`)と CONST速度 `_vt` に最大1加速ステップ分の段差がある(Panasonic方式の仕様と思われる)。

## 再現
`cd test && make`
