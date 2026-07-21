# 再開メモ (2026-07-21時点)

別PCでこの続きをやるための最小限のセットアップ手順。プロジェクトの背景・設計・
現状の数値は `CLAUDE.md` に詳しくまとめてあるので、作業前に必ず目を通すこと。

## セットアップ

```
git clone https://github.com/naoyatoku/cut.git
cd cut
pip install -r requirements.txt
python zhop/zhop_planner.py            # PNG生成のみ
python zhop/zhop_planner.py --show     # + hop_plan_3d.html をブラウザで自動オープン
```

`--show` は matplotlibのGUIバックエンド(Tk/Qt)には依存しない(過去に環境依存で
動かない事例があり、Plotly + ブラウザ表示に切り替え済み)。`plotly` が無い場合は
警告を出してPNG出力のみ継続する。

## 出力ファイル(すべて `zhop/` 配下、gitには含めていない生成物もあり)

- `hop_commands.csv` — RTEX向け1ms周期位置指令(git管理下、デモ値で再現可能)
- `hop_plan.png` / `hop_plan_3d.png` — 静止画可視化(git管理下)
- `hop_plan_3d.html` — インタラクティブ3D(**gitignore対象**、`--show`実行時に毎回生成される)

## 次にやること(優先度は特になし、CLAUDE.mdの「次にやること」と同期)

1. Z_safeオーバーシュートによる先行下降(下降テールの短縮)
2. S字(ジャーク制限)プロファイル化(`make_trapezoid` を差し替えるだけの構造)
3. 1ウェハ分シーケンスのライン毎時間一覧

## 直近のセッションで詰めた内容(要点だけ)

- 加速時定数は各軸100ms(ユーザー指定)で確定 → `T_ACC=0.1`
- Zホップは「カット終端のオーバートラベル走行中にZ上昇を先行開始」「次ライン助走中に
  Z下降を継続」の両方を実装済み。**ただし禁止円内(=加工中)はX等速・YZ完全固定**という
  制約があり、Z上昇/Xの減速は刃が禁止円を退出した後のマージン区間(オーバートラベル10mm
  − XYマージン3mm = 7mm)でのみ許される。この制約を`compute_cut_approach`が
  `CUT_FEED_SPEED`(仮値50mm/s)で可視化・判定している(`z_head_ok`, `decel_ok`)
- カット送り速度(`CUT_FEED_SPEED`)とパルス変換係数(`PULSES_PER_MM`)はまだ仮値。
  実機値が分かり次第差し替えが必要
- **z_head不足時の実機挙動(xy_delay)を実装済み**: z_headがt_availを超える場合(送りが
  速いなど)、エラーで止めるのではなく「超過分だけAに到着後、Z上昇待ちでXY出発を遅らせる」
  という実機の現実的な挙動をモデル化した。`plan_hop(..., t_avail=...)`に渡すと
  `HopPlan.xy_delay = max(0, z_head - t_avail)`が自動計算される。現状のデモ値
  (送り50mm/s)ではz_head_ok=OKなのでxy_delay=0(従来と同じ)。送り100mm/sの
  テストでxy_delay=23.2msが正しく発動し、連続性・安全チェックとも正常なことを確認済み

## 未解決の確認事項

- B側(次ラインの助走)についても、Z下降テール(98ms)が助走時間に収まるかの明示チェックは
  未実装(A側の`z_head_ok`/`xy_delay`に相当するものがB側にはまだ無い)。対称的に追加できる
- カット送り速度・パルス変換係数の実機値
