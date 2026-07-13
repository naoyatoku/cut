# -*- coding: utf-8 -*-
"""
ダイシング装置向け Zホップ軌道プランナ(プロトタイプ)

1ラインカット終了点 A から次ライン開始点 B への移動を、
ウェハ上空 z_safe 以下への進入禁止制約のもとで時間最短化する。

方式:
  - XY は A→B 直線を各軸 vmax/amax から導いた経路プロファイル(台形)で移動
  - Z は「上昇 → z_safe 待機 → 下降」の台形プロファイル
  - 衝突制約「XYが禁止円内にいる間は z >= z_safe」を満たす範囲で
    3軸の動作を最大限オーバーラップさせる(閉形式)
  - 出力: 1ms周期(RTEX)の各軸位置指令列 CSV + 可視化 PNG

パルス変換は PULSES_PER_MM の係数を掛けるだけ(現状は仮値)。
"""

import math
import csv
import os
import sys
from dataclasses import dataclass

import numpy as np
import matplotlib
if "--show" not in sys.argv:   # --show 付きで実行するとウィンドウ表示(3Dを回転できる)
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================
# パラメータ(装置に合わせてここを書き換える)
# ============================================================

# 軸性能 [mm/s]
VMAX = {"x": 1000.0, "y": 200.0, "z": 50.0}
# 各軸とも 100ms で最高速到達(ユーザー指定の加速時定数)
T_ACC = 0.100
AMAX = {ax: v / T_ACC for ax, v in VMAX.items()}  # X=10000, Y=2000, Z=500 mm/s^2

# パルス変換係数 [pulse/mm] ※仮値。実機の電子ギア比に合わせる
PULSES_PER_MM = {"x": 1000.0, "y": 1000.0, "z": 1000.0}

RTEX_DT = 0.001  # 指令周期 1ms

# ワーク(ウェハ)ジオメトリ
WAFER_CENTER = (0.0, 0.0)   # XY中心 [mm]
WAFER_RADIUS = 150.0        # 300mmウェハ [mm]
XY_MARGIN = 3.0             # 刃・フランジ等のXY方向余裕(保守的に一様膨張) [mm]
Z_SAFE = 1.5                # ウェハ上面(z=0)からの退避高さ [mm]

# シナリオ: 片方向カットの戻り(右端で終了 → 左端の次ライン開始へ)
CUT_DEPTH = -0.8            # カット中の刃最下点 [mm](ウェハ上面基準)
OVERTRAVEL = 10.0           # ウェハ縁からのオーバートラベル [mm]
INDEX_PITCH = 5.0           # ラインピッチ [mm]

# ============================================================
# 台形速度プロファイル
# ============================================================

@dataclass
class Trapezoid:
    """距離 dist を vmax/amax 制約下で最短時間で移動する台形(または三角)プロファイル"""
    dist: float
    v: float      # 実到達速度
    a: float
    ta: float     # 加速時間
    tc: float     # 定速時間

    @property
    def T(self) -> float:
        return 2.0 * self.ta + self.tc

    def s(self, t: float) -> float:
        """時刻 t での移動距離"""
        if t <= 0.0:
            return 0.0
        if t >= self.T:
            return self.dist
        if t < self.ta:
            return 0.5 * self.a * t * t
        s_acc = 0.5 * self.a * self.ta * self.ta
        if t < self.ta + self.tc:
            return s_acc + self.v * (t - self.ta)
        td = self.T - t
        return self.dist - 0.5 * self.a * td * td

    def vel(self, t: float) -> float:
        if t <= 0.0 or t >= self.T:
            return 0.0
        if t < self.ta:
            return self.a * t
        if t < self.ta + self.tc:
            return self.v
        return self.a * (self.T - t)

    def time_at_s(self, s: float) -> float:
        """移動距離 s に達する時刻(逆関数)"""
        if s <= 0.0:
            return 0.0
        if s >= self.dist:
            return self.T
        s_acc = 0.5 * self.a * self.ta * self.ta
        if s < s_acc:
            return math.sqrt(2.0 * s / self.a)
        if s < s_acc + self.v * self.tc:
            return self.ta + (s - s_acc) / self.v
        return self.T - math.sqrt(2.0 * (self.dist - s) / self.a)


def make_trapezoid(dist: float, vmax: float, amax: float) -> Trapezoid:
    dist = abs(dist)
    if dist < 1e-12:
        return Trapezoid(0.0, 0.0, amax, 0.0, 0.0)
    v_tri = math.sqrt(dist * amax)
    if v_tri <= vmax:  # 三角プロファイル
        ta = v_tri / amax
        return Trapezoid(dist, v_tri, amax, ta, 0.0)
    ta = vmax / amax
    tc = (dist - vmax * ta) / vmax
    return Trapezoid(dist, vmax, amax, ta, tc)


# ============================================================
# ホップ計画
# ============================================================

@dataclass
class HopPlan:
    A: tuple
    B: tuple
    xy_prof: Trapezoid
    xy_dir: tuple          # 単位ベクトル
    up_prof: Trapezoid     # zA -> z_safe
    down_prof: Trapezoid   # z_safe -> zB
    z_head: float          # Z上昇の先行開始時間(A到達の何s前に開始するか)
    t_down_start: float    # Z下降開始時刻(t=0 は XY開始 = A出発)
    t_enter: float         # 禁止円進入時刻(交差なしは -1)
    t_exit: float
    total: float           # A出発からZ下降完了までの時間
    total_sequential: float

    @property
    def down_tail(self) -> float:
        """XY完了後にZ下降が残る時間(次ラインの助走中に下降を継続できる)"""
        return max(0.0, self.t_down_start + self.down_prof.T - self.xy_prof.T)


def plan_hop(A, B, center, r_forbid, z_safe) -> HopPlan:
    ax_, ay, az = A
    bx, by, bz = B

    # --- XY 直線経路のプロファイル ---
    dx, dy = bx - ax_, by - ay
    L = math.hypot(dx, dy)
    if L < 1e-12:
        ux, uy = 1.0, 0.0
        xy_prof = make_trapezoid(0.0, VMAX["x"], AMAX["x"])
    else:
        ux, uy = dx / L, dy / L
        # 経路方向の実効 vmax/amax(遅い軸に律速される)
        v_path = min(VMAX["x"] / max(abs(ux), 1e-12),
                     VMAX["y"] / max(abs(uy), 1e-12))
        a_path = min(AMAX["x"] / max(abs(ux), 1e-12),
                     AMAX["y"] / max(abs(uy), 1e-12))
        xy_prof = make_trapezoid(L, v_path, a_path)

    # --- Z プロファイル ---
    up_prof = make_trapezoid(z_safe - az, VMAX["z"], AMAX["z"])
    down_prof = make_trapezoid(z_safe - bz, VMAX["z"], AMAX["z"])

    # --- 禁止円との交差区間(経路パラメータ s) ---
    cx, cy = center
    ox, oy = ax_ - cx, ay - cy
    b_half = ox * ux + oy * uy
    c = ox * ox + oy * oy - r_forbid * r_forbid
    disc = b_half * b_half - c
    s_enter = s_exit = None
    if disc > 0.0 and L > 1e-12:
        s1 = -b_half - math.sqrt(disc)
        s2 = -b_half + math.sqrt(disc)
        if s2 > 0.0 and s1 < L:
            s_enter = max(s1, 0.0)
            s_exit = min(s2, L)

    seq = up_prof.T + xy_prof.T + down_prof.T

    if s_enter is None:
        # 交差なし: Z直行 + XY を完全オーバーラップ
        z_direct = make_trapezoid(bz - az, VMAX["z"], AMAX["z"])
        total = max(xy_prof.T, z_direct.T)
        return HopPlan(A, B, xy_prof, (ux, uy), z_direct,
                       make_trapezoid(0.0, VMAX["z"], AMAX["z"]),
                       0.0, total, -1.0, -1.0, total, seq)

    t_enter_rel = xy_prof.time_at_s(s_enter)
    t_exit_rel = xy_prof.time_at_s(s_exit)

    # Z上昇の先行開始(ユーザー承認済み):
    # カット終端のオーバートラベル走行中(刃はウェハ外)に Z上昇を始めてよいので
    # XY開始遅延は常にゼロ。z_head = A到達の何s前に上昇開始が必要か。
    # ※カット終端の x=+r_forbid→A の走行時間が z_head 以上あることは要確認
    #   (カット送り速度に依存)。不足する場合はその分だけXY開始が遅れる。
    z_head = max(0.0, up_prof.T - t_enter_rel)

    # Z下降は円退出と同時に開始。XY完了後に残る下降テール(down_tail)は
    # 次ラインの助走(オーバートラベル)中に継続してよい(ユーザー承認済み)
    t_down_start = t_exit_rel
    total = max(xy_prof.T, t_down_start + down_prof.T)

    return HopPlan(A, B, xy_prof, (ux, uy), up_prof, down_prof,
                   z_head, t_down_start, t_enter_rel, t_exit_rel,
                   total, seq)


def sample_plan(plan: HopPlan, dt: float = RTEX_DT):
    """1ms周期で各軸位置をサンプリング

    t=0 が XY開始(A出発)。z_head > 0 の場合は t<0 の行(Z先行上昇、XY増分ゼロ)
    を含む。この区間の dz はカット終端の指令列に合成して使う。
    """
    n_head = int(math.ceil(plan.z_head / dt))
    n_main = int(math.ceil(plan.total / dt)) + 1
    t = np.arange(-n_head, n_main + 1) * dt
    ax_, ay, az = plan.A
    bx, by, bz = plan.B
    ux, uy = plan.xy_dir

    x = np.empty_like(t)
    y = np.empty_like(t)
    z = np.empty_like(t)
    for i, ti in enumerate(t):
        s = plan.xy_prof.s(ti)
        x[i] = ax_ + ux * s
        y[i] = ay + uy * s
        if plan.t_enter < 0:
            # 交差なし: up_prof に zA->zB 直行が入っている
            sz = plan.up_prof.s(ti)
            z[i] = az + math.copysign(sz, bz - az) if plan.up_prof.dist > 0 else az
        else:
            if ti < plan.t_down_start:
                z[i] = az + plan.up_prof.s(ti + plan.z_head)
            else:
                z[i] = (az + plan.up_prof.dist) - plan.down_prof.s(ti - plan.t_down_start)
    # 終端を厳密に一致させる
    x[-1], y[-1], z[-1] = bx, by, bz
    return t, x, y, z


def write_csv(path, t, x, y, z):
    """RTEX向け: 位置[mm]・1ms増分[mm]・パルス位置を出力"""
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["t_ms",
                    "x_mm", "y_mm", "z_mm",
                    "dx_mm", "dy_mm", "dz_mm",
                    "x_pulse", "y_pulse", "z_pulse"])
        for i in range(len(t)):
            dx = x[i] - x[i - 1] if i > 0 else 0.0
            dy = y[i] - y[i - 1] if i > 0 else 0.0
            dz = z[i] - z[i - 1] if i > 0 else 0.0
            w.writerow([f"{t[i]*1000:.0f}",
                        f"{x[i]:.6f}", f"{y[i]:.6f}", f"{z[i]:.6f}",
                        f"{dx:.6f}", f"{dy:.6f}", f"{dz:.6f}",
                        round(x[i] * PULSES_PER_MM["x"]),
                        round(y[i] * PULSES_PER_MM["y"]),
                        round(z[i] * PULSES_PER_MM["z"])])


def plot_plan(plan: HopPlan, t, x, y, z, r_forbid, out_png):
    for name in ("MS Gothic", "Yu Gothic", "Meiryo"):
        try:
            matplotlib.rcParams["font.family"] = name
            break
        except Exception:
            pass
    matplotlib.rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(2, 2, figsize=(13, 8))

    # 側面図 (X-Z)
    a = axes[0][0]
    cx = WAFER_CENTER[0]
    a.add_patch(plt.Rectangle((cx - WAFER_RADIUS, -0.9), 2 * WAFER_RADIUS, 0.9,
                              color="0.75", label="ウェハ"))
    a.axhline(Z_SAFE, ls="--", color="gray", lw=1)
    a.plot(x, z, "-", color="tab:green", lw=2)
    a.plot([plan.A[0], plan.B[0]], [plan.A[2], plan.B[2]], "o", color="tab:green")
    a.set_xlabel("X [mm]"); a.set_ylabel("Z [mm]")
    a.set_title("側面図 X-Z")
    a.grid(alpha=0.3)

    # 上面図 (X-Y)
    a = axes[0][1]
    th = np.linspace(0, 2 * np.pi, 200)
    a.plot(WAFER_CENTER[0] + WAFER_RADIUS * np.cos(th),
           WAFER_CENTER[1] + WAFER_RADIUS * np.sin(th), color="0.6")
    a.plot(WAFER_CENTER[0] + r_forbid * np.cos(th),
           WAFER_CENTER[1] + r_forbid * np.sin(th), "--", color="0.6", lw=1)
    a.plot(x, y, "-", color="tab:green", lw=2)
    a.set_xlabel("X [mm]"); a.set_ylabel("Y [mm]")
    a.set_title("上面図 X-Y(破線=膨張フットプリント)")
    a.set_aspect("equal"); a.grid(alpha=0.3)

    # 位置時系列
    a = axes[1][0]
    a.plot(t * 1000, x, label="X")
    a.plot(t * 1000, y, label="Y")
    a.plot(t * 1000, z * 50, label="Z ×50")
    if plan.t_enter >= 0:
        a.axvspan(plan.t_enter * 1000, plan.t_exit * 1000, color="orange", alpha=0.15,
                  label="円内通過区間")
    if plan.z_head > 0:
        a.axvspan(-plan.z_head * 1000, 0, color="tab:blue", alpha=0.10,
                  label="Z先行上昇(カット終端走行中)")
    a.set_xlabel("t [ms]"); a.set_ylabel("位置 [mm]")
    a.set_title("位置指令(1ms周期)")
    a.legend(); a.grid(alpha=0.3)

    # 速度時系列(1ms差分)
    a = axes[1][1]
    a.plot(t[1:] * 1000, np.diff(x) / RTEX_DT, label="Vx")
    a.plot(t[1:] * 1000, np.diff(y) / RTEX_DT, label="Vy")
    a.plot(t[1:] * 1000, np.diff(z) / RTEX_DT * 10, label="Vz ×10")
    a.set_xlabel("t [ms]"); a.set_ylabel("速度 [mm/s]")
    a.set_title("速度プロファイル")
    a.legend(); a.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def plot_plan_3d(plan: HopPlan, t, x, y, z, r_forbid, out_png, show=False):
    """禁止円筒(z_safe以下)・ウェハ円板・実軌跡の3D表示

    Z軸は物理スケールだと潰れて見えないため誇張表示(軸ラベルに明記)。
    """
    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(projection="3d")
    cx, cy = WAFER_CENTER

    # ウェハ円板 (z=0)
    th = np.linspace(0.0, 2.0 * np.pi, 80)
    rr = np.linspace(0.0, WAFER_RADIUS, 2)
    TH, RR = np.meshgrid(th, rr)
    ax.plot_surface(cx + RR * np.cos(TH), cy + RR * np.sin(TH),
                    np.zeros_like(TH), color="0.6", alpha=0.5,
                    linewidth=0, shade=False)

    # 禁止円筒の側面 (0 <= z <= z_safe, 半径 r_forbid)
    ZZ = np.linspace(0.0, Z_SAFE, 2)
    THc, ZC = np.meshgrid(th, ZZ)
    ax.plot_surface(cx + r_forbid * np.cos(THc), cy + r_forbid * np.sin(THc),
                    ZC, color="orange", alpha=0.25, linewidth=0, shade=False)
    # 円筒の上縁 (z = z_safe)
    ax.plot(cx + r_forbid * np.cos(th), cy + r_forbid * np.sin(th),
            Z_SAFE * np.ones_like(th), "--", color="darkorange", lw=1.5,
            label=f"禁止円筒上縁 z_safe={Z_SAFE}mm")

    # 軌跡: 円内通過区間を色分け
    if plan.t_enter >= 0:
        inside = (t >= plan.t_enter) & (t <= plan.t_exit)
        ax.plot(x[~inside & (t < plan.t_enter)], y[~inside & (t < plan.t_enter)],
                z[~inside & (t < plan.t_enter)], "-", color="tab:green", lw=2,
                label="軌跡(円外)")
        ax.plot(x[inside], y[inside], z[inside], "-", color="tab:red", lw=2,
                label="軌跡(円筒上空通過)")
        ax.plot(x[t > plan.t_exit], y[t > plan.t_exit], z[t > plan.t_exit],
                "-", color="tab:green", lw=2)
    else:
        ax.plot(x, y, z, "-", color="tab:green", lw=2, label="軌跡")

    # 始点・終点
    ax.scatter(*plan.A, color="tab:blue", s=50, label=f"A(カット終了)")
    ax.scatter(*plan.B, color="tab:purple", s=50, label=f"B(次ライン開始)")

    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_zlabel("Z [mm](誇張表示)")
    ax.set_title(f"Zホップ3D軌跡  合計 {plan.total*1000:.1f} ms")
    # XYは等縮尺、Zのみ誇張
    ax.set_box_aspect((1.0, 1.0, 0.35))
    ax.legend(loc="upper left", fontsize=9)
    ax.view_init(elev=28, azim=-75)

    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    if show:
        plt.show()
    plt.close(fig)


# ============================================================
# メイン: 片方向カットの戻りシナリオ
# ============================================================

def main():
    r_forbid = WAFER_RADIUS + XY_MARGIN

    # カット終了点(右端オーバートラベル位置、カット深さ)
    x_end = WAFER_CENTER[0] + WAFER_RADIUS + OVERTRAVEL
    x_start = WAFER_CENTER[0] - WAFER_RADIUS - OVERTRAVEL
    y_line = 0.0  # ウェハ中央のライン(最悪ケース: 通過距離最大)

    A = (x_end, y_line, CUT_DEPTH)
    B = (x_start, y_line + INDEX_PITCH, CUT_DEPTH)

    plan = plan_hop(A, B, WAFER_CENTER, r_forbid, Z_SAFE)
    t, x, y, z = sample_plan(plan)

    out_dir = os.path.dirname(os.path.abspath(__file__))
    write_csv(os.path.join(out_dir, "hop_commands.csv"), t, x, y, z)
    plot_plan(plan, t, x, y, z, r_forbid, os.path.join(out_dir, "hop_plan.png"))
    plot_plan_3d(plan, t, x, y, z, r_forbid,
                 os.path.join(out_dir, "hop_plan_3d.png"),
                 show=("--show" in sys.argv))

    # 安全チェック: 円内で z >= z_safe が守れているか
    rr = np.hypot(x - WAFER_CENTER[0], y - WAFER_CENTER[1])
    violation = np.any((rr < r_forbid - 1e-9) & (z < Z_SAFE - 1e-6))

    print(f"A = {A}  ->  B = {B}")
    print(f"XY距離        : {plan.xy_prof.dist:.1f} mm")
    print(f"XY移動時間    : {plan.xy_prof.T*1000:.1f} ms (開始遅延なし)")
    print(f"Z上昇/下降    : {plan.up_prof.T*1000:.1f} / {plan.down_prof.T*1000:.1f} ms")
    if plan.z_head > 0:
        print(f"Z先行上昇     : A到達の {plan.z_head*1000:.1f} ms 前に開始"
              f"(カット終端オーバートラベル走行がこの時間以上あること)")
    if plan.t_enter >= 0:
        print(f"円内通過区間  : {plan.t_enter*1000:.1f} - {plan.t_exit*1000:.1f} ms")
    if plan.down_tail > 0:
        print(f"Z下降テール   : XY完了後 {plan.down_tail*1000:.1f} ms"
              f"(次ラインの助走中に下降継続可)")
    print(f"安全チェック  : {'NG! 円内で z < z_safe' if violation else 'OK(円内は常に z >= z_safe)'}")
    print(f"── A→B所要(Z下降完了まで)   : {plan.total*1000:.1f} ms")
    print(f"── 実効ライン間ロス(XY律速): {plan.xy_prof.T*1000:.1f} ms "
          f"※Z先行上昇と下降テールを前後のカット走行に重ねた場合")
    print(f"── 参考: 逐次実行           : {plan.total_sequential*1000:.1f} ms")
    print(f"出力: hop_commands.csv ({len(t)}行, t<0はZ先行上昇分), hop_plan.png, hop_plan_3d.png")


if __name__ == "__main__":
    main()
