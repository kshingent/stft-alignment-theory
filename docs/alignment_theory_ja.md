# 2つの目盛間の整列理論とSTFTにおける目盛間の整列の定式化

## 1. 目盛間の整列理論
本章では、物理量（例えば時間）の軸上にある2つの独立した目盛（ $t_{1}$ と $t_{2}$）を数学的に定義し、その整合性を評価するための基礎概念を規定する。

### 1.1 時間軸の数学的定義
物理時間軸 $\mathbb{R}$ 上の2つの数列を以下の通り定義する。
*   **入力側の目盛 $t_{1}[k_{1}]$**: $k_{1} \Delta t_{1} + t_{1, \mathrm{off}}$（$k_{1} \in \mathbb{Z}$ はサンプルのインデックス、$\Delta t_{1}$ はサンプリング間隔、$t_{1, \mathrm{off}}$ は時刻のオフセット）。
*   **出力側の目盛 $t_{2}[k_{2}]$**: $k_{2} \Delta t_{2} + t_{2, \mathrm{off}}$（$k_{2} \in \mathbb{Z}$ はサンプルのインデックス、$\Delta t_{2}$ はサンプリング間隔、$t_{2, \mathrm{off}}$ は時刻のオフセット）。

### 1.2 無次元化と変数関係
2つの目盛の関係性を記述するため、以下の無次元量を定義する。
*   **目盛1のサンプリング間隔に対する目盛2のサンプリング間隔のスケール $s_{2 \leftarrow 1}$**: 目盛2のサンプリング間隔 $\Delta t_{2}$ と目盛1のサンプリング間隔 $\Delta t_{1}$ の比（$\Delta t_{2} = s_{2 \leftarrow 1} \Delta t_{1}$）。
*   **目盛1のサンプリング間隔で正規化した、目盛2の時刻オフセット $\Delta \kappa_{{2 \leftarrow 1}}$**: 2つの時間目盛の起点間の時間差を $\Delta t_{1}$ で正規化した値

$$
\Delta \kappa_{{2 \leftarrow 1}} \coloneqq \frac{t_{2, \mathrm{off}} - t_{1, \mathrm{off}}}{\Delta t_{1}} (\in \mathbb{R})
$$

これらを用いると、インデックス空間における関係式は以下のように導出される。

$$
\frac{t_{2}[k_{2}] - t_{1}[k_{1}]}{\Delta t_{1}} = (s_{2 \leftarrow 1} k_{2} - k_{1}) + \Delta \kappa_{{2 \leftarrow 1}}
$$

つまり左辺が0になるインデックスの関係式は以下のようになる。

$$
k_{1} = s_{2 \leftarrow 1} k_{2} + \Delta \kappa_{2 \leftarrow 1}
$$

この方程式を整列（Alignment）方程式という。

### 1.3 特殊状態：整列と合同
*   **調和的 (Harmonic)**: $s_{2 \leftarrow 1} \in \mathbb{Z}$。このとき、「目盛2は目盛1に」調和的であるという。$t_{1} \mid t_{2}$ と書く。逆は成り立たない。
*   **整列 (Alignment)**: 調和的かつ $\Delta \kappa_{2 \leftarrow 1} \in \mathbb{Z}$。$\Delta \kappa_{2 \leftarrow 1} \in \mathbb{Z}$であることを明示的に $\Delta \kappa_{{2 \leftarrow 1}} = \Delta k_{{2 \leftarrow 1}}$ と書く。このとき、目盛2の全要素が目盛1のいずれかの要素と物理時間軸上で完全に一致する。$t_{1} \supseteq t_{2}$ と書く。
*   **合同 (Congruent)**: 調和的なうえでさらに $s_{2 \leftarrow 1} = 1$ （すなわち $\Delta t_{1} = \Delta t_{2}$）である状態。目盛全体にわたって固定された位相関係を保持するための前提条件となる。$t_{1} \sim t_{2}$ と書く。
*   **一致 (Coincidence)**: 整列かつ合同な状態。$t_{1} \equiv t_{2}$ と書く。
*   **恒等 (Identical)**: 一致かつインデックスも等しい状態。$t_{1} = t_{2}$ と書く。

### 1.4 3つ以上の目盛間の関係性
目盛1と2にある関係性があり、目盛2と3にある関係性が成り立っているとする。このとき目盛1と3の間にどのような関係性が成り立つか整理する。

調和的、整列、合同、一致、恒等のすべてにおいて推移律が成り立つ。また、合同、一致、恒等は反射律と対称律も成り立つから同値関係である。

$t_{1}$ が $t_{2}$ に対して調和的（ $t_{2} \mid t_{1}$ ）かつ $t_{3}$ が $t_{2}$ に対して整列（ $t_{2} \mid t_3$ ）しているとする。この時、$t_{1}$ と $t_{3}$ の関係性は直ちには決まらないが、

$$
s_{3 \leftarrow 1}= \frac{\Delta t_3}{ \Delta t_{1} } = \frac{\Delta t_3}{ \Delta t_{2} }\frac{\Delta t_{2}}{ \Delta t_{1} } = \frac{s_{3 \leftarrow 2}}{s_{1 \leftarrow 2}}
$$

という関係が成り立つことと、$s_{1 \leftarrow 2}= \frac{\Delta t_{1}}{\Delta t_{2}} \in \mathbb{Z},s_{3 \leftarrow 2}= \frac{\Delta t_3}{\Delta t_{2} } \in \mathbb{Z}$ であることから、$s_{1 \leftarrow 2}$ が $s_{3 \leftarrow 2}$ の約数であれば、$s_{3 \leftarrow 1} \in \mathbb{Z}$ すなわち、$t_{3}$ が $t_{1}$ に対して調和的（ $t_{1} \mid t_{3}$ ）である。


## 2. STFTへの整列理論の適用
本章では、第1章で定義した整列概念を信号処理における短時間フーリエ変換（STFT）に適用し、入力信号とSTFTの各窓における時間目盛の関係性を定式化する。

### 2.1 記号と変数名の定義
STFTの各パラメータおよび目盛を記述するため、以下の記号を定義する。

| 数理記号 | 定義・説明 |
| :--- | :--- |
| $\mathrm{in}$ | 入力時間目盛を表す添え字 |
| $\mathrm{out}$ | STFTの時間目盛（窓のインデックス）を表す添え字 |
| $\mathrm{w}[k_{\mathrm{out}}]$ | $k_{\mathrm{out}}$ 番目の窓を表す添え字 |
| $h$ | 窓のホップサイズ（サンプル数） |
| $M$ | 窓関数のサンプル数 |
| $m_{mid}$ | 窓の中心インデックス（$M // 2$） |

### 2.2 STFTにおける窓目盛と入力目盛の関係
STFTにおいて、各窓の時間目盛と入力信号の時間目盛は、常に「一致」の関係にある。すなわち、任意の窓 $k_{\mathrm{out}}$ において以下が成立する。

$$
s_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}} = 1, \quad \Delta \kappa_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}} = \Delta k_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}}
$$

一致関係は同値関係であることから、STFTにおいて、窓の時間目盛と別の窓の時間目盛は、常に「一致」関係にある。さらに、STFTにおける各窓の起点は、ホップサイズ $h$ に応じた間隔 $\Delta t_{\mathrm{out}} = h \Delta t_{\mathrm{in}}$ の等差数列として並ぶ。したがって、各窓のオフセット時間は次のように記述できる。

$$
t_{\mathrm{w}[k_{\mathrm{out}}], \mathrm{off}} = k_{\mathrm{out}} \Delta t_{\mathrm{out}} + t_{\mathrm{w}[0], \mathrm{off}}
$$

これを各窓の配置を規定する独立した時間目盛とみなし、新たに「オフセット目盛（$\mathrm{woff}$）」として定義する。

$$
t_{\mathrm{woff}}[k_{\mathrm{out}}] \coloneqq k_{\mathrm{out}} \Delta t_{\mathrm{out}} + t_{\mathrm{woff, off}}
$$

ここで、初期オフセットは $t_{\mathrm{woff, off}} \coloneqq t_{\mathrm{w}[0], \mathrm{off}}$ である。すなわち、$\Delta \kappa_{\mathrm{woff} \leftarrow \mathrm{w}[0]} = 0 $であるから、$t_{\mathrm{w}[0]} \supseteq t_{\mathrm{woff}} $である。推移律から $t_{\mathrm{w}[k_{\mathrm{out}}]} \supseteq t_{\mathrm{woff}} $も直ちに成り立つ。

### 2.3 STFTにおける入力目盛と出力目盛の関係
オフセット目盛を用いると入力目盛と出力目盛の関係を簡単に理解することができる。オフセット目盛 $t_{\mathrm{woff}}$ と出力時間目盛 $t_{\mathrm{out}}$ の目盛間隔は同じであるから「合同」である。入力時間目盛 $t_{\mathrm{in}}$ に対してオフセット目盛 $t_{\mathrm{woff}}$ は「整列」している。したがって、入力時間目盛に対して出力時間目盛は「調和的」である。これが「整列」であるかどうかは、$\Delta \kappa_{\mathrm{out} \leftarrow \mathrm{in}}$ が整数かどうかで決まる。

### 2.4 STFTにおける入力目盛、窓の目盛、出力目盛のインデックス間の関係
オフセット目盛を介した関係性から、入力目盛、窓の目盛、出力目盛のインデックス間の関係式を導出する。この方程式はSTFTのインデックス間の対応関係を簡潔に示した式であり、計算を行う上で非常に重要である。
整列方程式 $k_{1} = s_{2 \leftarrow 1} k_{2} + \Delta \kappa_{2 \leftarrow 1}$ から

$$
k_{\mathrm{in}} = k_{\mathrm{w}[k_{\mathrm{out}}]} + \Delta k_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}}
$$

である。ただし、窓の目盛と入力時間目盛が「一致」関係にあるため、 $\kappa \rightarrow k$ とした。右辺第二項を定義に基づいて展開して

$$
\Delta k_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}} = \frac{t_{\mathrm{w}[k_{\mathrm{out}}], \mathrm{off}} - t_{\mathrm{in}, \mathrm{off}}}{\Delta t_{\mathrm{in}}}
$$

となる。 $t_{\mathrm{w}[k_{\mathrm{out}}], \mathrm{off}} = t_{\mathrm{woff}}[k_{\mathrm{out}}]$ だから、代入すると以下を得る。

$$
\Delta k_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}} = \frac{ t_{\mathrm{woff}}[k_{\mathrm{out}}] - t_{\mathrm{in}, \mathrm{off}}}{\Delta t_{\mathrm{in}}}
$$

さらに $t_{\mathrm{woff}}[k_{\mathrm{out}}]$ を展開して

$$
\Delta k_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}} = \frac{ k_{\mathrm{out}} \Delta t_{\mathrm{out}} + t_{\mathrm{woff, off}} - t_{\mathrm{in}, \mathrm{off}}}{\Delta t_{\mathrm{in}}}
$$

となる。右辺分母の第一項と第二項は時間オフセットのインデックスの差分に対応しているから定義に基づいてひとまとめにすると

$$
\Delta k_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}} = h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}
$$

となる。ここでも、オフセット目盛が入力時間目盛に対して整列していることを利用して、$\kappa \rightarrow k$ とした。整列方程式に代入して

$$
k_{\mathrm{in}} = k_{\mathrm{w}[k_{\mathrm{out}}]} + h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}
$$

となる。この式をSTFTにおけるインデックス対応方程式と呼ぶことにする。

## 3. STFTにおけるデータ境界Indexの整理
STFT（短時間フーリエ変換）のデータ境界におけるインデックス（$k_{\mathrm{out}}$）を定義する。インデックス対応方程式から $k_{\mathrm{in}}(k_{\mathrm{w}[k_{\mathrm{out}}]}, k_{\mathrm{out}})$と書くことにする。

| 状態 | 定義の説明 | 提案記号 $k_{\mathrm{out,*}}$ | 数式条件 |
| :--- | :--- | :--- | :--- |
| **完全外（左）の最大** | 窓の右端がデータの開始点にすら達していない最大値 | $k_{\mathrm{out, pre\_max}}$ | $k_{\mathrm{in}}(M-1, k_{\mathrm{out}}) < 0$ |
| **一部内（左）の最小** | 窓の右端がデータ領域に入り始めた瞬間 | $k_{\mathrm{out, pre\_max}}$ | $k_{\mathrm{in}}(M-1, k_{\mathrm{out}}) \ge 0$ |
| **完全内（左）の最小** | 窓の左端がデータ領域に完全に入った最小値 | $k_{\mathrm{out, full\_min}}$ | $k_{\mathrm{in}}(0, k_{\mathrm{out}}) \ge 0$ |
| **完全内（右）の最大** | 窓の右端がデータ領域を超えない最大値 | $k_{\mathrm{out, full\_max}}$ | $k_{\mathrm{in}}(M-1, k_{\mathrm{out}}) < N$ |
| **一部内（右）の最大** | 窓の左端がまだデータ領域内に留まっている最大値 | $k_{\mathrm{out, end\_edge}}$ | $k_{\mathrm{in}}(0, k_{\mathrm{out}}) < N$ |
| **完全外（右）の最小** | 窓の左端がデータ領域を完全に通り過ぎた最小値 | $k_{\mathrm{out, post\_min}}$ | $k_{\mathrm{in}}(0, k_{\mathrm{out}}) \ge N$ |

$k_{\mathrm{out,*}}$ の具体的な計算式を導出する。

### 2. 各境界インデックスの計算

#### ① 完全外（左）の最大：$k_{\mathrm{out, pre\_max}}$
条件：$k_{\mathrm{in}}(M-1, k_{\mathrm{out}}) < 0$

$$
(M-1) + h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} < 0
$$

$$
h k_{\mathrm{out}} < -(M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}
$$

$$
k_{\mathrm{out}} < \frac{-(M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h}
$$

この条件を満たす**最大**の整数 $k_{\mathrm{out}}$ は：

$$k_{\mathrm{out, pre\_max}} = \left\lfloor \frac{-(M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \right\rfloor
$$

#### ② 一部内（左）の最小：$k_{\mathrm{out, start\_edge}}$
条件：$k_{\mathrm{in}}(M-1, k_{\mathrm{out}}) \ge 0$

$$(M-1) + h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} \ge 0
$$

$$
h k_{\mathrm{out}} \ge -(M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}
$$

この条件を満たす**最小**の整数 $k_{\mathrm{out}}$ は：

$$
k_{\mathrm{out, start\_edge}} = \left\lceil \frac{-(M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \right\rceil
$$

#### ③ 完全内（左）の最小：$k_{\mathrm{out, full\_min}}$
条件：$k_{\mathrm{in}}(0, k_{\mathrm{out}}) \ge 0$

$$
0 + h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} \ge 0
$$

$$
k_{\mathrm{out}} \ge -\frac{\Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h}
$$

この条件を満たす**最小**の整数 $k_{\mathrm{out}}$ は：

$$
k_{\mathrm{out, full\_min}} = \left\lceil -\frac{\Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \right\rceil
$$

#### ④ 完全内（右）の最大：$k_{\mathrm{out, full\_max}}$
条件：$k_{\mathrm{in}}(M-1, k_{\mathrm{out}}) < N$

$$
(M-1) + h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} < N
$$

$$
h k_{\mathrm{out}} < N - 1 - (M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} = N - M - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}
$$

この条件を満たす**最大**の整数 $k_{\mathrm{out}}$ は：

$$
k_{\mathrm{out, full\_max}} = \left\lfloor \frac{N - M - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \right\rfloor
$$

#### ⑤ 一部内（右）の最大：$k_{\mathrm{out, end\_edge}}$
条件：$k_{\mathrm{in}}(0, k_{\mathrm{out}}) < N$

$$
0 + h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} < N
$$

$$
h k_{\mathrm{out}} < N - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}
$$

この条件を満たす**最大**の整数 $k_{\mathrm{out}}$ は：

$$
k_{\mathrm{out, end\_edge}} = \left\lfloor \frac{N - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \right\rfloor
$$

#### ⑥ 完全外（右）の最小：$k_{\mathrm{out, post\_min}}$
条件：$k_{\mathrm{in}}(0, k_{\mathrm{out}}) \ge N$

$$
0 + h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} \ge N
$$

$$
h k_{\mathrm{out}} \ge N - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}
$$

この条件を満たす**最小**の整数 $k_{\mathrm{out}}$ は：

$$
k_{\mathrm{out, post\_min}} = \left\lceil \frac{N - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \right\rceil
$$

---

### 3. 計算結果まとめ

| 境界の名称 | 計算式（$\Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} = \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}$ とする） |
| :--- | :--- |
| $k_{\mathrm{out, pre\_max}}$ | $\lfloor \frac{-(M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \rfloor$ |
| $k_{\mathrm{out, start\_edge}}$ | $\lceil \frac{-(M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \rceil$ |
| $k_{\mathrm{out, full\_min}}$ | $\lceil -\frac{\Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \rceil$ |
| $k_{\mathrm{out, full\_max}}$ | $\lfloor \frac{N - M - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \rfloor$ |
| $k_{\mathrm{out, end\_edge}}$ | $\lfloor \frac{N - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \rfloor$ |
| $k_{\mathrm{out, post\_min}}$ | $\lceil \frac{N - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \rceil$ |

インデックス対応方程式から窓 $\mathrm{w}[k_{\mathrm{out}}]$ の両端に対応する $k_{\mathrm{in}}(M-1, k_{\mathrm{out,*}}), k_{\mathrm{in}}(0, k_{\mathrm{out,*}})$ も直ちに計算できる。