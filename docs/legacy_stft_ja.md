# 旧`scipy.signal.stft`における入力時間目盛のインデックス  $k_2$  範囲と境界処理の理解
本資料は、旧 `scipy.signal.stft`（以下、旧実装）における出力時間軸決定ロジックを、STFTにおける目盛間の整列理論および `scipy/signal/_spectral.py` の内部実装に基づき定式化したものである。特に、整数演算による「非対称なアンカー（Asymmetric Anchoring）」が物理的な時間アライメントに与える影響と、処理長  $N_{\mathrm{padded}}$  および総窓数  $P$  の派生関係を厳密に定義することを目的とする。

## 2. 指示関数  $\delta_\mathrm{b}$  と離散アンカーポイント  $k_{\mathrm{off}}^{(1,w)}$  の定式化
旧実装では、境界条件 `boundary` の設定によって窓の初期配置（離散インデックス空間上の位置）を決定する。

### 2.1 指示関数  $\delta_\mathrm{b}$
引数 `boundary` に基づき、拡張処理の有無を制御する指示関数を次のように定義する。

- $\mathrm{boundary} = \mathrm{None}$ の場合：
$\delta_\mathrm{b} = 1$

- $\mathrm{boundary} \in \{\text{'zeros'}, \text{'reflect'}, \dots\}$ の場合：
  $\delta_\mathrm{b} = 0$

### 2.2 窓起点オフセット  $k_{\mathrm{off}}^{(1,w)}$
入力信号の起点に対する、初期フレーム（ $k_2=0$ ）における窓の左端（ $m=0$ ）の正規化オフセット  $k_{\mathrm{off}}^{(1,w)}$  は`scipy.signal.stft`においては以下のように定められる。  $$k_{\mathrm{off}}^{(1,w)} = -(1 - \delta_\mathrm{b}) \cdot (M // 2) \tag{4}$$  ここで  $M$  は窓長（`nperseg`）である。
実装上の制約（量子化制約） :  $M // 2$  は整数切り捨て除算である。この「離散アンカーポイント」の決定に整数演算を用いる設計が、奇数長窓において物理的な中心（ $M/2$ ）との間に 時間的量子化誤差 を生じさせる原因となる。

## 3. 最終処理長  $N_{\mathrm{padded}}$  の導出プロセス
旧実装は、実際のFFT計算に先立ち、信号長を「境界拡張」と「フィッティングパディング」の二段階で確定させる。

### 3.1 境界拡張後の長さ  $N_{\mathrm{ext}}$
`boundary` が None でない場合、窓の中心を信号の両端に合わせるため、入力信号  $N$  に対して両端に  $M//2$  ずつの拡張が行われる。  $$N_{\mathrm{ext}} = N + (1 - \delta_\mathrm{b}) \cdot 2 \cdot (M // 2) \tag{5}$$
### 3.2 パディング後の最終処理長（Processing Length）  $N_{\mathrm{padded}}$
`padded` = True の場合、拡張された信号  $N_{\mathrm{ext}}$  をホップサイズ  $h$  で整数個のセグメントに収めるため、右端に  $N_{\mathrm{add}}$  サンプルが追加される。内部ロジックに基づく  $N_{\mathrm{add}}$  の算出式は以下の通りである。  $$N_{\mathrm{add}} = (-(N_{\mathrm{ext}} - M) \pmod h) \pmod h \tag{6}$$  これにより、最終的な処理長  $N_{\mathrm{padded}}$  が決定される。  $$N_{\mathrm{padded}} = N_{\mathrm{ext}} + N_{\mathrm{add}} \tag{7}$$

## 4. 出力時間インデックス  $k_2$  の範囲と物理時刻のマッピング
確定した  $N_{\mathrm{padded}}$  に基づき、出力フレーム数  $P$  および各インデックスの物理時刻を導出する。

### 4.1 総フレーム数  $P$  の決定
$N_{\mathrm{padded}}$  の定義  $(N_{\mathrm{padded}} - M) \pmod h = 0$  より、生成される総フレーム数  $P$  は以下の整数演算で一意に定まる。  $$P = (N_{\mathrm{padded}} - M) // h + 1 \tag{8}$$  したがって、出力時間インデックス集合は  $k_2 \in \{0, 1, \dots, P-1\}$  となる。

### 4.2 物理時刻  $t_2k_2$  への写像
各インデックス  $k_2$  に対応する物理時刻は、以下の物理オフセット  $t_{2,\mathrm{off}}$  を基準に計算される。  $$t_{2,\mathrm{off}} = t_{1,\mathrm{off}} + \delta_\mathrm{b} \cdot (M/2) \cdot \Delta t_1 \tag{9}$$   $$t_2k_2 = k_2 (h \Delta t_1) + t_{2,\mathrm{off}}$$  ここで、式(9)の  $M/2$  は浮動小数点演算である。
