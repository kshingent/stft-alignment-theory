# 理論と実装の検証用コード

```python
import numpy as np
from scipy.signal import ShortTimeFFT
import math

def verify_stft_alignment(M, h, n):
    """
    SciPy ShortTimeFFT の実装値とアライメント理論定義の一致を検証する。
    """
    # --- 1. 理論に基づく定数と関数の定義 ---
    m_mid = M // 2
    k_off = 0  # ShortTimeFFT の構造的制約（起点固定）

    # 窓の左端 L(p) と 右端 R(p) のインデックス関数
    L = lambda p: h * p + k_off - m_mid
    R = lambda p: h * p + k_off - m_mid + M - 1

    # --- 2. 8つの指標の理論定義式 ---
    # p_min: 窓の右端が信号の開始(0)に最初にかかるスライス
    p_min_th = math.ceil(-(k_off + M - 1 - m_mid) / h)
    
    # p_max: 窓の左端が信号の終端(n-1)を完全に超える最初のスライス
    p_max_th = math.ceil((n + m_mid - k_off) / h)
    
    # p_lb: プリパディングの影響がなくなる（左端が0以上になる）最初のスライス
    p_lb_th = math.ceil((m_mid - k_off) / h)
    
    # p_ub: ポストパディングの影響が始まる（右端がn-1を超える）最初のスライス
    # 修正版：SciPyの判定基準 R(p) >= n に基づく導出
    p_ub_th = math.floor((n - M + m_mid - k_off) / h) + 1

    # サンプルインデックス系の指標
    k_min_th = L(p_min_th)
    k_max_th = R(p_max_th - 1) + 1
    k_lb_th = R(p_lb_th - 1) + 1  # プリパディングを抜けた最初のサンプル
    k_ub_th = L(p_ub_th)         # ポストパディングが始まった最初のサンプル

    # --- 3. ShortTimeFFT 実装値の取得 ---
    SFT = ShortTimeFFT(np.ones(M), hop=h, fs=1)

    p_min_impl = SFT.p_min
    p_max_impl = SFT.p_max(n)
    k_lb_impl, p_lb_impl = SFT.lower_border_end
    k_ub_impl, p_ub_impl = SFT.upper_border_begin(n)
    k_min_impl = SFT.k_min
    k_max_impl = SFT.k_max(n)

    # --- 4. 結果の照合 ---
    results = {
        "p_min": (p_min_th, p_min_impl),
        "p_max": (p_max_th, p_max_impl),
        "p_lb": (p_lb_th, p_lb_impl),
        "p_ub": (p_ub_th, p_ub_impl),
        "k_min": (k_min_th, k_min_impl),
        "k_max": (k_max_th, k_max_impl),
        "k_lb": (k_lb_th, k_lb_impl),
        "k_ub": (k_ub_th, k_ub_impl),
    }

    all_match = True
    print(f"--- Case: M={M}, h={h}, n={n} ---")
    for key, (th, impl) in results.items():
        match = (th == impl)
        all_match &= match
        print(f"{key:5}: Theory={th:4}, Impl={impl:4} -> {'OK' if match else 'NG'}")
    
    return all_match

# 複数のパターンで実行検証（偶数・奇数・長尺）
test_cases = [
    (6, 2, 50), (9, 3, 50), (50, 10, 1000),
    (8, 4, 64), (7, 3, 64), (31, 15, 256), (63, 16, 512)
]

success_count = 0
for M, h, n in test_cases:
    if verify_stft_alignment(M, h, n):
        success_count += 1

print(f"\nVerification Results: {success_count}/{len(test_cases)} cases passed.")

```

# 実行結果
```bash
--- Case: M=6, h=2, n=50 ---
p_min: Theory=  -1, Impl=  -1 -> OK
p_max: Theory=  27, Impl=  27 -> OK
p_lb : Theory=   2, Impl=   2 -> OK
p_ub : Theory=  24, Impl=  24 -> OK
k_min: Theory=  -5, Impl=  -5 -> OK
k_max: Theory=  55, Impl=  55 -> OK
k_lb : Theory=   5, Impl=   5 -> OK
k_ub : Theory=  45, Impl=  45 -> OK
--- Case: M=9, h=3, n=50 ---
p_min: Theory=  -1, Impl=  -1 -> OK
p_max: Theory=  18, Impl=  18 -> OK
p_lb : Theory=   2, Impl=   2 -> OK
p_ub : Theory=  16, Impl=  16 -> OK
k_min: Theory=  -7, Impl=  -7 -> OK
k_max: Theory=  56, Impl=  56 -> OK
k_lb : Theory=   8, Impl=   8 -> OK
k_ub : Theory=  44, Impl=  44 -> OK
--- Case: M=50, h=10, n=1000 ---
p_min: Theory=  -2, Impl=  -2 -> OK
p_max: Theory= 103, Impl= 103 -> OK
p_lb : Theory=   3, Impl=   3 -> OK
p_ub : Theory=  98, Impl=  98 -> OK
k_min: Theory= -45, Impl= -45 -> OK
k_max: Theory=1045, Impl=1045 -> OK
...
k_lb : Theory=  48, Impl=  48 -> OK
k_ub : Theory= 465, Impl= 465 -> OK

Verification Results: 7/7 cases passed.
Output is truncated. View as a scrollable element or open in a text editor. Adjust cell output settings...
```