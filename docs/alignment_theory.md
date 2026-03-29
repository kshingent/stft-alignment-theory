# Alignment Theory Between Two Scales and Formulation of Scale Alignment in STFT

## 1. Alignment Theory Between Scales
This chapter mathematically defines two independent scales ( $t_1$ and $t_2$) on the axis of a physical quantity (e.g., time) and establishes the foundational concepts for evaluating their consistency.

### 1.1 Mathematical Definition of the Time Axis
Two sequences on the physical time axis $\mathbb{R}$ are defined as follows.
*   **Input scale $t_1[k_1]$**: $k_1 \Delta t_1 + t_{1, \mathrm{off}}$ (where $k_1 \in \mathbb{Z}$ is the sample index, $\Delta t_1$ is the sampling interval, and $t_{1, \mathrm{off}}$ is the time offset).
*   **Output scale $t_2[k_2]$**: $k_2 \Delta t_2 + t_{2, \mathrm{off}}$ (where $k_2 \in \mathbb{Z}$ is the sample index, $\Delta t_2$ is the sampling interval, and $t_{2, \mathrm{off}}$ is the time offset).

### 1.2 Non-Dimensionalization and Variable Relationships
To describe the relationship between the two scales, the following dimensionless quantities are defined.
*   **Scale $s_{2 \leftarrow 1}$ of scale 2's sampling interval relative to scale 1's sampling interval**: the ratio of $\Delta t_2$ to $\Delta t_1$ (i.e., $\Delta t_2 = s_{2 \leftarrow 1} \Delta t_1$).
*   **Normalized time offset of scale 2 relative to scale 1, $\Delta \kappa_{{2 \leftarrow 1}}$**: the time difference between the origins of the two scales, normalized by $\Delta t_1$:
$$
\Delta \kappa_{{2 \leftarrow 1}} \coloneqq \frac{t_{2, \mathrm{off}} - t_{1, \mathrm{off}}}{\Delta t_1} (\in \mathbb{R})
$$

Using these quantities, the relationship in index space is derived as follows.
$$\frac{t_2[k_2] - t_1[k_1]}{\Delta t_1} = (s_{2 \leftarrow 1} k_2 - k_1) + \Delta \kappa_{{2 \leftarrow 1}}$$
The index relationship for which the left-hand side equals zero is then:
$$
k_1 = s_{2 \leftarrow 1} k_2 + \Delta \kappa_{2 \leftarrow 1}
$$

This equation is called the **alignment equation**.

### 1.3 Special States: Alignment and Congruence
*   **Harmonic**: $s_{2 \leftarrow 1} \in \mathbb{Z}$. In this case, "scale 2 is harmonic with respect to scale 1." Written as $t_{1} \mid t_{2}$. The converse does not necessarily hold.
*   **Alignment**: Harmonic and $\Delta \kappa_{2 \leftarrow 1} \in \mathbb{Z}$. When $\Delta \kappa_{2 \leftarrow 1} \in \mathbb{Z}$, it is written explicitly as $\Delta \kappa_{{2 \leftarrow 1}} = \Delta k_{{2 \leftarrow 1}}$. In this state, every element of scale 2 exactly coincides with some element of scale 1 on the physical time axis. Written as $t_1 \supseteq t_2$.
*   **Congruent**: Harmonic and further satisfying $s_{2 \leftarrow 1} = 1$ (i.e., $\Delta t_1 = \Delta t_2$). This is a prerequisite for maintaining a fixed phase relationship across the entire scale. Written as $t_{1} \sim t_{2}$.
*   **Coincidence**: Both aligned and congruent. Written as $t_{1} \equiv t_{2}$.
*   **Identical**: Coincident and with equal indices. Written as $t_{1} = t_{2}$.

### 1.4 Relationships Among Three or More Scales
Suppose scale 1 and scale 2 have a certain relationship, and scale 2 and scale 3 also have a certain relationship. We examine what relationship holds between scale 1 and scale 3.

The transitive law holds for all of: harmonic, alignment, congruent, coincidence, and identical. Furthermore, congruent, coincidence, and identical also satisfy reflexivity and symmetry, and are therefore equivalence relations.

Suppose $t_{1}$ is harmonic with respect to $t_{2}$ ( $t_{2} \mid t_{1}$ ) and $t_{3}$ is aligned with respect to $t_{2}$ ( $t_2 \mid t_3$ ). The relationship between $t_{1}$ and $t_{3}$ is not immediately determined, but from the relation
$$
s_{3 \leftarrow 1}= \frac{\Delta t_3}{ \Delta t_1 } = \frac{\Delta t_3}{ \Delta t_2 }\frac{\Delta t_2}{ \Delta t_1 } = \frac{s_{3 \leftarrow 2}}{s_{1 \leftarrow 2}}
$$
and the fact that $s_{1 \leftarrow 2}= \frac{\Delta t_1}{\Delta t_2} \in \mathbb{Z}$ and $s_{3 \leftarrow 2}= \frac{\Delta t_3}{\Delta t_2 } \in \mathbb{Z}$, if $s_{1 \leftarrow 2}$ is a divisor of $s_{3 \leftarrow 2}$, then $s_{3 \leftarrow 1} \in \mathbb{Z}$, meaning $t_{3}$ is harmonic with respect to $t_{1}$ ( $t_{1} \mid t_{3}$ ).


## 2. Application of Alignment Theory to STFT
This chapter applies the alignment concepts defined in Chapter 1 to the Short-Time Fourier Transform (STFT) in signal processing, and formulates the relationships between the time scales of the input signal and each window of the STFT.

### 2.1 Notation and Variable Definitions
The following symbols are defined to describe the STFT parameters and scales.

| Symbol | Definition / Description |
| :--- | :--- |
| $\mathrm{in}$ | Subscript denoting the input time scale |
| $\mathrm{out}$ | Subscript denoting the STFT time scale (window index) |
| $\mathrm{w}[k_{\mathrm{out}}]$ | Subscript denoting the $k_{\mathrm{out}}$-th window |
| $h$ | Hop size of the window (number of samples) |
| $M$ | Number of samples in the window function |
| $m_{mid}$ | Center index of the window ($M // 2$) |

### 2.2 Relationship Between Window Scale and Input Scale in STFT
In STFT, the time scale of each window and the time scale of the input signal are always in a "coincidence" relationship. That is, for any window $k_{\mathrm{out}}$, the following holds:
$$s_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}} = 1, \quad \Delta \kappa_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}} = \Delta k_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}}$$

Since coincidence is an equivalence relation, in STFT, the time scale of any window is always in a "coincidence" relationship with the time scale of any other window. Furthermore, the origin of each window in STFT is arranged as an arithmetic sequence with spacing $\Delta t_{\mathrm{out}} = h \Delta t_{\mathrm{in}}$ according to the hop size $h$. Therefore, the offset time of each window can be written as:
$$
t_{\mathrm{w}[k_{\mathrm{out}}], \mathrm{off}} = k_{\mathrm{out}} \Delta t_{\mathrm{out}} + t_{\mathrm{w}[0], \mathrm{off}}
$$
This is regarded as an independent time scale governing the placement of each window, and is newly defined as the **"offset scale" ($\mathrm{woff}$)**:
$$t_{\mathrm{woff}}[k_{\mathrm{out}}] \coloneqq k_{\mathrm{out}} \Delta t_{\mathrm{out}} + t_{\mathrm{woff, off}}$$
where the initial offset is $t_{\mathrm{woff, off}} \coloneqq t_{\mathrm{w}[0], \mathrm{off}}$. That is, since $\Delta \kappa_{\mathrm{woff} \leftarrow \mathrm{w}[0]} = 0$, we have $t_{\mathrm{w}[0]} \supseteq t_{\mathrm{woff}}$. By the transitive law, $t_{\mathrm{w}[k_{\mathrm{out}}]} \supseteq t_{\mathrm{woff}}$ also holds immediately.

### 2.3 Relationship Between Input Scale and Output Scale in STFT
Using the offset scale, the relationship between the input scale and the output scale can be understood simply. Since the offset scale $t_{\mathrm{woff}}$ and the output time scale $t_{\mathrm{out}}$ have the same spacing, they are "congruent." The offset scale $t_{\mathrm{woff}}$ is "aligned" with the input time scale $t_{\mathrm{in}}$. Therefore, the output time scale is "harmonic" with respect to the input time scale. Whether this constitutes "alignment" depends on whether $\Delta \kappa_{\mathrm{out} \leftarrow \mathrm{in}}$ is an integer.

### 2.4 Index Relationships Among Input, Window, and Output Scales in STFT
From the relationships mediated by the offset scale, we derive the index correspondence equations among the input scale, window scale, and output scale. This equation concisely expresses the correspondence between indices in STFT and is of fundamental importance in computation.

From the alignment equation $k_1 = s_{2 \leftarrow 1} k_2 + \Delta \kappa_{2 \leftarrow 1}$:

$$
k_{\mathrm{in}} = k_{\mathrm{w}[k_{\mathrm{out}}]} + \Delta k_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}}
$$

Here, since the window scale and the input time scale are in a "coincidence" relationship, $\kappa$ is replaced by $k$. Expanding the second term on the right-hand side using the definition:

$$
\Delta k_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}} = \frac{t_{\mathrm{w}[k_{\mathrm{out}}], \mathrm{off}} - t_{\mathrm{in}, \mathrm{off}}}{\Delta t_{\mathrm{in}}}
$$

Since $t_{\mathrm{w}[k_{\mathrm{out}}], \mathrm{off}} = t_{\mathrm{woff}}[k_{\mathrm{out}}]$, substituting gives:

$$
\Delta k_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}} = \frac{ t_{\mathrm{woff}}[k_{\mathrm{out}}] - t_{\mathrm{in}, \mathrm{off}}}{\Delta t_{\mathrm{in}}}
$$

Further expanding $t_{\mathrm{woff}}[k_{\mathrm{out}}]$:

$$
\Delta k_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}} = \frac{ k_{\mathrm{out}} \Delta t_{\mathrm{out}} + t_{\mathrm{woff, off}} - t_{\mathrm{in}, \mathrm{off}}}{\Delta t_{\mathrm{in}}}
$$

The first and second terms in the numerator on the right-hand side correspond to the index difference of the time offsets, so grouping them by definition:

$$
\Delta k_{\mathrm{w}[k_{\mathrm{out}}] \leftarrow \mathrm{in}} = h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}
$$

Here again, $\kappa$ is replaced by $k$ using the fact that the offset scale is aligned with the input time scale. Substituting into the alignment equation:

$$
k_{\mathrm{in}} = k_{\mathrm{w}[k_{\mathrm{out}}]} + h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}
$$

This equation is referred to as the **index correspondence equation in STFT**.

## 3. Summary of Data Boundary Indices in STFT
We define the boundary indices ($k_{\mathrm{out}}$) at the data boundaries of the STFT (Short-Time Fourier Transform). Using the index correspondence equation, we write $k_{\mathrm{in}}(k_{\mathrm{w}[k_{\mathrm{out}}]}, k_{\mathrm{out}})$.

| State | Description | Proposed symbol $k_{\mathrm{out,*}}$ | Condition |
| :--- | :--- | :--- | :--- |
| **Fully outside (left), maximum** | Maximum index where the right edge of the window has not yet reached the start of the data | $k_{\mathrm{out, pre\_max}}$ | $k_{\mathrm{in}}(M-1, k_{\mathrm{out}}) < 0$ |
| **Partially inside (left), minimum** | Minimum index where the right edge of the window just enters the data region | $k_{\mathrm{out, start\_edge}}$ | $k_{\mathrm{in}}(M-1, k_{\mathrm{out}}) \ge 0$ |
| **Fully inside (left), minimum** | Minimum index where the left edge of the window is completely inside the data region | $k_{\mathrm{out, full\_min}}$ | $k_{\mathrm{in}}(0, k_{\mathrm{out}}) \ge 0$ |
| **Fully inside (right), maximum** | Maximum index where the right edge of the window does not exceed the data region | $k_{\mathrm{out, full\_max}}$ | $k_{\mathrm{in}}(M-1, k_{\mathrm{out}}) < N$ |
| **Partially inside (right), maximum** | Maximum index where the left edge of the window is still inside the data region | $k_{\mathrm{out, end\_edge}}$ | $k_{\mathrm{in}}(0, k_{\mathrm{out}}) < N$ |
| **Fully outside (right), minimum** | Minimum index where the left edge of the window has completely passed the data region | $k_{\mathrm{out, post\_min}}$ | $k_{\mathrm{in}}(0, k_{\mathrm{out}}) \ge N$ |

We derive the explicit formulas for each $k_{\mathrm{out,*}}$.

### 2. Computation of Each Boundary Index

#### ① Fully outside (left), maximum: $k_{\mathrm{out, pre\_max}}$
Condition: $k_{\mathrm{in}}(M-1, k_{\mathrm{out}}) < 0$
$$(M-1) + h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} < 0$$
$$h k_{\mathrm{out}} < -(M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}$$
$$k_{\mathrm{out}} < \frac{-(M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h}$$
The **maximum** integer $k_{\mathrm{out}}$ satisfying this condition is:
$$k_{\mathrm{out, pre\_max}} = \left\lfloor \frac{-(M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \right\rfloor$$

#### ② Partially inside (left), minimum: $k_{\mathrm{out, start\_edge}}$
Condition: $k_{\mathrm{in}}(M-1, k_{\mathrm{out}}) \ge 0$
$$(M-1) + h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} \ge 0$$
$$h k_{\mathrm{out}} \ge -(M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}$$
The **minimum** integer $k_{\mathrm{out}}$ satisfying this condition is:
$$k_{\mathrm{out, start\_edge}} = \left\lceil \frac{-(M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \right\rceil$$

#### ③ Fully inside (left), minimum: $k_{\mathrm{out, full\_min}}$
Condition: $k_{\mathrm{in}}(0, k_{\mathrm{out}}) \ge 0$
$$0 + h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} \ge 0$$
$$k_{\mathrm{out}} \ge -\frac{\Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h}$$
The **minimum** integer $k_{\mathrm{out}}$ satisfying this condition is:
$$k_{\mathrm{out, full\_min}} = \left\lceil -\frac{\Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \right\rceil$$

#### ④ Fully inside (right), maximum: $k_{\mathrm{out, full\_max}}$
Condition: $k_{\mathrm{in}}(M-1, k_{\mathrm{out}}) < N$
$$(M-1) + h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} < N$$
$$h k_{\mathrm{out}} < N - 1 - (M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} = N - M - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}$$
The **maximum** integer $k_{\mathrm{out}}$ satisfying this condition is:
$$k_{\mathrm{out, full\_max}} = \left\lfloor \frac{N - M - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \right\rfloor$$

#### ⑤ Partially inside (right), maximum: $k_{\mathrm{out, end\_edge}}$
Condition: $k_{\mathrm{in}}(0, k_{\mathrm{out}}) < N$
$$0 + h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} < N$$
$$h k_{\mathrm{out}} < N - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}$$
The **maximum** integer $k_{\mathrm{out}}$ satisfying this condition is:
$$k_{\mathrm{out, end\_edge}} = \left\lfloor \frac{N - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \right\rfloor$$

#### ⑥ Fully outside (right), minimum: $k_{\mathrm{out, post\_min}}$
Condition: $k_{\mathrm{in}}(0, k_{\mathrm{out}}) \ge N$
$$0 + h k_{\mathrm{out}} + \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}} \ge N$$
$$h k_{\mathrm{out}} \ge N - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}$$
The **minimum** integer $k_{\mathrm{out}}$ satisfying this condition is:
$$k_{\mathrm{out, post\_min}} = \left\lceil \frac{N - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \right\rceil$$

---

### 3. Summary of Results

| Boundary name | Formula (with $\Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}$) |
| :--- | :--- |
| $k_{\mathrm{out, pre\_max}}$ | $\lfloor \frac{-(M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \rfloor$ |
| $k_{\mathrm{out, start\_edge}}$ | $\lceil \frac{-(M-1) - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \rceil$ |
| $k_{\mathrm{out, full\_min}}$ | $\lceil -\frac{\Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \rceil$ |
| $k_{\mathrm{out, full\_max}}$ | $\lfloor \frac{N - M - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \rfloor$ |
| $k_{\mathrm{out, end\_edge}}$ | $\lfloor \frac{N - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \rfloor$ |
| $k_{\mathrm{out, post\_min}}$ | $\lceil \frac{N - \Delta k_{\mathrm{woff} \leftarrow \mathrm{in}}}{h} \rceil$ |

From the index correspondence equation, $k_{\mathrm{in}}(M-1, k_{\mathrm{out,*}})$ and $k_{\mathrm{in}}(0, k_{\mathrm{out,*}})$, corresponding to both edges of window $\mathrm{w}[k_{\mathrm{out}}]$, can also be computed immediately.
