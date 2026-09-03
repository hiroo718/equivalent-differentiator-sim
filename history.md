# `func_c2d()` の来歴 / History of `func_c2d()`

本リポジトリの `func_c2d()` は、一見すると小さなユーティリティ関数だが、その背景には比較的長い来歴がある。

*At first glance, `func_c2d()` in this repository may look like a small utility function, but it carries a relatively long history behind it.*

源流の一つは、指導教官の研究室で用いられていた外乱オブザーバに遡る。外乱オブザーバでは微分器を用いるが、微分によるノイズの増幅を抑えるため、微分器には1次のローパスフィルタが組み合わされていた。この微分器と1次ローパスフィルタからなる伝達関数を部分分数分解によって等価変換し、等価微分器として定式化する、という知見が得られていた。この考え方は、後に山崎の博士論文にも採用されている。

*One of its roots dates back to the disturbance observer used in my academic advisor's laboratory. Although a disturbance observer utilizes a differentiator, a first-order low-pass filter was combined with it to suppress noise amplification caused by differentiation. An insight was gained by applying partial fraction expansion to equivalently transform the transfer function of this combined differentiator and low-pass filter, formulating it as an "equivalent differentiator." This concept was later incorporated into Yamazaki's doctoral dissertation.*

一方、現在の `func_c2d()` が直接生まれたきっかけは、より実務的な事情だった。会社の開発環境では MATLAB や Control System Toolbox を利用できなかったため、連続時間系の状態方程式を自前で離散化する必要があった。

*On the other hand, the direct catalyst for creating the current `func_c2d()` stemmed from a more practical scenario. Because MATLAB and the Control System Toolbox were unavailable in the company's development environment, it was necessary to discretize continuous-time state-space equations from scratch.*

そこで、森下巌先生・上野敏行先生による『ディジタル制御工学』を繰り返し読み返し、連続時間の状態方程式

*Returning to first principles, I revisited Digital Control Engineering by Prof. Iwao Morishita and Prof. Toshiyuki Ueno, reflecting on the fundamental continuous-time state equation:*

$$\dot{x} = Ax + Bu$$

を離散化する際の

*and its discretized counterparts:*

$$A_d = e^{AT}, \qquad B_d = \int_0^T e^{A\lambda} \, d\lambda \, B$$

という式に立ち返った。そして、これらを求めるための行列指数関数をマクローリン展開し、NumPyだけで直接計算する実装を作った。

*I then implemented a solver using the Maclaurin series expansion of the matrix exponential, computing everything directly using only NumPy.*

実装では、行列指数関数の打ち切り級数を $E$ として一度計算し、それを $A_d$ と $B_d$ の双方で共有するようにした。これにより、同じ $AT$ の冪乗をそれぞれについて個別に計算することを避けている。その後、Gene F. Franklin、J. David Powell、Michael L. Workman による *Digital Control of Dynamic Systems*（第2版、Addison-Wesley, 1990）も購入し、離散制御に関する理解を補強していった。

*In the implementation, the truncated series of the matrix exponential was computed once as $E$ and shared between $A_d$ and $B_d$. This avoided calculating powers of $AT$ twice independently. Later, I acquired Digital Control of Dynamic Systems (2nd ed., Addison-Wesley, 1990) by Gene F. Franklin, J. David Powell, and Michael L. Workman to further solidify my understanding of digital control systems.*

さらに数年を経て、当初は伝達関数として定式化していた等価微分器についても、本リポジトリで用いている簡潔な状態空間表現へと改めて書き直すことになった。これは、MATLABでハイパスフィルタをスマートに実装しようと考え、紙と鉛筆で状態方程式と出力方程式として定式化した成果である。

*Years later, the equivalent differentiator—originally formulated in terms of transfer functions—was rewritten into the clean state-space representation used in this repository. This resulted from an effort to elegantly implement a high-pass filter in MATLAB, deriving its state equation and output equation by hand with pen and paper.*

振り返れば、`func_c2d()` はもともと「手元にライブラリがない」という実務上の制約から生まれた実装だった。しかしその過程で、教科書の式を自分で展開して行列指数関数からの離散化を理解し、さらに後年には等価微分器の理論そのものを状態空間表現から再構成するところまでたどり着いた。

*Looking back, `func_c2d()` was originally born out of a practical constraint: the lack of available libraries. Yet along the way, expanding textbook formulas manually led to a deep understanding of matrix exponential-based discretization, which eventually culminated years later in reconstructing the theory of equivalent differentiators from a state-space perspective.*

小さなユーティリティ関数ではあるが、そこには「既存のライブラリがなくても、必要な理論をゼロから組み立て直せる」ということを実際に試した道標が刻まれている。

*Though it is just a small utility function, it serves as a signpost demonstrating that even without existing libraries, one can rebuild and master the underlying theory from scratch.*
