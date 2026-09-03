# Equivalent Differentiator Simulation (SISO State-Space Model)

[English]
In practical simulations and Model-Based Development (MBD), an ideal differentiator (G(s)=s) is non-proper (physically unrealizable) and carries the risk of infinitely amplifying high-frequency noise. This repository provides a mathematical model of an **"Equivalent Differentiator (Approximate Differentiator)"** with a first-order lag element to solve this problem, along with its discretization simulation implemented in Python.

[日本語]
実務のシミュレーションやモデルベース開発（MBD）において、理想微分（G(s)=s）は非プロパー（物理的に実現不可能）であり、高周波ノイズを無限大に増幅するリスクを孕んでいます。本リポジトリでは、この課題を解決するために1次遅れ要素を付加した**「等価微分器（近似微分器）」**の数理モデルの構築、およびPythonによる離散化シミュレーションの実装例をまとめています。

---

## 1. Mathematical Model (数理モデル)

[English]
The transfer function G(s) of an equivalent differentiator with a first-order lag of time constant τ is expressed as follows:
[日本語]
時定数 τ の1次遅れを持つ等価微分器の伝達関数 G(s) は以下の通りです。

\[G(s) = \frac{s}{\tau s + 1}\]

[English] Expanding this into a Single-Input Single-Output (SISO) state-space model yields the following state and output equations.
[日本語] これを1入力1出力（SISO）の状態空間モデルに展開すると、以下の状態方程式および出力方程式が得られます。

### State Equation (状態方程式)
\[\dot{x}(t) = -\frac{1}{\tau} x(t) + u(t)\]

### Output Equation (出力方程式)
\[y(t) = -\frac{1}{\tau^2} x(t) + \frac{1}{\tau} u(t)\]

*(Note: The configuration of coefficients may vary depending on the choice of state variables, but the system remains equivalent. / ※状態変数の取り方により係数の配置は変わりますが、システムとしては等価です)*

---

## 2. Comparison & Discussion: Differentiator vs. HPF
## (ハイパスフィルタとの比較考察)

[English]
Comparing the transfer function of the high-pass filter and that of the equivalent differentiator reveals that their gains differ by a factor of the filter time constant $\tau$.Alternative phrasing:"If you compare the transfer functions of the high-pass filter and the equivalent differentiator, you will notice that their gains differ by the filter time constant $\tau$."
[日本語]
ハイパスフィルタの伝達関数と等価微分器の伝達関数を比較すると、フィルタ時定数τ分ゲインが異なることに気が付くだろう。

---

## 3. Python Implementation (離散化シミュレーション実装)

[English]
This production-ready script directly discretizes (C2D conversion) continuous-time state-space matrices using a custom 6th-order Taylor expansion without relying on built-in libraries like `scipy.signal`. It evaluates step responses and power spectral density (PSD) via Welch's method.

[日本語]
ライブラリ（`scipy.signal`のC2D等）に依存せず、連続系の状態方程式からC2D変換（離散化）のアルゴリズム（テイラー展開6次項まで）を自作し、ステップ応答およびウェルチ法（Welch's method）による周波数解析を行うガチ仕様のコードです。

## History
For the theoretical background and story behind `func_c2d()`, see [History of func_c2d()](/HISTORY.md).

## 🛠️ Credits (クレジット)

[English]
This code was basically written by Hiroo Yamazaki, with comments and the README created with the support of AI.

[日本語]
このコードは基本的に山崎が作成し、コメント文とREADMEの作成はAIのサポートを得ました。
```python
# 等価微分器 
# equivalent-differentiator-sim
# s/(Ts+1)
# by Dr. Hiroo Yamazaki
# Aug. 12th 2026
# Aug. 30th 2026 updated

from scipy import signal
from scipy.fft import fft, fftfreq
from scipy.integrate import odeint
from scipy.linalg import *
from matplotlib.pyplot import figure, plot, grid, show, subplot,\
    xlabel, ylabel, ylim, semilogx, xlim, legend
from numpy import *

def func_c2d(A,B,dt):
    """
    Discretize continuous state-space matrices (A, B) using Taylor expansion.
    """
    AT=A*dt
    n=len(A)
    E=eye(n)+1/2*AT+1/(2*3)*AT.dot(AT)+1/(2*3*4)*AT.dot(AT).dot(AT)\
    +1/(2*3*4*5)*AT.dot(AT).dot(AT).dot(AT)\
    +1/(2*3*4*5*6)*AT.dot(AT).dot(AT).dot(AT).dot(AT)
    ad=eye(n)+E.dot(AT)
    bd=E.dot(dt).dot(B)
    return [ad, bd]

# Simulation parameters
N=2**12
fm=10
tau=1/fm
dt=tau/5

# Continuous-time state-space models (4D structure for extensibility)
A=-1/tau*eye(4)
B= 1/tau*eye(4)
C=A
D=B

# C2D Conversion
[ad, bd]=func_c2d(A,B,dt)

tmax=(N-1)*dt
t=arange(0,tmax,dt)
x=zeros((4,1))
u=zeros((4,1))
y=C.dot(x)+D.dot(u)
dat=zeros((len(t),10))

# Time-step loop
for i in arange(0,len(t)):
    # Apply step input between 1.0s and 2.0s
    if i*dt<=2.0 and i*dt >1.0:
        u[0]=1.0
    else:
        u[0]=0.0
        
    y=C.dot(x)+D.dot(u)
    x=ad.dot(x)+bd.dot(u)

    dat[i,0]=u[0]
    dat[i,1]=y[0]
    dat[i,2]=y[1]
    dat[i,3]=y[2]
    
# Frequency analysis using Welch's method
n = len(dat)
n = 4*1024
fs = 1/dt
window='hann'
nperseng=n//4
noverlap=2*256
nfft=1024

f, P = signal.welch(dat[:,1], fs, window, nperseg=nperseng, noverlap=noverlap, nfft=nfft)
f1, P1 = signal.welch(dat[:,0], fs, window, nperseg=nperseng, noverlap=noverlap, nfft=nfft)

# Plot results
figure(1)
subplot(211)
plot(t,dat[:,0],'-k',t,dat[:,1],'-r')
grid()
ylabel('input[-]')

subplot(212)
plot(t,dat[:,1],'-r')
grid()
ylabel('output[-]')
xlabel('time[s]')

figure(2)
semilogx(f,20*log10(abs(P)),'-r',f1,20*log10(abs(P1/P1)),'-k')
xlim([1e-1,20*1e3])
grid()
xlabel("Frequency[Hz]")
ylabel("Power/frequency[dB/Hz]")
show()
```
