# 等価微分器 
# equivalent-differentiator-sim
# s/(Ts+1)
# by Dr. Hiroo Yamazaki
# Aug. 12th 2026
# Aug. 30th 2026 updated

# -*- coding : utf-8 -*-

from scipy import signal
from scipy.fft import fft, fftfreq
from scipy.integrate import odeint
from scipy.linalg import *
from matplotlib.pyplot import figure, plot, grid, show, subplot,\
    xlabel, ylabel, ylim, semilogx, xlim, legend
from numpy import *

def func_c2d(A,B,dt):
    AT=A*dt
    n=len(A)
    E=eye(n)+1/2*AT+1/(2*3)*AT.dot(AT)+1/(2*3*4)*AT.dot(AT).dot(AT)\
    +1/(2*3*4*5)*AT.dot(AT).dot(AT).dot(AT)\
    +1/(2*3*4*5*6)*AT.dot(AT).dot(AT).dot(AT).dot(AT)
    ad=eye(n)+E.dot(AT)
    bd=E.dot(dt).dot(B)
    return [ad, bd]

N=2**12
fm=10
tau=1/fm
dt=tau/5
A=-1/tau*eye(4)
B= 1/tau*eye(4)
C=A
D=B
[ad, bd]=func_c2d(A,B,dt)
tmax=(N-1)*dt
t=arange(0,tmax,dt)
x=zeros((4,1))
u=zeros((4,1))
#u[0]=random.randn(1)
y=C.dot(x)+D.dot(u)
dat=zeros((len(t),10))
for i in arange(0,len(t)):
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
    
n = len(dat)
n= 4*1024
fs = 1/dt
window='hann'
nperseng=n//4
#nperseng=2
noverlap=2*256
#noverlap=1
nfft=1024
f, P = signal.welch(dat[:,1], fs, window, nperseg=nperseng,noverlap=noverlap,nfft=nfft)
f1, P1 = signal.welch(dat[:,0], fs, window, nperseg=nperseng,noverlap=noverlap,nfft=nfft)

figure(1)
subplot(211)
plot(t,dat[:,0],'-k',t,dat[:,1],'-r')
grid()
ylabel('input[-]')
subplot(212)
plot(t,dat[:,1],'-r')
grid()
ylabel('ourtput[-]')
xlabel('time[s]')
figure(2)
semilogx(f,20*log10(abs(P)),'-r',f1,20*log10(abs(P1/P1)),'-k')
#ylim([-80, 80])
xlim([1e-1,20*1e3])
grid()
xlabel("Frequency[Hz]")
ylabel("Power/frequency[dB/Hz]")
show()
