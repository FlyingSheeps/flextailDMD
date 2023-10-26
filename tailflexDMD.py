# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 15:57:46 2023

@author: ko-ko
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 23:30:43 2023

@author: nakamura_ryosuke
"""

import scipy.integrate as sol
import numpy as np
import matplotlib.pyplot as plt
import control as mt



"""定数"""
#環境定数
rho = 1.2
U0 = 8.0
g = 9.8
L = 5.0

"""数値計算用設定"""
#分割数
N = 20
#分割位置
x_list = np.linspace(0,L,N)
dx = x_list[1] - x_list[0]

"""尾翼パラメータ"""
#尾翼面積
St = 1.48
#尾翼揚力傾斜
at = 2.0*np.pi
#尾翼重量
mT = 1.4



sigma = 1.2/L
#曲げ剛性（弾性率*断面二次モーメント）
EJ_list = (54000-9000)*(L-x_list)/L + 9000
#疑似構造減衰
dump = 10*10**(-5)



"""これはたわみねじれがあるけど，制御が入ってないやつ"""
"""動力学部分"""
def sim(t,x):
    
    """変数を展開"""
    du_list= x[0:N]
    u_list = x[N:2*N]
    
    """出力用変数"""
    utt_list = np.zeros_like(u_list)
    
    """片持ち境界条件"""
    u_list[0] = 0.0
    u_list[1] = 0.0
    
    """たわみを二階微分して，剛性をかけて曲げ弾性モーメントにする"""
    Me_list = np.zeros_like(u_list)
    for n in range(N-2):
        up1 = u_list[n+2]
        u0 = u_list[n+1]
        um1 = u_list[n]
        dup1 = du_list[n+2]
        du0 = du_list[n+1]
        dum1 = du_list[n]
        Me_list[n+1] = EJ_list[n+1]*(up1-2.0*u0+um1)/(dx**2) + EJ_list[n+1]*dump*(dup1-2.0*du0+dum1)/(dx**2)
    
    """尾翼にかかる揚力を計算する"""
    Lift = 0.5*rho*U0**2*St*at*( -(u_list[-1]-u_list[-2])/dx - du_list[-1]/U0 )
    
    """曲げモーメントの二回微分で分布力を計算"""
    dF_list = np.zeros_like(u_list)
    for n in range(N-2):
        Mep1 = Me_list[n+2]
        Me0 = Me_list[n+1]
        Mem1 = Me_list[n]
        dF_list[n+1] = -(Mep1-2.0*Me0+Mem1)/(dx**2)
    
    """尾翼のついている部分は境界条件を入れる"""
    dF_list[N-1] = sigma/mT*( Lift + (Me_list[-1]-Me_list[-2])/dx )
    
    """加速度に変換"""
    utt_list = dF_list/sigma
    """境界条件満たすために調整"""
    utt_list[0] = 0.0
    utt_list[1] = 0.0
    print(t)

    return np.hstack((utt_list,du_list))



"""時間"""
Tmax = 0.6
t_span = [0.0,Tmax]
t_eval = np.linspace(*t_span,100)
"""初期値"""
"""du,u,q,theta"""
u0 = np.zeros_like(x_list)
du0 = np.zeros_like(x_list)
du0[-1] = 1.0
init = np.hstack((du0,u0))
flight = sol.solve_ivp(sim,t_span,init,method='RK45',t_eval=t_eval)

"""結果プロット"""
plt.plot(flight.t,flight.y[2*N-1]) 
plt.xlabel("time[s]")
plt.ylabel("deflection of tailtip [m]",fontsize=15)
plt.grid(visible=True)
plt.savefig("tailtip_tawami.pdf")
plt.clf()

plt.plot(flight.t,flight.y[N-1]) 
plt.xlabel("time[s]")
plt.ylabel("deflectionspeed of tailtip [m]",fontsize=15)
plt.grid(visible=True)
plt.savefig("tailtip_tawamispeed.pdf")
plt.clf()

for n in range(100):
    plt.plot(x_list,flight.y[N:2*N,n]) 
    plt.grid(visible=True)
    plt.savefig("video/tailmovie_"+str(n)+".png")
    plt.clf()

print(flight.y.shape)

"""DMDでのモード抽出"""
X = flight.y[:,0:-1]
Y = flight.y[:,1:]

U2,sig,V = np.linalg.svd(X)
plt.plot(sig)
plt.savefig("sig.pdf")
plt.clf()

r = 6
U = U2[:,:r]
sig = np.diag(sig)[:r,:r]
V = V.conj().T[:,:r]

Atil = np.dot(np.dot(np.dot(U.conj().T, Y), V), np.linalg.inv(sig))
mu,W = np.linalg.eig(Atil)

#固有値を連続システムに変換したものを出力（最も絶対値の小さな複素数の組が第一次モード）
eigc = np.log(mu)/(t_eval[1]-t_eval[0]) 
print(eigc)

#DMDによって得られたモードを出力する
Phi = np.dot(np.dot(np.dot(Y, V), np.linalg.inv(sig)), W)
for n in range(r):
    plt.plot(Phi[:,n])
    plt.savefig("Phi/Phi"+str(n)+".pdf")
    plt.clf()

"""モードのボード線図を書く"""
eig1 = eigc[0]
eig2 = eigc[1]

tf = mt.tf([eig1*eig2],[1,-eig1-eig2,eig1*eig2])
mt.bode(tf)
plt.savefig("tailbode.pdf")
plt.clf()