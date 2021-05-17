import numpy as np
import matplotlib.pyplot as plt
import math

#Исходные данных
N1 = [0,1,2]
N2 = [0,1,2]
M1 = [0,1]
M2 = [0,1]
N = max(N1)+max(N2)
M = max(M1)+max(M2)
K = 2
lam1 = 1
lam2 = 5
mu1 = 0.1
mu2 = 0.1
delta = 0.02

def state_space(N1,N2,M1,M2):
    X = []
    for N_1 in N1:
        for N_2 in N2:
            for n1 in range(len(N1)):
                for n2 in range(len(N2)):
                    if  n1 <= N_1 and n2 <= N_2 and n1+n2 <= N and N_1+N_2 == max(N1):
                        X.append([n1,n2,0,0,N_1,N_2])

    for N_1 in N1:
        for N_2 in N2:
            for n1 in range(len(N1)):
                for n2 in range(len(N2)):
                        for M_1 in M1:
                            for M_2 in M2:
                                for m1 in range(len(M1)):
                                    for m2 in range(len(M2)):
#Первая очередь пуста и вторая очередь не пуста
                                        if 0<=n1<N_1 and n2 == N_2 and 0<m2<=M_2 and N_1 + N_2 == max(N1) and n1+n2<=max(N1):
                                            X.append([n1+1,n2,0,m2,N_1,N_2])
                                        elif 0<n1<=N_1 and n2 == N_2 and 0<m2<=M_2 and N_1 + N_2 == max(N1):
                                            X.append([n1-1,n2,0,m2,N_1,N_2])
                                        elif 0<=n1<=N_1 and n2 == N_2 and 0 <m2<=M_2 and N_1+N_2 == max(N1) and n1 + n2 <= max(N1):
                                            X.append([n1,n2,0,m2,N_1,N_2])
                                        elif 0<=n1<= N_1 and 0<=n2<=N_2 and 0 < m2 <= M_2 and M_2 <= N_1 - n1 and N_1 + N_2 == max(N1):
                                            X.append([n1,n2,0,m2,N_1+1,N_2-1])
#Первая очередь не пуста и вторая очередь пуста                                           
                                        elif n1 == N_1 and 0<=n2<N_2 and 0<m1<=M_1 and N_1+N_2==max(N1):
                                            X.append([n1,n2+1,m1,0,N_1,N_2])
                                        elif n1 == N_1 and 0<n2<=N_2 and 0<m1<=M_1 and N_1+N_2==max(N1):
                                            X.append([n1,n2-1,m1,0,N_1,N_2])
                                        elif n1==N_1 and 0<=n2<=N_2 and 0<m1<=M_1 and N_1+N_2 == max(N1):
                                            X.append([n1,n2,m1,0,N_1,N_2])
                                        elif 0<=n1<=N_1 and 0<=n2<=N_2 and 0<m1<=M_1 and M_1<=N_2-n2 and N_1+N_2 == max(N1):
                                            X.append([n1,n2,m1,0,N_1-1,N_2+1])
#Оба очередь не пуста
                                        elif n1 == N_1 and n2 == N_2 and m1+m2==M and N_1 + N_2 == max(N1) and n1+n2==max(N1): 
                                            X.append([n1,n2,m1,m2,N_1,N_2])
    X_list = []                      
    for i in X:
        if i not in X_list:
            X_list.append(i)
    return X_list

def infinitesimal_generator(lam1,lam2,mu1,mu2,delta):
    a = np.array([[-(lam1+lam2),lam2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,lam1,0.,0.,0.,0.,0.,0.,0.,0.,0.],
    [mu2,-(lam1+lam2+mu2),0.,0.,0.,0.,0.,0.,0.,0.,lam1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
    [0.,mu2,-(lam1+lam2+mu2),0.,0.,0.,0.,0.,0.,0.,0.,lam1,0.,lam2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
    [0.,0.,0.,-(lam1+lam2),lam2,lam1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
    [0.,0.,0.,mu2,-(lam1+lam2+mu2),0.,lam1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,lam2,0.,0.,0.,0.,0.],
    [0.,0.,0.,mu1,0.,-(lam1+lam2+mu1),lam2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
    [0.,0.,0.,0.,mu1,mu2,-(lam1+lam2+mu1+mu2),0.,0.,0.,0.,0.,0.,0.,0.,0.,lam2,lam1,0.,0.,0.,0.,0.,0.,0.],
    [0.,0.,0.,0.,0.,0.,0.,-(lam1+lam2),lam1,0.,0.,0.,0.,0.,lam2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
    [0.,0.,0.,0.,0.,0.,0.,mu1,-(lam1+lam2+mu1),lam1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,lam2,0.,0.,0.],
    [0.,0.,0.,0.,0.,0.,0.,0.,mu1,-(lam1+lam2+mu1),0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,lam2,0.,lam1],
    [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-(lam2+mu2+delta),lam2,0.,0.,0.,mu2,0.,0.,0.,0.,0.,0.,0.,0.,0.],
    [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,mu2,-(lam2+mu2),lam2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
    [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,mu2,-mu2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
    [0.,0.,mu2,0.,0.,0.,0.,0.,0.,0.,0.,0.,lam1,-(lam1+mu2),0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
    [0.,0.,0.,0.,delta,0.,0.,0.,0.,0.,0.,0.,0.,0.,-(lam1+delta),0.,0.,0.,0.,0.,0.,lam1,0.,0.,0.],
    [0.,0.,0.,0.,0.,delta,0.,0.,0.,0.,lam2,0.,0.,0.,0.,-(lam2+delta),0.,0.,0.,0.,0.,0.,0.,0.,0.],
    [0.,0.,0.,0.,0.,0.,mu2,0.,0.,0.,0.,0.,0.,0.,0.,0.,-(lam1+mu1+mu2),0.,lam1,mu1,0.,0.,0.,0.,0.],
    [0.,0.,0.,0.,0.,0.,mu1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-(lam2+mu1+mu2),lam2,0.,mu2,0.,0.,0.,0.],
    [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,mu1,mu2,-(mu1+mu2),0.,0.,0.,0.,0.,0.],
    [0.,0.,delta,0.,mu2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,lam1,0.,0.,-(lam1+mu2+delta),0.,0.,0.,0.,0.],
    [0.,0.,0.,0.,0.,mu1,0.,0.,0.,delta,0.,0.,0.,0.,0.,0.,0.,lam2,0.,0.,-(lam2+mu1+delta),0.,0.,0.,0.],
    [0.,0.,0.,0.,0.,0.,delta,0.,0.,0.,0.,0.,0.,0.,mu1,0.,0.,0.,0.,0.,0.,-(lam1+mu1+delta),lam1,0.,0.],
    [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,mu1,-(lam1+mu1),lam1,0.],
    [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,mu1,-mu1,0.],
    [0.,0.,0.,0.,0.,0.,0.,0.,0.,mu1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,lam2,-(lam2+mu1)]
    ])
    return a

def pr_infinitesimal_generator(lam1,lam2,mu1,mu2,delta):
    a = infinitesimal_generator(lam1,lam2,mu1,mu2,delta)
    n = len(state_space(N1,N2,M1,M2))
    b = np.zeros(n)
    a[:,0] = np.ones(n)
    b[0] = 1.
    a = a.T
    p = np.linalg.solve(a, b)
    return p

def UTIL(lam1,lam2,mu1,mu2,delta):
    p = pr_infinitesimal_generator(lam1,lam2,mu1,mu2,delta)
    L1 = (p[5]+p[6]+p[8]+p[16]+p[17]+p[18]+p[20])+2*(p[9]+p[22]+p[23]+p[24])
    L2 = (p[1]+p[4]+p[6]+p[10]+p[16]+p[17]+p[18]+p[19])+2*(p[2]+p[11]+p[12]+p[13])
    UTIL = (L1+L2)/N
    return UTIL

def alpha(lam1,lam2,mu1,mu2,delta):
    p = pr_infinitesimal_generator(lam1,lam2,mu1,mu2,delta)
    a1 = (p[7]+p[8]+p[10]+p[14]+p[21]+p[22]+p[23]+p[24]) + (N1[1]/N)*(p[3]+p[4]+p[5]+p[6]+p[16]+p[17]+p[18]+p[19]+p[20])
    a2 = (p[0]+p[1]+p[2]+p[10]+p[11]+p[12]+p[13]+p[15]) + (N1[1]/N)*(p[3]+p[4]+p[5]+p[6]+p[16]+p[17]+p[18]+p[19]+p[20])
    alpha = (a1+a2)/K
    return alpha

def beta(lam1,lam2,mu1,mu2,delta):
    p = pr_infinitesimal_generator(lam1,lam2,mu1,mu2,delta)
    beta = p[10] + p[14] + p[15] + p[19] + p[20] + p[21]
    return beta  

delta1 = np.arange(0.01,10,0.01)
UTIL = [UTIL(lam1, lam2, mu1,mu2,delta) for delta in delta1]
alpha = [alpha(lam1,lam2,mu1,mu2,delta) for delta in delta1]
beta = [beta(lam1,lam2,mu1,mu2,delta) for delta in delta1]
plt.figure(figsize=(10, 7), dpi=80)
plt.xticks(ticks=np.arange(0.01,11,1))
plt.yticks(ticks=np.arange(0.0004,0.6004,0.1))
plt.ylabel('Коэффициент эффективности модели')
plt.xlabel('Интервал нарезки $\Delta$, s')
plt.title('График зависимости коэффициентов эффективности модели от интервалов нарезки')
plt.plot(delta1, UTIL,'r-', label='UTIL')
plt.plot(delta1, alpha,'b-', label=r'$\alpha$')
plt.plot(delta1,beta,'y-',label=r'$\beta$')
plt.legend(loc='best')
plt.grid(color='black', ls='--')
plt.show()

