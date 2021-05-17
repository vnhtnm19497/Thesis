import numpy as np
import matplotlib.pyplot as plt 
import math

class Erlang_1:
    def __init__(self,lamda,b,theta,gama,r,C):
        self.lamda = lamda #интенсивность 1-ого рода
        self.b = b #минимальное требование к ресурсу
        self.theta = theta #размер файла
        self.gama = gama #длина перебаваемого блока данных
        self.r = r  #длина очереди
        self.C = C #число единиц канального ресурса

        self.N = int(C//b) #максимальное число заявок, находящихся на обслуживании ресурсом

    def pr_infinitesimal_generator(self):
        prob = []
        p_01 = sum([(self.lamda * self.theta / self.C) ** n for n in range(1,self.N + 1)])
        p_02 = (self.theta / self.C) ** self.N * sum([self.lamda ** n / math.prod([self.C / self.theta + i * self.gama for i in range(1, n - self.N + 1)]) for n in range(self.N + 1, self.N + self.r + 1)])
        p0 = (1+p_01+p_02)**(-1)
        prob.append(p0)
        return prob

    def coefficient(self):
        a1 = self.pr_infinitesimal_generator()
        UTIL = 1 - a1[0]
        return UTIL
    
    def pr_block(self):
        a2 = self.pr_infinitesimal_generator()
        p0 = a2[0]
        p_Nr = p0*(self.theta/self.C) ** self.N * (self.lamda **(self.N+self.r)) / math.prod([self.C / self.theta + i*self.gama for i in range(1,self.r+1)])
        return p_Nr

    def mean_queue(self):
        a3 = self.pr_infinitesimal_generator()
        p0 = a3[0]
        Lq = p0 * ((self.theta/self.C)**self.N) * sum([i*(self.lamda**(self.N+i))/(math.prod([(self.C/self.theta + j*self.gama) for j in range(1,i+1)])) for i in range(1,self.r+1)])
        return Lq
    
    def mean_service(self):
        a4 = self.pr_infinitesimal_generator()
        p0 = a4[0]
        sum_1 = sum([i*p0*(((self.lamda * self.theta)/self.C)**i) for i in range(self.N+1)])
        sum_2 = self.N * (self.theta/self.C)**self.N * sum([p0*i*(self.lamda**(self.N+i)) / math.prod([self.C / self.theta + j*self.gama for j in range(1,i+1)]) for i in range(1,self.r+1)]) 
        Lser = sum_1 + sum_2
        return Lser

    def mean_system(self):
        Lq = self.mean_queue()
        Lser = self.mean_service()
        Lsys = Lq + Lser 
        return Lsys

    def mean_time_queue(self):
        Lq = self.mean_queue()
        Wq = Lq / self.lamda
        return Wq
    
    def mean_time_service(self):
        Lser = self.mean_service()
        Wser = Lser / self.lamda
        return Wser

    def mean_time_system(self):
        Wq = self.mean_time_queue()
        Wser = self.mean_time_service()
        Wsys = Wq + Wser
        return Wsys
