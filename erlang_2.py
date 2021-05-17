import numpy as np
import matplotlib.pyplot as plt
import math

class Erlang_2:
    def __init__(self,M,eps,b,theta,gama,r,C):
        self.M = M #количество пользователей
        self.eps = eps #интенсивнось 2-ого рода
        self.b = b #минимальное требование к ресурсу
        self.theta = theta #размер файла
        self.gama = gama #длина перебаваемого блока данных
        self.r = r  #длина очереди
        self.C = C #число единиц канального ресурса

        self.N = int(C//b) #максимальное число заявок, находящихся на обслуживании ресурсом

    def pr_infinitesimal_generator(self):
        prob = []
        if (self.M >= self.N + self.r):
            p01_1 = sum([math.perm(self.M, n)*(((self.eps*self.theta)/self.C)**n) for n in range(1, self.N+1)])
            p01_2 = sum([math.perm(self.M,n)*(self.eps**n)/(math.prod([(self.C/self.theta + i*self.gama) for i in range(1,n-self.N+1)])) for n in range(self.N+1,self.N+self.r+1)])
            p01 = (1+p01_1+p01_2)**(-1)
            prob.append(p01)
        elif (self.M >= self.N and self.M < self.N+self.r):
            p02_1 = sum([math.perm(self.M, n)*(((self.eps*self.theta)/self.C)**n) for n in range(1, self.N+1)])
            p02_2 = sum([math.perm(self.M,n)*(self.eps**n)/(math.prod([(self.C/self.theta + i*self.gama) for i in range(1,n-self.N+1)])) for n in range(self.N+1,self.M+1)])
            p02 = (1+p02_1+p02_2)**(-1)
            prob.append(p02)
        else:
            p03_1 = sum([math.perm(self.M, n)*(((self.eps*self.theta)/self.C)**n) for n in range(1, self.M)])
            p03 = (1+p03_1)**(-1) 
            prob.append(p03)
        return prob
    
    def coefficient(self):
        a1 = self.pr_infinitesimal_generator()
        p0 = a1[0]
        UTIL = 1 - p0
        return UTIL
    
    def pr_block(self):
        a2 = self.pr_infinitesimal_generator()
        p0 = a2[0]
        if (self.M >= self.N + self.r):
           prod_1 = ((self.theta/self.C)**(self.N))*(math.perm(self.M, self.N+self.r))
           prod_2 = ((self.eps)**(self.N+self.r))/(math.prod([self.C/self.theta + i*self.gama for i in range(1,self.r+1)]))
           p_block = prod_1 * prod_2 * p0
           return p_block
        else:
            return 0
    
    def mean_service(self):
        a3 = self.pr_infinitesimal_generator()
        p0 = a3[0]
        if (self.M >= self.N + self.r):
            sum1_1 = sum([p0*i*math.perm(self.M,i)*((self.eps*self.theta/self.C)**i) for i in range(self.N+1)])
            sum1_2 = self.N*sum([p0*math.perm(self.M,i)*((self.theta/self.C)**self.N)*((self.eps)**i/math.prod([self.C/self.theta + j*self.gama for j in range(1,i-self.N+1)])) for i in range(self.N+1,self.N+self.r+1)])
            Lser_1 = sum1_1 + sum1_2
            return Lser_1
        elif (self.M >= self.N and self.M < self.N + self.r):
            sum2_1 = sum([p0*i*((self.eps*self.theta/self.C)**i)*math.perm(self.M,i) for i in range(self.N+1)])
            sum2_2 = self.N*sum([p0*math.perm(self.M,i)*((self.theta/self.C)**self.N)*((self.eps)**i/math.prod([self.C/self.theta + j*self.gama for j in range(1,i-self.N+1)])) for i in range(self.N+1,self.M+1)])
            Lser_2 = sum2_1 + sum2_2
            return Lser_2
        else:
            Lser_3 = sum([p0*i*((self.eps*self.theta/self.C)**i)*math.perm(self.M,i) for i in range(self.M+1)])
            return Lser_3

    def mean_system(self):
        a4 = self.pr_infinitesimal_generator()
        p0 = a4[0]
        if (self.M >= self.N + self.r):
            sum1_1 = sum([p0*i*math.perm(self.M,i)*(self.eps*self.theta/self.C)**i for i in range(self.N+1)])
            sum1_2 = sum([p0*i*math.perm(self.M,i)*((self.theta/self.C)**self.N)*((self.eps**i)/(math.prod([self.C/self.theta + j*self.gama for j in range(1,i-self.N+1)]))) for i in range(self.N+1,self.N+self.r+1)])
            Lsys_1 = sum1_1 + sum1_2
            return Lsys_1
        elif (self.M >= self.N and self.M < self.N+self.r):
           sum2_1 = sum([p0*i*math.perm(self.M,i)*(self.eps*self.theta/self.C)**i for i in range(self.N+1)])
           sum2_2 = sum([p0*i*math.perm(self.M,i)*((self.theta/self.C)**self.N)*((self.eps**i)/(math.prod([self.C/self.theta + j*self.gama for j in range(1,i-self.N+1)]))) for i in range(self.N+1,self.M+1)])
           Lsys_2 = sum2_1 + sum2_2
           return Lsys_2
        else:
            Lsys_3 = sum([p0*i*((self.eps*self.theta/self.C)**i)*math.perm(self.M,i) for i in range(self.M+1)])
            return Lsys_3
        
    def mean_queue(self):
        Lsys = self.mean_system()
        Lser = self.mean_service()
        Lq = Lsys - Lser
        return Lq

    def mean_time_service(self):
        Lser = self.mean_service()
        a5 = self.pr_infinitesimal_generator()
        p0 = a5[0]
        if (self.M >= self.N + self.r) and (self.M >= self.N and self.M < self.N+self.r):
            sum1 = sum([p0*(self.M-n)*self.eps*math.perm(self.M,n)*((self.eps*self.theta/self.C)**n) for n in range(self.N+1)])
            sum2 = sum([p0*(self.M-n)*self.eps*math.perm(self.M,n)*((self.theta/self.C)**self.N)*((self.eps**n)/math.prod([self.C/self.theta + i*self.gama for i in range(1,n-self.N+1)]))for n in range(self.N+1, self.M)])
            dominator1 = sum1 + sum2
            Wser1 = Lser/dominator1
            return Wser1
        else:
            dominator2 = sum([p0*(self.M-n)*self.eps*math.perm(self.M,n)*((self.eps*self.theta/self.C)**n) for n in range(self.M)])
            Wser2 = Lser/dominator2
            return Wser2

    def mean_time_system(self):
        Lsys = self.mean_system()
        a5 = self.pr_infinitesimal_generator()
        p0 = a5[0]
        if (self.M >= self.N + self.r) and (self.M >= self.N and self.M < self.N+self.r):
            sum1 = sum([p0*(self.M-n)*self.eps*math.perm(self.M,n)*((self.eps*self.theta/self.C)**n) for n in range(self.N+1)])
            sum2 = sum([p0*(self.M-n)*self.eps*math.perm(self.M,n)*((self.theta/self.C)**self.N)*((self.eps**n)/math.prod([self.C/self.theta + i*self.gama for i in range(1,n-self.N+1)]))for n in range(self.N+1, self.M)])
            dominator1 = sum1 + sum2
            Wsys1 = Lsys/dominator1
            return Wsys1
        else:
            dominator2 = sum([p0*(self.M-n)*self.eps*math.perm(self.M,n)*((self.eps*self.theta/self.C)**n) for n in range(self.M)])
            Wsys2 = Lsys/dominator2
            return Wsys2

    def mean_time_queue(self):
        Wsys = self.mean_time_system()
        Wser = self.mean_time_service()
        Wq = Wsys - Wser
        return Wq


            
            




        
            