import numpy as np


class Simaan2008():
    def __init__(self):
        self.start_t = 0 
        self.step   = 0.0001
        self.end_t   = 3
            
            
        self.T = list(np.arange(self.start_t,self.end_t,self.step)) #
        self.n = len(self.T)

        #Cardiovascular System
        self.HR   = 60
        self.Emax = 2.5 #elastance amplitude
        self.Emin = 0.06
        self.tc   = 60/self.HR #cardiac cycle

        self.t_max = 0.2 + 0.15*self.tc; #max time cardiac cycle
        self.E = self.elastance(self.T) #normalized elastance


            
        self.Rs  = 1.0000 #vascular resistence
        self.Rm  = 0.0050 #mitral resistence
        self.Cae = 4.4000 #artic elastance
        self.Ra  = 0.0010 #aortic resistence
        self.Rc  = 0.0398 
        self.Cs  = 1.3300 #Elastance system
        self.Cao = 0.0800 #Elastance aortic
        self.Ls  = 0.0005 

        self.Vo = 10 #Initial Pressure

        self.Pao = np.zeros_like(self.T) #aortic pressure
        self.Qa  = np.zeros_like(self.T) #aortic flow
        self.Vve = np.zeros_like(self.T) #left ventricle
        self.Pas = np.zeros_like(self.T) #aortic pressure system
        self.Pae = np.zeros_like(self.T) #left artic pressure
        self.Pve = np.zeros_like(self.T) #left ventricle pressure
        self.Dm_ = np.zeros_like(self.T) 
        self.Da_ = np.zeros_like(self.T) 


        #Initial Conditions
        self.Pao[0] =  90
        self.Qa[0]  =   0
        self.Vve[0] = 140 #Left ventricular volume
        self.Pas[0] =  90
        self.Pae[0] =  10

        self.Pve[0] = self.E[0]* (self.Vve[0] - self.Vo) #Left ventricular pressure

        self.x = np.transpose([self.Pao[0], self.Qa[0], self.Vve[0], self.Pas[0], self.Pae[0]])

        #Initial States of diodes
        self.Dm = 0
        self.Da = 0



        #Parameters
    def parameters(self,start_t, end_t, step, HRs, Emaxs):
 
        self.start_t = start_t 
        self.step   = step
        self.end_t   = end_t
            

        #time scale used
        self.T = list(np.arange(self.start_t,self.end_t,self.step)) #
        self.n = len(self.T)

        #Cardiovascular System
        self.HR   = HRs
        self.Emax = Emaxs #elastance amplitude
        self.Emin = 0.06
        self.tc   = 60/self.HR #cardiac cycle

        self.t_max = 0.2 + 0.15*self.tc; #max time cardiac cycle
        self.E = self.elastance(self.T) #normalized elastance


            
        self.Rs  = 1.0000 #vascular resistence
        self.Rm  = 0.0050 #mitral resistence
        self.Cae = 4.4000 #artic elastance
        self.Ra  = 0.0010 #aortic resistence
        self.Rc  = 0.0398 
        self.Cs  = 1.3300 #Elastance system
        self.Cao = 0.0800 #Elastance aortic
        self.Ls  = 0.0005 

        self.Vo = 10 #Initial Pressure

        self.Pao = np.zeros_like(self.T) #aortic pressure
        self.Qa  = np.zeros_like(self.T) #aortic flow
        self.Vve = np.zeros_like(self.T) #left ventricle
        self.Pas = np.zeros_like(self.T) #aortic pressure system
        self.Pae = np.zeros_like(self.T) #left artic pressure
        self.Pve = np.zeros_like(self.T) #left ventricle pressure
        self.Dm_ = np.zeros_like(self.T) 
        self.Da_ = np.zeros_like(self.T) 


        #Initial Conditions
        self.Pao[0] =  90
        self.Qa[0]  =   0
        self.Vve[0] = 140 #Left ventricular volume
        self.Pas[0] =  90
        self.Pae[0] =  10

        self.Pve[0] = self.E[0]* (self.Vve[0] - self.Vo) #Left ventricular pressure

        self.x = np.transpose([self.Pao[0], self.Qa[0], self.Vve[0], self.Pas[0], self.Pae[0]])

        #Initial States of diodes
        self.Dm = 0
        self.Da = 0
        return

    def elastance(self, t):

        tn = np.asarray(t)%self.tc/self.t_max
        En = 1.55 * np.power(np.asarray(tn)/.7, 1.9) / (1 + np.power(np.asarray(tn)/.7, 1.9)) \
        / (1 + np.power(np.asarray(tn)/1.17, 21.9))
        return (self.Emax-self.Emin)*En + self.Emin
    
    def rugkut4(self,step,A,x,B,i):
        
        Ar = np.array(A)
        xr = np.array(x)
        Br = np.array(B)

        xdot = np.matmul(Ar, xr) + Br
        kx1 = step*xdot

        x1 = x + 0.5*kx1
        xdot = np.matmul(Ar, x1) + Br
        kx2 = step*xdot

        x1 = x + 0.5*kx2
        xdot = np.matmul(Ar, x1) + Br
        kx3 = step*xdot

        x1 = x + kx3
        xdot = np.matmul(Ar, x1) + Br
        kx4 = step*xdot

        value = np.asarray(x + (kx1 + 2*kx2 + 2*kx3 + kx4)/6)

        return value
