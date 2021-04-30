import numpy as np
import mpmath as mp
import scipy
import scipy.constants as const

class Schrodinger:

    def __init__(self, vBias,Temp, fermilevel,metalWorkFunc, barrier, noElements = 10):
        self.noElements = noElements
        self.noNodes = noElements+1
        self.elementLength = barrier['length']/self.noElements
        self.nodeXcooords = np.arange(0,self.noNodes,1)*self.elementLength
        # print(self.nodeXcooords)


        self.V = vBias
        self.eV = self.V*const.e
        self.barL = barrier['length']
        self.er = barrier['relative permitivity']
        self.barAffinity = barrier['electron affinity']
        self.T = Temp
        self.Ef = fermilevel
        self.mWF = metalWorkFunc

        # self.barEVArr = np.linspace((self.Ef+(self.mWF-self.barAffinity)),(self.Ef-self.eV+self.mWF-self.barAffinity),self.noNodes)
        # print(self.barEVArr)
        x= [0,self.barL]
        y = [(self.Ef+self.mWF-self.barAffinity),(self.Ef-self.eV+self.mWF-self.barAffinity)]
        self.barEvEq = np.poly1d(np.polyfit(x, y, 1))
       
        

        #constants
        self.e0 = const.physical_constants['vacuum electric permittivity'][0]
        self.lam0 = np.sqrt((const.hbar**2)/(2*const.e*const.m_e*self.V))
        self.lam = self.barL/self.lam0
        self.phig = self.eV/const.physical_constants['Hartree energy in eV'][0]
        self.Tbar = const.Boltzmann*self.T/self.eV

        self.gammaConst = 9/(4*const.pi) * (self.lam**2)/np.sqrt(2*self.phig)*self.Tbar

        # self.phix=self.phiX()
        print(self.probability(1e-9)        )

    def gamma1(self):
        func= lambda n : mp.log(1+mp.exp((-n/self.eV)/self.Tbar))*self.probability(n/self.eV)
        integral = mp.quad(func,[-mp.inf,mp.inf])
        return self.gammaConst*integral

    def gamma2(self):
        func= lambda n : mp.log(1+mp.exp(((-n/self.eV)+1)/self.Tbar))*self.probability(n/self.eV)
        integral = mp.quad(func,[-mp.inf,mp.inf])
        return self.gammaConst*integral

    def gammaNet(self):
        return self.gamma1()-self.gamma2()

    def probability(self,x):
        root1 = mp.findroot(lambda n : x-self.phiX(n), -10000, verbose=True)
        # print(root1)
        root2 = mp.findroot(lambda n : x-self.phiX(n), 10000, verbose=True)
        # print(root2)
        integral = mp.quad(lambda n: mp.sqrt(2*const.m_e*( self.phiX(n)-x ) ), [root1,root2])
        return mp.exp(-2/const.hbar * integral)

    def eVarr(self):
        pass

    def phiX(self,x):
        # phi = np.zeros(self.noNodes)
        phi = self.Ef + self.mWF -self.barAffinity +self.phiImg(x)+self.barEvEq(x) +self.phixc(x) #TODO: impliment PhiXC
        return phi
    
    # def phiImg(self):
    #     phi = np.zeros(self.noNodes)
    #     for index, x in enumerate(self.nodeXcooords):

    #         summationeq = lambda n: (n*self.barL/(n**2*self.barL**2 - x**2) - 1/n*self.barL )

    #         phi[index] = ((-const.e**2)/(8*const.pi*self.er*self.e0)) * (0.5*x + mp.nsum(summationeq , [1,mp.inf], ignore=True) )
    #     return phi

    def phiImg(self, x):
        summationeq = lambda n: (n*self.barL/(n**2*self.barL**2 - x**2) - 1/n*self.barL )
        return ((-const.e**2)/(8*const.pi*self.er*self.e0)) * (0.5*x + mp.nsum(summationeq , [1,mp.inf], ignore=True) )


    def phiXc(self,x):
        ex = (-3/4)*(3/2*const.pi)**(2/3)*(1/self.seitzR(x))
        return

    def seitzR(self,x):
        a0 = const.physical_constants['Bohr radius'][0]

        # rs = (3/(4*const.pi*self.eDensity(x)))**3

        # (rs *a0)**3 = 3/(4*const.pi*self.eDensity(x))
        rs = mp.cbrt(3/(4*const.pi*self.eDensity(x)))/a0
        return rs

    def eDensity(self,x):
        # return 1e-20
        pass

if __name__=='__main__':

    barrier = {
        'relative permitivity': 1,#9.7,
        'effective mass': const.m_e,#0.5*const.m_e,
        'electron affinity': 0,#3*const.e,
        'length': 5e-9,
    }

    Schrodinger(vBias=10, metalWorkFunc=5*const.e, Temp=10, barrier=barrier, fermilevel=5*const.e)




