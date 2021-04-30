import numpy as np
import mpmath as mp
from CurvatureCalc import ComputeCurvature 


class FDMCantilever:
    def __init__(self, cantileverSettings, noElements={'Number of Elements':100},geomertyMods={}, coatings={},coatingMods={},mechanical={}):

        self.cantileverSettings = cantileverSettings
        self.cantilever = {
            'w'  : cantileverSettings["Width"],
            'l'  : cantileverSettings["Length"],
            't'  : cantileverSettings["Thickness"],
            'E'  : cantileverSettings["Youngs Modulus"],
            'vs' : cantileverSettings["Possion's Ratio"],
            'angle':[cantileverSettings["angle cantilever"], cantileverSettings["angle degrees const"],cantileverSettings["angle rad"]],
            'angle change': [cantileverSettings["angle change"],cantileverSettings["changing angle function"]],
            'Part Angle': cantileverSettings["Part Angle Cantilever"]
        }
        self.mechanical  = {"Force":0,"Offset Force":False}if not mechanical else mechanical


        self.noElements = noElements['Number of Elements']
        self.noNodes = self.noElements+1
        self.elementLength = self.cantilever["l"]/self.noElements

        #setup node arrays
        self.nodeXCoords = np.arange(0,self.noNodes,1)*self.elementLength #ne_x? #TODO: variable mesh
        self.nodeCentroidCoords = self.nodeXCoords-(self.elementLength/2) #ne_xc #TODO: variable mesh
        self.nodeCentroidCoords[0] = 0

        self.cantilever['angle change'][1] = lambda: cantileverSettings["changing angle function"](x=self.nodeXCoords)
        
        self.nodeAlongCantileverCoords = self.nodeXCoords if not cantileverSettings["Part Angle Cantilever"]["Part Angle"] else np.zeros(self.noNodes)   #n_A

        self.baseElementXLenArr = np.ones(self.noNodes)*self.elementLength    #Base Element X-Length Array #Le_xa
        self.baseElementXLenArr[0] = 0
        self.baseCantileverLen = self.cantilever['l']

        #pristmatic widths
        self.widthArr = np.ones(self.noNodes)*self.cantilever['w']
        
        #TODO: Variable mesh

        #Cantilever cutouts /tapers
        if geomertyMods:
            self.modify(modDict=geomertyMods, arr=self.widthArr)

        
        self.element2ndMomentOfArea = self.widthArr*self.cantilever['t']**3/12
        self.element2ndMomentOfArea[0] = np.Infinity

        #Coating
        if coatings:
            if coatings['Coating']:
                self.coatings = coatings
                self.coatingWidthArr=np.asarray([self.coatings['width'] if round(self.nodeCentroidCoords[i],8)>self.coatings['startX'] and round(self.nodeCentroidCoords[i],8)<=self.coatings['endX'] else 0 for i in range(0,self.noNodes,1)])   
        
                if coatingMods:
                    self.modify(modDict=coatingMods, arr=self.coatingWidthArr)

            #TODO: Composite Mechanical <-- coating changes + underCoating goes in here
            self.compositeMechanicalCalc()


        self.partAngledAngleArr = np.zeros(self.noNodes)    #PAr
        self.angleCantilever()

        # self.mechanicalInputCalc()
        # self.mechanicalOutCalc()

    def doMechanical(self):
        self.mechanicalInputCalc()
        # self.mechanicalOutCalc()
        return self.mechanicalOutCalc()
        

    def modify(self, modDict, arr):
        if 'Taper' in modDict.keys() :
            if  modDict['Taper']['Taper']:
                self.addTaper(arr=arr, params=modDict['Taper']['parameters'])
        
        if 'Gaps/Width' in modDict.keys() :
            if modDict['Gaps/Width']["Gap/Width Change"]:
                params = modDict['Gaps/Width']
                params = {x: modDict['Gaps/Width'][x] for x in modDict['Gaps/Width'] if x not in {'Gap/Width Change'}} #eliminate True from dict 

                for key in params.keys():
                    width = params[key]
                    if not width['active']: 
                        continue 
                    self.addWidthChange(arr=arr,params=width)

        # if 'coating change' in modDict.keys() :   #TODO: IMPLEMENT coating changes
        #     if modDict['coating change']["coating change"]:
        #         params = modDict['Gaps/Width']
        #         params.pop('Gap/Width Change')

        #         for key in params.keys():
        #             width = params[key]
        #             if not width['active']: 
        #                 continue 
        #             self.addWidthChange(arr=arr,params=width)

            


    def addTaper(self, arr, params):
        TM=(self.cantilever['w']-params['end width'])/(params['start']-params['end'])
        TS_N = next(index-1 for index,val in enumerate(self.nodeCentroidCoords) if val>=params['start'])

        for i in range(TS_N,self.noNodes,1):

            if self.nodeCentroidCoords[i]>params['start'] and self.nodeCentroidCoords[i]<=params['end']:
                arr[i]=(TM * (self.nodeCentroidCoords[i]-params['start']) ) + self.cantilever['w']
                
            elif self.nodeCentroidCoords[i]>params['end']:
                arr[i]=params['end width']

    def addWidthChange(self, arr, params):
        GM=(params['end change']-params['start change'])/(params['endX']-params['startX'])
        GS_N=next(index-1 for index,val in enumerate(self.nodeCentroidCoords) if val>=params['startX'])
        GE_N=next((index-1 for index,val in enumerate(self.nodeCentroidCoords) if val>=params['endX']),self.noNodes)

        for i in range (GS_N,GE_N,1):
            if self.nodeCentroidCoords[i]>params['startX'] and self.nodeCentroidCoords[i]<=params['endX']:
                arr[i]=arr[i]-params['start change']-(GM*(self.nodeCentroidCoords[i]-params['startX']))

    def angleCantilever(self):
        self.angledArr = self.cantilever['angle change'][1]() if self.cantilever['angle change'][0] else np.zeros(self.noNodes)
        
        if self.cantilever['Part Angle']['Part Angle']:
            PE = self.cantilever['Part Angle']['end']
            PS = self.cantilever['Part Angle']['start']
            PAd = self.cantilever['Part Angle']['degrees'] 

            self.baseCantileverLen = (self.cantilever['l'] - (PE-PS)) + ((PE-PS)/np.cos(np.deg2rad(PAd))) #Total Cantilever Along-Length (m)
            self.partAngledEndAlongLength = PS+((PE-PS)/np.cos(np.deg2rad(PAd)))  #Part Angled End Along-Length (m)
            #Part Angled Substrate Element X-Length (m)
            PS_N = next(i for i in range(1,self.noNodes,1) if self.nodeCentroidCoords[i-1]<=PS and self.nodeCentroidCoords[i]>PS)
            PE_N = next(i if self.nodeCentroidCoords[i-1]<=PE and self.nodeCentroidCoords[i]>PE else self.noNodes for i in range(PS_N,self.noNodes,1))
            Pe_x=(self.partAngledEndAlongLength-self.nodeXCoords[PS_N-1])/(PE_N-PS_N)
            #Part Angled Element X-Length Array
            self.baseElementXLenArr = np.asarray([Pe_x if self.nodeCentroidCoords[i]>PS and self.nodeCentroidCoords[i]<=PE else self.baseElementXLenArr[i] for i in range(0,self.noNodes,1)])

            self.partAngledAngleArr = np.asarray([np.deg2rad(PAd) if self.nodeCentroidCoords[i]>=PS and self.nodeCentroidCoords[i]<=PE else 0 for i in range(0,self.noNodes,1)]) 

            for i in range(1,self.noNodes,1):
                self.nodeAlongCantileverCoords[i] =  self.nodeAlongCantileverCoords[i-1]+self.baseElementXLenArr[i]

    def mechanicalInputCalc(self):
        if self.mechanical["Offset Force"]:
            PS = self.cantilever['Part Angle']['start']
            PE = self.cantilever['Part Angle']['end']
            OF_Px=self.mechanical['Force application X point']
            
            OF_PA = ((PS+((OF_Px-PS)/np.cos(self.partAngledAngleArr))) if self.cantilever['Part Angle']['Part Angle'] and OF_Px>PS and OF_Px<=PE else OF_Px)
            self.offsetForceZeroArr=[0 for i in range(self.noNodes) if self.nodeXCoords[i]>OF_Px] 

        F_Px=(OF_Px if self.mechanical["Offset Force"] else self.cantilever['l'])
        #Force Application Along-Point
        F_PA=(OF_PA if self.mechanical["Offset Force"] else self.baseCantileverLen)
        #Force Application Point Array
        self.forceAppPointArr=np.asarray([F_PA if self.cantilever['Part Angle']['Part Angle'] and self.nodeCentroidCoords[i]>self.cantilever['Part Angle']['start'] and self.nodeCentroidCoords[i]<=self.cantilever['Part Angle']['end'] else F_Px for i in range(0,self.noNodes)])
        

    def mechanicalOutCalc(self):

        angleConst = np.cos(self.partAngledAngleArr+self.cantilever['angle'][2]+self.angledArr)**2
        #Spring Constants
        #Spring Constant component of Y displacement
        springUY = (self.baseElementXLenArr**2/(self.cantilever['E']*self.element2ndMomentOfArea)*((self.baseElementXLenArr/3)  \
                            +((self.forceAppPointArr-self.nodeAlongCantileverCoords)/2)))                                       \
                            *angleConst
       
        #Spring Constant Component for Rotational Y Displacement
        springRXY = (self.baseElementXLenArr**2/(self.cantilever['E']*self.element2ndMomentOfArea)*((self.baseElementXLenArr/2)  \
                            +self.forceAppPointArr-self.nodeAlongCantileverCoords))                                       \
                            *angleConst
       
        if self.mechanical['Offset Force']:
            springUY[self.noNodes-len(self.offsetForceZeroArr):self.noNodes]=self.offsetForceZeroArr
            springRXY[self.noNodes-len(self.offsetForceZeroArr):self.noNodes]=self.offsetForceZeroArr

        #Total Spring Constant
        springRUY=np.zeros(self.noNodes)
        for i in range(1,self.noNodes):
            # if VM == "Y" and Le_xa[i] != SM_Le_x: #TODO: variable mesh
            #     k_R_UY[i]=(np.sum(k_RXY[0:VM_N])/SM_Le_x*RM_Le_x)+np.sum(k_RXY[VM_N:i])
            # else:
            springRUY[i]=np.sum(springRXY[0:i]) 

        springT=springUY+springRUY
        springT[0]= np.Infinity
        springE = 1/springT
        springE[0] = np.Infinity

        springAlong=np.zeros(self.noNodes)
        for i in range(1,self.noNodes+1,1):
            ke = 1/springE
            ke[i:self.noNodes]=0
            springAlong[0]=0
            if np.sum(ke) == 0:                                                 
                springAlong[i-1]=0                                                     #Remove division by zero error
            else:
                springAlong[i-1]=1/np.sum(ke)  
        
        
        #Displacements
        F = self.mechanical['Force']
        #Y displacement
        dispY = ((F*self.baseElementXLenArr**3/(3*self.cantilever['E']*self.element2ndMomentOfArea))            \
                    +(F*(self.forceAppPointArr-self.nodeAlongCantileverCoords)*self.baseElementXLenArr**2       \
                    /(2*self.cantilever['E']*self.element2ndMomentOfArea)))                                      \
                    *angleConst**2
        #Rotational displacement
        dispR=((F*self.baseElementXLenArr**2/(2*self.cantilever['E']*self.element2ndMomentOfArea))                  \
                    +(F*(self.forceAppPointArr-self.nodeAlongCantileverCoords)*self.baseElementXLenArr           \
                    /(self.cantilever['E']*self.element2ndMomentOfArea)))                                       \
                    *angleConst**2

        if self.mechanical['Offset Force']:
            dispY[self.noNodes-len(self.offsetForceZeroArr):self.noNodes]=self.offsetForceZeroArr
            dispR[self.noNodes-len(self.offsetForceZeroArr):self.noNodes]=self.offsetForceZeroArr

        
        #rotational Y displacement
        dispRUY = np.zeros(self.noNodes)
        dispUY= np.zeros(self.noNodes)
        dispUY= np.tan(dispR)*self.baseElementXLenArr


        # dispUY = np.array(dispUY,dtype=float)
        for i in range(1,self.noNodes):
            # if VM == "Y" and Le_xa[i] != SM_Le_x: #TODO: variable mesh
            #     R_UY[i]=(np.sum(UY[0:VM_N])/SM_Le_x*RM_Le_x)+np.sum(UY[VM_N:i])
            # else:
            dispRUY[i]=np.sum(dispUY[0:i])   

        #Displacement Y total
        dispYTotal = dispY+dispRUY
        #Y displacement along cantilever
        dispUYAlong = np.zeros(self.noNodes)
        #Along cantilever Rotational Displacement
        dispRXYAlong = np.zeros(self.noNodes)
        for i in range(1,self.noNodes):
            dispUYAlong[i]=np.sum(dispYTotal[0:i+1])
            dispRXYAlong[i]=np.sum(dispR[0:i+1])


        curvature = np.zeros(self.noNodes)
        do = ComputeCurvature()
        for i in range(1,self.noNodes-1,1):
            curvature[i] = do.curvature(
                xx=np.r_[self.nodeAlongCantileverCoords[i-1],self.nodeAlongCantileverCoords[i],self.nodeAlongCantileverCoords[i+1]], 
                yy=np.r_[dispUYAlong[i-1],dispUYAlong[i],dispUYAlong[i+1]]
                )
            
        #Euler-bernouilli 
        strainArr =-curvature*(self.cantilever['t']/2) 
        deformedLenAtSurfaceArr = (strainArr*self.baseElementXLenArr) + self.baseElementXLenArr
        
        self.mechOut = {
            'spring along' : springAlong,
            'spring total': 1/np.sum(1/springE),
            'disp along Y': dispUYAlong,
            'disp along R': dispRXYAlong,
            'disp end Y' : np.sum(dispYTotal),
            'strain along Y' : strainArr,
            'strain total' : np.sum(strainArr),
        }
        # print(1/curvature)
        print(self.mechOut['strain along Y'])
        print('spring const: ',self.mechOut['spring total'])
        print('mechanical disp: ',self.mechOut['disp end Y'])
        # print('disp: ',self.mechOut['disp along Y'])
        return self.mechOut


        

    def compositeMechanicalCalc(self):
        pass



    def derivitave1(self,eq,boundries,interval):
        t=np.linspace(interval[0],interval[-1], self.noElements+1)
        h = (interval[-1]-interval[0])/self.noElements

        """Input matrix A for [A][y]=[b]"""
        A[0,0] = 1
        A[self.noElements,self.noElements] = 1
        for i in range(1,self.noElements):
            A[i,i-1]= 1
            A[i,i] = 0
            A[i,i+1]= 1


        """output matrix b for [A][y]=[b]"""
        b = b + (eval(eq) if type(eq) is str else eq) * 2*h 

        for index,boundry in enumerate(boundries):
            if boundry is not None:
                b[int(index/h)] = boundry 
        

        print(A)
        print(b)
        return t, np.linalg.solve(A,b)



    def derivitave2(self, eq, boundries, xAxis=None,interval=None):
        A = np.zeros((self.noNodes,self.noNodes))
        b = np.zeros(self.noNodes)

        if interval is not None:
            t=np.linspace(interval[0],interval[-1], self.noNodes)   
        if xAxis is not None:
            t = xAxis
        h = (interval[-1]-interval[0])/self.noElements


        # if interval is None:
        #     [i for i, j in enumerate(boundries) if j is not None]

        """Input matrix A for [A][y]=[b]"""
        A[0,0] = 1
        A[self.noElements,self.noElements] = 1
        for i in range(1,self.noElements):
            A[i,i-1]= 1
            A[i,i] = -2
            A[i,i+1]= 1


        """output matrix b for [A][y]=[b]"""
        b = b + (eval(eq) if type(eq) is str else eq) * h**2 

        for index,boundry in enumerate(boundries):
            if boundry is not None:
                b[int(index/h)] = boundry 
        

        print(A)
        print(b)
        return t, np.linalg.solve(A,b)

class FDMCantileverSensitivity:
    def __init__(self, variableDict, variable, arr=None, lims=None, steps=100, **kwarg):
        if (arr is None and lims is None) or (arr is not None and lims is not None):
            raise NotImplementedError('Array OR limits with number steps, not both (or neither)!')

        if lims is not None:
            arr = np.linspace(lims[0],lims[1],steps)
        print(f"Sweeping variable {variableDict} {variable} in {arr} ")


        if arr is not None:
            out = []
            for i,val in enumerate(arr):
                kwarg[variableDict][variable] = val
                obj = FDMCantilever(**kwarg)
                out.append([val,obj.doMechanical()['disp end Y']])




        # do = FDMCantilever(**kwarg)

if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt

    settings = {
        "Number of Elements" :121,
    }

    settings.update({
        "Cantilever Width" : 120e-6,
        "Cantilever Length": 150e-6,
        "Cantilever Thickness" : 400e-9,
        "Cantilever Possion's Ratio" : 0.263,
        "Cantilever Youngs Modulus" : 2.4162e11,
    })
    settings.update({
        "Cantilever Taper" :{   "Taper" : True,
                                "start" : 90e-6,
                                "end"   : settings["Cantilever Length"],
                                "end width": 0e-6,
                            }
    })

    settings.update({
        "Part Angle Cantilever":{  "PartAngle": False,
                                    "start"  : 139e-6, 
                                    "End"    : settings["Cantilever Length"], 
                                    "degrees": 56, 
                                    "radians": np.deg2rad(56)
                                }
    })



    do = FDMCantilever(cantileverSettings=settings)
    eq = '-9'
    interval = [0,5]
    boundries = [None] * (interval[-1]+1)
    boundries[0]=0
    boundries[5]=50

    constants = None
    # print(boundries)
    # t,y = do.derivitave2(eq=eq, boundries=boundries, interval=interval)
    # print(t,y)



    # plt.figure(figsize=(10,8))
    # plt.plot(t,y)
    # # plt.plot(5, 50, 'ro')
    # plt.xlabel('time (s)')
    # plt.ylabel('altitude (m)')
    # plt.show()