import numpy as np
import mpmath as mp


class FDMCantilever:
    def __init__(self, cantileverSettings, noElements=100,geomertyMods={}, coatings={},coatingMods={}):

        self.cantilever = {
            'w'  : cantileverSettings["Width"],
            'l'  : cantileverSettings["Length"],
            't'  : cantileverSettings["Thickness"],
            'E'  : cantileverSettings["Youngs Modulus"],
            'vs' : cantileverSettings["Possion's Ratio"],
        }
        

        self.noElements = noElements
        self.noNodes = self.noElements+1
        self.elementLength = self.cantilever["l"]/self.noElements

        #setup node arrays
        self.nodeXCoords = np.arange(0,self.noNodes,1)*self.elementLength #TODO: variable mesh
        self.nodeCentroidCoords = self.nodeXCoords-(self.elementLength/2) #TODO: variable mesh
        self.nodeCentroidCoords[0] = 0
        
        self.nodeAlongCantileverCoords = self.nodeXCoords if not cantileverSettings["Part Angle Cantilever"]["PartAngle"] else np.zeros(self.noNodes)

        self.baseElementXLenArr = np.ones(self.noNodes)*self.cantilever["l"]
        self.baseElementXLenArr[0] = 0
        self.baseCantileverXLen = self.cantilever['w']

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

        

    def modify(self, modDict, arr):
        if 'Taper' in modDict.keys() :
            if  modDict['Taper']['Taper']:
                self.addTaper(arr=arr, params=modDict['Taper']['parameters'])
        
        if 'Gaps/Width' in modDict.keys() :
            if modDict['Gaps/Width']["Gap/Width Change"]:
                params = modDict['Gaps/Width']
                params.pop('Gap/Width Change')

                for key in params.keys():
                    width = params[key]
                    if not width['active']: 
                        continue 
                    self.addWidthChange(arr=arr,params=width)
        

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


    

    def derivitave1(self,eq,boundries,interval):
        t=np.linspace(interval[0],interval[-1], self.noElements+1)
        h = (interval[-1]-interval[0])/self.noElements



    def derivitave2(self, eq, boundries, interval):
        A = np.zeros((self.noNodes,self.noNodes))
        b = np.zeros(self.noNodes)

        t=np.linspace(interval[0],interval[-1], self.noNodes)
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