import numpy as np
import mpmath as mp


class FDMCantilever:
    def __init__(self, settingsList, cantileverModifications={}):

        self.cantilever = {
            'w'  : settingsList["Cantilever Width"],
            'l'  : settingsList["Cantilever Length"],
            't'  : settingsList["Cantilever Thickness"],
            'E'  : settingsList["Cantilever Youngs Modulus"],
            'vs' : settingsList["Cantilever Possion's Ratio"],
        }

        self.noElements = settingsList["Number of Elements"]
        self.noNodes = self.noElements+1
        self.elementLength = self.cantilever["l"]/self.noElements



        self.nodeXCoords = np.arange(0,self.noNodes,1)*self.elementLength #TODO: variable mesh
        self.nodeCentroidCoords = self.nodeXCoords-(self.elementLength/2) #TODO: variable mesh
        self.nodeCentroidCoords[0] = 0
        
        self.nodeAlongCantileverCoords = self.nodeXCoords if not settingsList["Part Angle Cantilever"]["PartAngle"] else np.zeros(self.noNodes)

        self.baseElementXLenArr = np.ones(self.noNodes)*self.cantilever["l"]
        self.baseElementXLenArr[0] = 0
        self.baseCantileverXLen = self.cantilever['w']

        #pristmatic widths
        self.widthArr = np.ones(self.noNodes)*self.cantilever['w']
        
        #TODO: Variable mesh

        if cantileverModifications:
            self.modifyGeometry(cantileverModifications)

        self.element2ndMomentOfArea = self.widthArr*self.cantilever['t']**3/12
        self.element2ndMomentOfArea[0] = np.Infinity


    def modifyGeometry(self,modDict):
        if 'Taper' in modDict.keys() :
            if  modDict['Taper']['Taper']:
                taper = modDict['Taper']['parameters']

                TM=(self.cantilever['w']-taper['end width'])/(taper['start']-taper['end'])

                # TS_N = ([i for i in range(0,self.noElements) if self.nodeCentroidCoords[i]<=taper['start'] and self.nodeCentroidCoords[i+1]>taper['start']])[0]
                # print(TS_N)
                TS_N = next(index-1 for index,val in enumerate(self.nodeCentroidCoords) if val>=taper['start'])



                for i in range(TS_N,self.noNodes,1):

                    if self.nodeCentroidCoords[i]>taper['start'] and self.nodeCentroidCoords[i]<=taper['end']:
                        self.widthArr[i]=(TM * (self.nodeCentroidCoords[i]-taper['start']) ) + self.cantilever['w']
                        
                    elif self.nodeCentroidCoords[i]>taper['end']:
                        self.widthArr[i]=taper['end width']

        
        if 'Gaps/Width' in modDict.keys() :
            if modDict['Gaps/Width']["Gap/Width Change"]:
                params = modDict['Gaps/Width']
                params.pop('Gap/Width Change')

                for key in params.keys():
                    width = params[key]
                    if not width['active']: 
                        continue 
                    GM=(width['end change']-width['start change'])/(width['endX']-width['startX'])

                    GS_N=([i for i in range(0,self.noElements,1) if self.nodeCentroidCoords[i]<=width['startX'] and self.nodeCentroidCoords[i+1]>width['startX']])[0]
                    GE_N=([i if self.nodeCentroidCoords[i-1]<=width['endX'] and self.nodeCentroidCoords[i]>width['endX'] else self.noNodes for i in range(GS_N,self.noNodes,1)])[0] 
                    print(GS_N)
                    print(GE_N)
                    GS_N=next(index-1 for index,val in enumerate(self.nodeCentroidCoords) if val>=width['startX'])
                    GE_N=next((index-1 for index,val in enumerate(self.nodeCentroidCoords) if val>=width['endX']),self.noNodes)
                    print(GS_N)
                    print(GE_N)

                    for i in range (GS_N,GE_N,1):
                        if self.nodeCentroidCoords[i]>width['startX'] and self.nodeCentroidCoords[i]<=width['endX']:
                            self.widthArr[i]=self.widthArr[i]-width['start change']-(GM*(self.nodeCentroidCoords[i]-width['startX']))




        


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



    do = FDMCantilever(settingsList=settings)
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