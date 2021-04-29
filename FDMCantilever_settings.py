
from FDMCantilever import FDMCantilever, FDMCantileverSensitivity

if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt

    elements = {
        "Number of Elements" :500,
    }

    """
    ----------------------------Cantilever settings---------------------------------
    """
    cantileverSettings = ({
        "Width" : 120e-6,
        "Length": 150e-6,
        "Thickness" : 400e-9,
        "Possion's Ratio" : 0.263,
        "Youngs Modulus" : 2.4162e11,
    })

    
    cantileverSettings.update({
        'angle cantilever': False,
        'angle degrees const' : 13,
        'angle change': False,
        'changing angle function' : lambda x: (1.4*x**0.4), #func as lambda

        "Part Angle Cantilever":{  "Part Angle": False,
                                    "start"  : 139e-6, 
                                    "end"    : cantileverSettings["Length"], 
                                    "degrees": 56, 
                                }
    })
    
    """
    ---------------------------Mechanical----------------
    """
    mechanical = {
        "Force": 1e-9,
        "Offset Force": False,
        "Force application X point": 106E-6
    }
    """
    ----------------------------Geometric Modifications--------------------------------
    """
    geomertyMods = {}  #can be of type 'gap' or 'taper'

    geoGapWidth ={
        "Gap/Width Change" : True
    }
    if geoGapWidth["Gap/Width Change"]:    #add as many gaps as you like
        geoGapWidth.update({

            "Change1" : {
                    "active"        :True,
                    "startX"        :8e-6,
                    "endX"          :20e-6,#cantileverSettings["Cantilever Length"],
                    "start change"  :100e-6,
                    "end change"    :1e-6,
                    },

            "Change2" : {
                    "active"        :False,
                    "startX"        :110e-6,
                    "endX"          :cantileverSettings["Length"],
                    "start change"  :119e-6,
                    "end change"    :0e-6,
                    },

        })



    geoTaper = {
        "Taper" : True,
        "parameters"    : { "start" : 90E-6,
                            "end"   : cantileverSettings["Length"],
                            "end width": 0e-6,
                        }
    }


    """
    ----------------------------Coating--------------------------------
    """
    
    coatings = {
        "Coating": False,
        "startX" : 0,
        "endX"   : cantileverSettings['Length'],
        "Thickness" : 70-9,
        "width" : cantileverSettings['Width'],
        "Possion's Ratio" : 0.45,
        "Youngs Modulus" : 64.39e+9,
    }

    coatingMods = {}  #can be of type 'gap' or 'taper'

    coatGapsWidths = {
        "Gap/Width Change" : False
    }
    coatGapsWidths.update({  #add as many gaps as you like
            "Change1" : {
                    "active"        :True,
                    "startX"        :8e-6,
                    "endX"          :20e-6,
                    "start change"  :10e-6,
                    "end change"    :10e-6,
                    },

            "Change2" : {
                    "active"        :True,
                    "startX"        :110e-6,
                    "endX"          :cantileverSettings["Length"],
                    "start change"  :119e-6,
                    "end change"    :0e-6,
                    },  
    })
    coatTaper = {
        "Taper" : False,
        "parameters"    : { "start" : 110e-6,
                            "end"   : cantileverSettings["Length"],
                            "end width": 0e-6,
                        }
    }


    if coatTaper["Taper"]:
        coatingMods["Taper"] = coatTaper

    if coatGapsWidths["Gap/Width Change"]:
        coatingMods['Gaps/Width'] = coatGapsWidths

    if geoTaper["Taper"]:
        geomertyMods["Taper"] = geoTaper

    if geoGapWidth["Gap/Width Change"]:
        geomertyMods['Gaps/Width'] = geoGapWidth

    cantileverSettings['angle rad'] = np.deg2rad(cantileverSettings['angle degrees const']) if cantileverSettings['angle cantilever'] else 0
    cantileverSettings['Part Angle Cantilever']['rad'] = np.deg2rad(cantileverSettings['Part Angle Cantilever']['degrees']) if cantileverSettings['Part Angle Cantilever']['Part Angle'] else 0

    # noElements = elements['Number of Elements']


    variableDict = "mechanical"
    variable = "Force"
    arr = [1e-9,2e-9]
    # do = FDMCantilever(cantileverSettings=cantileverSettings, noElements=elements, geomertyMods=geomertyMods, coatings=coatings,coatingMods=coatingMods, mechanical=mechanical)
    do = FDMCantileverSensitivity(variableDict=variableDict,variable=variable,arr=arr,cantileverSettings=cantileverSettings, noElements=elements, geomertyMods=geomertyMods, coatings=coatings,coatingMods=coatingMods, mechanical=mechanical)
    # eq = '-9'
    # interval = [0,5]
    # boundries = [None] * (interval[-1]+1)
    # boundries[0]=0
    # boundries[5]=50

    # print(boundries)
    # t,y = do.derivitave2(eq=eq, boundries=boundries, interval=interval)
    # print(t,y)



    # plt.figure(figsize=(10,8))
    # plt.plot(t,y)
    # # plt.plot(5, 50, 'ro')
    # plt.xlabel('time (s)')
    # plt.ylabel('altitude (m)')
    # plt.show()