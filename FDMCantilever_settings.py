
from FDMCantilever import FDMCantilever

if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt

    elements = {
        "Number of Elements" :121,
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
        "Part Angle Cantilever":{  "PartAngle": False,
                                    "start"  : 139e-6, 
                                    "End"    : cantileverSettings["Length"], 
                                    "degrees": 56, 
                                    "radians": np.deg2rad(56)
                                }
    })
    """
    ----------------------------Geometric Modifications--------------------------------
    """
    geomertyMods = {}  #can be of type 'gap' or 'taper'

    geoGapWidth ={
        "Gap/Width Change" : False
    }
    if geoGapWidth["Gap/Width Change"]:    #add as many gaps as you like
        geoGapWidth.update({

            "Change1" : {
                    "active"        :True,
                    "startX"        :8e-6,
                    "endX"          :20e-6,#cantileverSettings["Cantilever Length"],
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



    geoTaper = {
        "Taper" : False,
        "parameters"    : { "start" : 110e-6,
                            "end"   : cantileverSettings["Length"],
                            "end width": 0e-6,
                        }
    }

    if geoTaper["Taper"]:
        geomertyMods["Taper"] = geoTaper

    if geoGapWidth["Gap/Width Change"]:
        geomertyMods['Gaps/Width'] = geoGapWidth


    """
    ----------------------------Coating--------------------------------
    """
    
    coatings = {
        "Coating": True,
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

 
    noElements = elements['Number of Elements']



    do = FDMCantilever(cantileverSettings=cantileverSettings, noElements=noElements, geomertyMods=geomertyMods, coatings=coatings,coatingMods=coatingMods)
    eq = '-9'
    interval = [0,5]
    boundries = [None] * (interval[-1]+1)
    boundries[0]=0
    boundries[5]=50

    # print(boundries)
    # t,y = do.derivitave2(eq=eq, boundries=boundries, interval=interval)
    # print(t,y)



    # plt.figure(figsize=(10,8))
    # plt.plot(t,y)
    # # plt.plot(5, 50, 'ro')
    # plt.xlabel('time (s)')
    # plt.ylabel('altitude (m)')
    # plt.show()