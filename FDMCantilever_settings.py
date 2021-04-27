
from Numerical_FDM_solver import FDMCantilever

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
        "Part Angle Cantilever":{  "PartAngle": False,
                                    "start"  : 139e-6, 
                                    "End"    : settings["Cantilever Length"], 
                                    "degrees": 56, 
                                    "radians": np.deg2rad(56)
                                }
    })

    modifications = {}  #can be of type 'gap' or 'taper'

    gaps ={
        "Gap/Width Change" : True
    }
    if gaps["Gap/Width Change"]:    #add as many gaps as you like
        gaps.update({

            "Gap1" : {
                    "active"        :True,
                    "startX"        :8e-6,
                    "endX"          :10e-6,#settings["Cantilever Length"],
                    "start change"  :10e-6,
                    "end change"    :10e-6,
                    },

            "Gap2" : {
                    "active"        :True,
                    "startX"        :110e-6,
                    "endX"          :settings["Cantilever Length"],
                    "start change"  :119e-6,
                    "end change"    :0e-6,
                    },

        })



    taper = {
        "Taper" : False,
        "parameters"    : { "start" : 110e-6,
                            "end"   : settings["Cantilever Length"],
                            "end width": 0e-6,
                        }
    }

    if taper["Taper"]:
        modifications["Taper"] = taper

    if gaps["Gap/Width Change"]:
        modifications['Gaps/Width'] = gaps

 



    do = FDMCantilever(settingsList=settings, cantileverModifications=modifications)
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