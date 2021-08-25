
mass =  {
    "C"     : 12.000000, #C
    "H"     : 1.007825,  #H
    "O"     : 15.994915, #O
    "N"     : 14.003074, #N
    "S"     : 31.972072, #S
    "P"     : 30.973763, #P
    "N15"   : 15.000109, #N15
    "H2"    : 2.014102,  #H2
    "C13"   : 13.003355, #C13
    "Hplus" : 1.0072765, #Hplus
    "Cl"    : 34.96885,  #Chlorine
    "Br"    : 78.91833,  #Bromine
    "Se"    : 79.9165213 #Selenium
}

heavy_map = {
    "C"     : 1.10,  #C
    "H"     : 0.015, #H
    "O"     : 0.2,   #O
    "N"     : 0.37,  #N
    "S"     : 4.21,  #S
    "P"     : 0.0,   #P
    "N15"   : 100.0, #N15
    "H2"    : 100.0, #H2
    "C13"   : 100.0, #C13
    "Hplus" : 100.0, #Hplus
    "Cl"    : 49.16, #Chlorine
    "Br"    : 24.24, #Bromine
    "Se"    : 49.31  #Selenium
}

light_map = {}
for e in heavy_map.keys():
    heavy_map[e] /= 100.0
    light_map[e] = 1-heavy_map[e]

