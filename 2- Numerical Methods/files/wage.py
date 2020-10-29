import numpy as np


bracket_16t19  = 337*52 # multiply by 52 to get yearly wages
bracket_20t24  = 450*52
bracket_25t34  = 643*52
bracket_35t44  = 769*52
bracket_45t54  = 790*52
bracket_55t64  = 803*52
bracket_65plus = 605*52


wageprofile    = np.array([bracket_16t19, bracket_20t24,
                 bracket_25t34,bracket_35t44,
                 bracket_45t54,bracket_55t64,
                 bracket_65plus])

agebrackets    = np.array([16+19,20+24,25+34,35+44,45+54,55+64,65+90])/2


xi = agebrackets
yi = wageprofile