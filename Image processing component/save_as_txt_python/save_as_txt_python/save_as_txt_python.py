
import numpy as np
a=np.array([[7.20, 152.45],
         [4.76, 67.69],
         [6.41, 243.42],
         [20.63, 243.55],
         [13.89, 249.13],
         [26.66, 233.74],
         [35.71, 213.02],
         ])
np.savetxt('myfile.txt', a, fmt="%.2lf", delimiter=" ")



