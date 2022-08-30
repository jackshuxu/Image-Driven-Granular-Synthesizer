#-*-coding:utf-8-*-

import numpy as np
f = open(r"D:\python codes\R7.txt")
line = f.readline()
data_list = []
while line:
    num = list(map(float,line.split()))
    data_list.append(num)
    line = f.readline()
f.close()
data_array = np.array(data_list)
print (np.sum(data_array))
