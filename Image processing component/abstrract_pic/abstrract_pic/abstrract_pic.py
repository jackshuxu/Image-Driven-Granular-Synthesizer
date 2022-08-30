from PIL import Image, ImageStat

import numpy
image = Image.open("D:\python codes\pic source\canny_7.jpg")
image = image.convert("RGB")
array = numpy.array(image) 

print(array.shape)
print(array)

# R
print(array[:,:,0],array[:,:,0].shape)
# G
print(array[:,:,1],array[:,:,1].shape)
# B
print(array[:,:,2],array[:,:,2].shape)
file0 = open(r'D:\python codes\R7.txt', 'w',encoding='UTF-8')
for i in range (265):
    for j in range(265):
        file0.write(str(array[i,j,0])+'\n')
file0.close()

file1 = open(r'D:\python codes\G1.txt', 'w',encoding='UTF-8')
for i in range (423):
    for j in range(750):
        file1.write(str(array[i,j,1])+'\n')
file1.close()

file2 = open(r'D:\python codes\B1.txt', 'w',encoding='UTF-8')
for i in range (423):
    for j in range(750):
        file2.write(str(array[i,j,2])+'\n')
file2.close()

def brightness( im_file ):
   im = Image.open(im_file).convert('L')
   stat = ImageStat.Stat(im)
   return stat.mean[0]

if __name__ == "__main__":
   b =  brightness('D:\python codes\pic source\gray_1.jpg')
