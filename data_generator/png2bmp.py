import fnmatch
import os
import pandas as pd
import numpy as np
import cv2
def ReadSaveAddr(Stra,Strb):
    # Stra:原始图片的位置  Strb:图片格式
    print("Read :",Stra,Strb)
    # 过滤所有符合要求的图片
    a_list = fnmatch.filter(os.listdir(Stra),Strb)
    print("Find = ",len(a_list))
    df = pd.DataFrame(np.arange(len(a_list)).reshape((len(a_list),1)),columns=['Addr'])
    df.Addr = a_list
    for i in range(len(a_list)):
        path = Stra+'/'+a_list[i]
        print(path)
        # 开始读取
        img = cv2.imread(path,cv2.IMREAD_GRAYSCALE)
        # 部分读取
        # temp = np.ones((255,255))
        # for m in range(255):
        #     for n in range(255):
        #         temp[m,n] = img[m,n]
        # print(temp.shape)
        t = a_list[i]
        t = t[:-4] #获取图片名
        t = '../Fieldsolver2d_hybrid_to_improve/input/pic/'+t+'.bmp'
        cv2.imwrite(t,img) #写入
    df.to_csv('Get.lst',columns=['Addr'],index=False,header=False)
    print("OK!")
# 读取该目录下所有后缀为.png的图片
ReadSaveAddr("../Fieldsolver2d_hybrid_to_improve/input/pic","*.png")