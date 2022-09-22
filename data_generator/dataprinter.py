from turtle import *
import re
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
from tkinter import *
from PIL import Image


import pandas as pd


class Conductor():
    def __init__(self, number, left_low, right_high):
        self.number = number
        self.left_low = left_low
        self.right_high = right_high


def read_input_file(file_location="../Fieldsolver2d_hybrid_to_improve/input/input_9.data"):
    file = open(file_location)
    file_lines = file.readlines()
    dielectric = 0
    ConductorList = []
    for i in range(len(file_lines)):
        # line = file_lines[i].split(" ")
        line = re.split(r'[ ]+', file_lines[i])
        if (line[0] == "dielectric"):
            dielectric = float(line[1][0:-1])
        elif (i == 0):
            continue
        elif (line[-1] == "\n"):
            continue
        elif (line[-1][-1:] == "\n" and len(line[-1]) > 2):
            number = line[1][3:]
            left_low = [float(line[2]), float(line[3])]
            right_high = [float(line[4]), float(line[5][0:-1])]
            ConductorList.append(Conductor(number, left_low, right_high))
        else:
            number = line[1][3:]
            left_low = [float(line[2]), float(line[3])]
            right_high = [float(line[4]), float(line[5])]
            ConductorList.append(Conductor(number, left_low, right_high))

    return dielectric, ConductorList, file_location[47:-5]


# def print_input_file(ConductorList):
#     # 画图
#     img = plt.imread("img.jpg")
#     h, w, c = img.shape
#     fig, ax = plt.subplots(1, 1)
#     ax.imshow(img)
#     currentAxis = fig.gca()
#
#     plt_list = []
#     for i in range(len(ConductorList)):
#         x = ConductorList[i].left_low[0] * 10000
#         y = ConductorList[i].left_low[1] * 10000
#         width = ConductorList[i].right_high[0] * 10000 - ConductorList[i].left_low[0] * 10000
#         hight = ConductorList[i].right_high[1] * 10000 - ConductorList[i].left_low[1] * 10000
#         plt_list.append([x.__round__(0), y.__round__(0), width.__round__(0), hight.__round__(0)])
#         # plt_list.append([x,y,width,hight])
#
#     for i in range(len(plt_list)):
#         rec = patches.Rectangle((plt_list[i][0], plt_list[i][1]), plt_list[i][2], plt_list[i][3])
#         currentAxis.add_patch(rec)
#         plt.show()
#
#     # plt.show()
#     print(plt_list)


def print_img(ConductorList):
    fig = plt.figure(figsize=(20, 10))

    ax = fig.add_subplot(111)
    plt.axis('on')
    plt_list = []
    for i in range(len(ConductorList)):
        x = ConductorList[i].left_low[0]
        y = ConductorList[i].left_low[1]
        width = ConductorList[i].right_high[0] - ConductorList[i].left_low[0]
        hight = ConductorList[i].right_high[1] - ConductorList[i].left_low[1]
        plt_list.append([x.__round__(4), y.__round__(4), width.__round__(4), hight.__round__(4)])

    print(plt_list)

    plt_list = np.array(plt_list)

    # todo 归一化 画到图像上
    # plt_list = 1/(1+np.exp(-plt_list))
    # plt_list = np.exp(plt_list)/np.sum(np.exp(plt_list))
    # plt_list[:, 1] = (plt_list[:, 1] - np.min(plt_list[:, 1])) / (np.max(plt_list[:, 1]) - np.min(plt_list[:, 1]))
    # plt_list[:, 2] = (plt_list[:, 2] - np.min(plt_list[:, 2])) / (np.max(plt_list[:, 2]) - np.min(plt_list[:, 2]))
    # plt_list[:, 3] = (plt_list[:, 3] - np.min(plt_list[:, 3])) / (np.max(plt_list[:, 3]) - np.min(plt_list[:, 3]))

    print(plt_list)

    for i in range(len(plt_list)):
        rect = plt.Rectangle((plt_list[i][0], plt_list[i][1]), plt_list[i][2], plt_list[i][3])
        ax.add_patch(rect)

    plt.show()


def tkinter_print(ConductorList,number):
    tk = Tk()
    canvas = Canvas(tk, width=4000, height=2000)
    canvas.pack()

    for i in range(len(ConductorList)):
        print(ConductorList[i].left_low,ConductorList[i].right_high)
        x1, y1 = (ConductorList[i].left_low[0]+10)*120, ConductorList[i].left_low[1]*100
        x2, y2 = (ConductorList[i].right_high[0]+10)*120, ConductorList[i].right_high[1]*100
        print(x1, y1, x2, y2)
        canvas.create_rectangle(x1, y1, x2, y2)
        canvas.update()
    canvas.postscript(file="../Fieldsolver2d_hybrid_to_improve/input/pic/"+number+".eps",colormode='color')
    img = Image.open("../Fieldsolver2d_hybrid_to_improve/input/pic/"+number+".eps")
    img.save("../Fieldsolver2d_hybrid_to_improve/input/pic/"+number+".png", "png")
    tk.mainloop()


if __name__ == '__main__':
    dielectric, ConductorList ,number= read_input_file()

    # print_input_file(ConductorList)
    # test_Conductor(ConductorList)
    # print_img(ConductorList)

    tkinter_print(ConductorList,number)
