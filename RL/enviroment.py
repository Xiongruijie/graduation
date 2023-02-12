import math
import random
import re
from time import sleep

import gym
import os
import cv2
import tempfile
import shutil

import numpy as np
import paddle
import torch
import torchvision
from matplotlib import pyplot as plt
from paddle.io import Dataset
from paddle.vision.transforms import Normalize


# end_index 输入
class CircuitDataset(Dataset):
    def __init__(self, start_index, end_index, pic_path='../Fieldsolver2d_hybrid_to_improve/input/pic/',
                 transform=Normalize(mean=[127.5], std=[127.5], data_format='CHW')):
        super(CircuitDataset, self).__init__()
        self.data_list = []
        for i in range(start_index, end_index):
            file_path = pic_path + i.__str__() + '.png'
            self.data_list.append(file_path)
        self.transform = transform

    def __getitem__(self, index):
        image_path = self.data_list[index]
        image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
        image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
        # image = image.astype('float32')
        # if self.transform is not None:
        #     image = self.transform(image)
        image = paddle.vision.transforms.to_tensor(image, data_format='CHW').unsqueeze(0)
        # transform1 = torchvision.transforms.Compose([
        #     torchvision.transforms.ToTensor()
        # ])

        # image = transform1(image)
        # # image = torch.tensor(image,dtype=torch.float32).unsqueeze(0)

        return image

    def __len__(self):
        return len(self.data_list)


class CircuitEnv(gym.Env):
    def __init__(self, local_index, init_start_index, init_end_index, prefer_accuracy=False, prefer_time=False,
                 prefer_memory=False):
        """

        :param local_index: 当前游标
        :param init_start_index: 游标起始位置
        :param init_end_index: 游标结束位置
        :param prefer_accuracy: 精确度偏好
        :param prefer_time: 时间偏好
        :param prefer_memory:  内存偏好
        """
        # 读取图片位置
        self.init_start_index = init_start_index
        self.init_end_index = init_end_index
        self.action_space = 3
        self.root = '../Fieldsolver2d_hybrid_to_improve/input/pic/'
        self.local_index = random.randint(self.init_start_index, self.init_end_index)
        # self.reward = None
        self.dataset = CircuitDataset(init_start_index, init_end_index)
        self.state = self.dataset.__getitem__(self.local_index)
        self.prefer_accuracy = prefer_accuracy
        self.prefer_time = prefer_time
        self.prefer_memory = prefer_memory

    def get_state(self):
        return self.dataset.__getitem__(self.local_index)

    def reward(self, action, index):
        ############################
        # action = paddle.tensor.to_tensor([0, 0, 1], dtype='float32')
        fdm = paddle.tensor.to_tensor([0], dtype='int64')
        bem = paddle.tensor.to_tensor([1], dtype='int64')
        frw = paddle.tensor.to_tensor([2], dtype='int64')

        choose_action = None
        ############################
        if paddle.equal(paddle.argmax(action), fdm):
            choose_action = 'fdm'
        elif paddle.equal(paddle.argmax(action), bem):
            choose_action = 'bem'
        elif paddle.equal(paddle.argmax(action), frw):
            choose_action = 'frw'
        output_path = '../Fieldsolver2d_hybrid_to_improve/output/' + choose_action + '/output_' + index.__str__() + '.out'
        compare_path = '../Fieldsolver2d_hybrid_to_improve/output/' + 'fdm' + '/output_' + index.__str__() + '.out'

        # 得分函数：
        # 输入：output_path = 选择的算法得出的计算结果文件位置
        #       compare_path = fdm基准结果文件位置
        accuracy, time, memory, net_number = judge_by_file(output_path, compare_path)
        accuracy_score = 0
        # 精确度得分：

        if self.prefer_accuracy is True:
            accuracy_increase = 5
        else:
            accuracy_increase = 1

        if accuracy < 0.50:
            accuracy_score += (0 * accuracy_increase)
        elif accuracy >= 0.50 and accuracy < 0.60:
            accuracy_score += (6 * accuracy_increase)
        elif accuracy >= 0.60 and accuracy < 0.70:
            accuracy_score += (7.5 * accuracy_increase)
        elif accuracy >= 0.70 and accuracy < 0.80:
            accuracy_score += (8.5 * accuracy_increase)
        elif accuracy >= 0.90 and accuracy < 0.95:
            accuracy_score += (9 * accuracy_increase)
        elif accuracy >= 0.95 and accuracy < 1:
            accuracy_score += (9.5 * accuracy_increase)
        elif accuracy == 1:
            accuracy_score += (9.6 * accuracy_increase)
        # 时间得分：
        if self.prefer_time is True:
            time_increase = 5
        else:
            time_increase = 1
        every_net_time = time / net_number
        time_score = math.log(((1 / every_net_time) * 1.0e5 * net_number) + 1, 2) * 2 * time_increase

        # 内存得分：
        if self.prefer_memory is True:
            memory_increase = 5
        else:
            memory_increase = 1

        every_net_memory = memory / net_number
        memory_score = math.log((1 / every_net_memory) * 1e4 + 1, 2) * 5 * memory_increase

        score = accuracy_score + memory_score + time_score

        return score

    def step(self, action):

        # 计算action 输入然后查找对应的解
        is_terminal = False
        action = paddle.to_tensor(action, dtype='float32')
        reward = self.reward(action, self.local_index)

        if self.local_index + 1 < self.init_end_index:
            next_state = self.dataset.__getitem__(self.local_index + 1)
            self.state = next_state
        else:
            is_terminal = True
            next_state = self.state
        self.local_index = self.local_index + 1

        return next_state, reward, is_terminal

    def reset(self, **kwargs):
        # self.local_index = self.init_start_index
        self.local_index = random.randint(0,88)
        return self.get_state()


def judge_by_file(output_path, compare_path):
    baseline_file = open(compare_path)
    selected_file = open(output_path)
    baseline_lines = baseline_file.readlines()
    selected_lines = selected_file.readlines()
    # 处理一下数据把前面的num0：去掉，把后面的换行符去掉
    baseline_data = re.split(r'[ ]+', baseline_lines[3])
    selected_data = re.split(r'[ ]+', selected_lines[3])
    del baseline_data[0:2]
    del selected_data[0:2]
    baseline_data_re = []
    selected_data_re = []
    for sub in baseline_data:
        baseline_data_re.append(sub.replace("\n", ""))
    for sub in selected_data:
        selected_data_re.append(sub.replace("\n", ""))

    for i in range(len(baseline_data_re)):
        baseline_data_re[i] = float(baseline_data_re[i][0:-2])
    for i in range(len(selected_data_re)):
        selected_data_re[i] = float(selected_data_re[i][0:-2])
    accuracy_list = []
    # for i in range(len(baseline_data_re)):
    #     accuracy_list.append(abs(baseline_data_re[i] - selected_data_re[i]) / baseline_data_re[i])
    accuracy_list.append(abs(baseline_data_re[0] - selected_data_re[0]) / (baseline_data_re[0] + 1e-8))
    accuracy = 1 - accuracy_list[0]

    # 时间、空间
    time = int(re.split(r'[ ]+', selected_lines[6])[3])
    memory = int(re.split(r'[ ]+', selected_lines[7])[3])

    baseline_file.close()
    selected_file.close()
    net_number = len(baseline_data_re)
    return accuracy, time, memory, net_number


def file_read(path):
    baseline_file = open(path)
    baseline_lines = baseline_file.readlines()
    baseline_data = re.split(r'[ ]+', baseline_lines[3])
    del baseline_data[0:2]
    time = int(re.split(r'[ ]+', baseline_lines[6])[3])
    memory = int(re.split(r'[ ]+', baseline_lines[7])[3])
    baseline_file.close()
    return time, memory


def paint_pic():
    fdm_timeList = []
    fdm_memoryList = []
    bem_timeList = []
    bem_memoryList = []
    frw_timeList = []
    frw_memoryList = []
    for i in range(10):
        fdm_data = '../Fieldsolver2d_hybrid_to_improve/output/' + 'fdm' + '/output_' + i.__str__() + '.out'
        bem_data = '../Fieldsolver2d_hybrid_to_improve/output/' + 'bem' + '/output_' + i.__str__() + '.out'
        frw_data = '../Fieldsolver2d_hybrid_to_improve/output/' + 'frw' + '/output_' + i.__str__() + '.out'
        t1, m1 = file_read(fdm_data)
        t2, m2 = file_read(bem_data)
        t3, m3 = file_read(frw_data)
        fdm_timeList.append(t1)
        fdm_memoryList.append(m1)
        bem_timeList.append(t2)
        bem_memoryList.append(m2)
        frw_timeList.append(t3)
        frw_memoryList.append(m3)
    return fdm_timeList, fdm_memoryList, bem_timeList, bem_memoryList, frw_timeList, frw_memoryList


def paint_bar(data, pic_name, pic_title):
    N = 10
    plt.bar(x=list(range(N)), height=data)
    plt.title(pic_title)
    plt.show()
    plt.savefig(pic_name)


def paint_bar_for_10():
    fdm_timeList, fdm_memoryList, bem_timeList, bem_memoryList, frw_timeList, frw_memoryList = paint_pic()
    paint_bar(fdm_memoryList, 'fdm_memoryList.jpg', 'fdm_memory_consume')
    paint_bar(fdm_timeList, 'fdm_timeList.jpg', 'fdm_time_consume')
    paint_bar(bem_timeList, 'bem_timeList.jpg', 'bem_time_consume')
    paint_bar(bem_memoryList, 'bem_memoryList.jpg', 'bem_memory_consume')
    paint_bar(frw_timeList, 'frw_timeList.jpg', 'frw_time_consume')
    paint_bar(frw_memoryList, 'frw_memoryList.jpg', 'frw_memory_consume')


if __name__ == '__main__':
    dataset = CircuitDataset(0, 99)
    for i in range(99):
        print(dataset.__getitem__(i).shape)
    # print(env.reward([],index=6))
    # out = '../Fieldsolver2d_hybrid_to_improve/output/bem/output_6.out'
    # cmp = '../Fieldsolver2d_hybrid_to_improve/output/fdm/output_6.out'
    # judge_by_file(out, cmp)
