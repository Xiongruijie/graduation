import gym
import os
import cv2
import tempfile
import shutil
import paddle
from paddle.io import Dataset
from paddle.vision.transforms import Normalize


class CircuitDataset(Dataset):
    def __init__(self, end_index, pic_path='../Fieldsolver2d_hybrid_to_improve/input/pic/',
                 transform=Normalize(mean=[127.5], std=[127.5], data_format='CHW')):
        super(CircuitDataset, self).__init__()
        self.data_list = []
        for i in range(end_index):
            file_path = pic_path + i.__str__() + '.png'
            self.data_list.append(file_path)
        self.transform = transform

    def __getitem__(self, index):
        image_path = self.data_list[index]
        image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
        image = image.astype('float32')
        if self.transform is not None:
            image = self.transform(image)
        return image

    def __len__(self):
        return len(self.data_list)


class CircuitEnv(gym.Env):
    def __init__(self, local_index, init_index):
        # 读取图片位置
        self.root = '../Fieldsolver2d_hybrid_to_improve/input/pic/'
        self.local_index = local_index
        # self.reward = None
        self.dataset = CircuitDataset(init_index)
        self.state = self.dataset.__getitem__(self.local_index)

    def reward(self, action, index):
        ############################
        action = paddle.tensor.to_tensor([0, 0, 1], dtype='float32')
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
        output_path = '../Fieldsolver2d_hybrid_to_improve/output/'+choose_action+'/output_'+index.__str__()+'.out'


        # return paddle.argmax(action)

    def step(self, action):
        next_state = self.dataset.__getitem__(self.local_index + 1)
        # 计算action 输入然后查找对应的解
        reward = self.reward(action, self.local_index)
        self.state = next_state
        self.local_index = self.local_index + 1
        return next_state, reward


if __name__ == '__main__':
    env = CircuitEnv(0, 10)
    print(env.reward([], 0))
