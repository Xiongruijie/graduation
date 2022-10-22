import random
import time
from collections import deque
import gym
import numpy as np
import paddle
import torch
from torch import nn
from torch.distributions import Normal

from RL import enviroment


class Actor(nn.Module):
    def __init__(self, action_dim=3, is_train=True):
        super(Actor, self).__init__()
        self.action_dim = action_dim
        self.is_train = is_train
        self.conv = nn.Sequential(
            nn.Conv2d(3, 256, 11, 4),
            nn.ReLU(),
            nn.MaxPool2d(3, 3),
            nn.Conv2d(256, 128, 3, 1, 1),
            nn.MaxPool2d(3, 3),
            nn.ReLU(),
            nn.Conv2d(128, 32, 3, 1, 1),
            nn.MaxPool2d(3, 3),
        )
        self.fc = nn.Sequential(
            nn.Linear(4896, 4096),  # 1* 3* 1031*1865
            # nn.Linear(10368,4096),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(4096, 2048),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(2048, 256),
            nn.ReLU(),
            nn.Linear(256, self.action_dim),
        )
        self.noisy = Normal(0, 0.2)

    def forward(self, obs):
        obs = obs.unsqueeze(0)
        feature = self.conv(obs)
        feature = feature.reshape(obs.shape[0], -1)
        output = self.fc(feature)
        return output

    def select_action(self, epsilon, state):
        state = torch.tensor(state, dtype=torch.float32)
        with torch.no_grad():
            action = self.forward(state)
            action = action + self.is_train * epsilon * self.noisy.sample(action.shape)
        return action.numpy()
        # return action.shape


# In[3]:


class Critic(nn.Module):
    def __init__(self, action_dim=3, ):
        super(Critic, self).__init__()
        self.action_dim = action_dim
        self.conv = nn.Sequential(
            nn.Conv2d(3, 256, 11, 4),
            nn.ReLU(),
            nn.MaxPool2d(3, 3),
            nn.Conv2d(256, 128, 3, 1, 1),
            nn.MaxPool2d(3, 3),
            nn.ReLU(),
            nn.Conv2d(128, 32, 3, 1, 1),
            nn.MaxPool2d(3, 3),
        )
        # 因为中间要加动作判定所以需要fc分层
        self.fc = nn.Sequential(
            nn.Linear(4896, 4096),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(4096, 2048),
            nn.ReLU(),
            nn.Dropout(0.5),
        )
        self.fc2 = nn.Sequential(
            nn.Linear(2048 + 3, 1),
        )

    def forward(self, obs, action):
        obs = obs.unsqueeze(0)
        feature = self.conv(obs)
        feature = feature.reshape([obs.shape[0], -1])
        output = self.fc(feature)

        output = paddle.concat((output, action), axis=1)
        output = self.fc2(output)

        return output


# class Memory(object):
#     def __init__(self, memory_size) -> None:
#         self.memory_size = memory_size
#         self.buffer = deque(maxlen=self.memory_size)
#
#     def add(self, experience) -> None:
#         self.buffer.append(experience)
#
#     def size(self):
#         return len(self.buffer)
#
#     def sample(self, batch_size: int, continuous: bool = False):
#         if batch_size > len(self.buffer):
#             batch_size = len(self.buffer)
#         if continuous:
#             rand = random.randint(0, len(self.buffer) - batch_size)
#             return [self.buffer[i] for i in range(rand, rand + batch_size)]
#         else:
#             indexes = np.random.choice(np.arange(len(self.buffer)), size=batch_size, replace=False)
#             return [self.buffer[i] for i in indexes]
#
#     def clear(self):
#         self.buffer.clear()


class DDPG(object):
    def __init__(self, state_dim, action_dim, replacement, memory_capacity=50000, gamma=0.9, lr_a=1e-4,
                 lr_c=3e-4, batch_size=32):
        super(DDPG, self).__init__()
        self.action_dim = action_dim
        self.memory_capacity = memory_capacity
        self.replacement = replacement
        self.t_replace_counter = 0
        self.gamma = gamma
        self.lr_a = lr_a
        self.lr_c = lr_c
        self.batch_size = batch_size
        self.state_dim = state_dim

        # 记忆库
        self.memory = np.zeros((memory_capacity, state_dim * 2 + action_dim + 1))
        self.pointer = 0
        # 定义 Actor 网络
        self.actor = Actor()
        self.actor_target = Actor()
        # 定义 Critic 网络
        self.critic = Critic()
        self.critic_target = Critic()
        # 定义优化器
        self.aopt = torch.optim.Adam(self.actor.parameters(), lr=lr_a)
        self.copt = torch.optim.Adam(self.critic.parameters(), lr=lr_c)
        # 选取损失函数
        self.mse_loss = nn.MSELoss()

    def sample(self):
        indices = np.random.choice(self.memory_capacity, size=self.batch_size)
        return self.memory[indices, :]

    def choose_action(self, state):
        state = torch.FloatTensor(state)
        action = self.actor(state)
        return action.detach().numpy()

    def learn(self):

        # soft replacement and hard replacement
        # 用于更新target网络的参数
        if self.replacement['name'] == 'soft':
            # soft的意思是每次learn的时候更新部分参数
            tau = self.replacement['tau']
            a_layers = self.actor_target.named_children()
            c_layers = self.critic_target.named_children()
            for al in a_layers:
                a = self.actor.state_dict()[al[0] + '.weight']
                al[1].weight.data.mul_((1 - tau))
                al[1].weight.data.add_(tau * self.actor.state_dict()[al[0] + '.weight'])
                al[1].bias.data.mul_((1 - tau))
                al[1].bias.data.add_(tau * self.actor.state_dict()[al[0] + '.bias'])
            for cl in c_layers:
                cl[1].weight.data.mul_((1 - tau))
                cl[1].weight.data.add_(tau * self.critic.state_dict()[cl[0] + '.weight'])
                cl[1].bias.data.mul_((1 - tau))
                cl[1].bias.data.add_(tau * self.critic.state_dict()[cl[0] + '.bias'])

        else:
            # hard的意思是每隔一定的步数才更新全部参数
            if self.t_replace_counter % self.replacement['rep_iter'] == 0:
                self.t_replace_counter = 0
                a_layers = self.actor_target.named_children()
                c_layers = self.critic_target.named_children()
                for al in a_layers:
                    al[1].weight.data = self.actor.state_dict()[al[0] + '.weight']
                    al[1].bias.data = self.actor.state_dict()[al[0] + '.bias']
                for cl in c_layers:
                    cl[1].weight.data = self.critic.state_dict()[cl[0] + '.weight']
                    cl[1].bias.data = self.critic.state_dict()[cl[0] + '.bias']

            self.t_replace_counter += 1

        # 从记忆库中采样bacth data
        bm = self.sample()
        bs = torch.FloatTensor(bm[:, :self.state_dim])
        ba = torch.FloatTensor(bm[:, self.state_dim:self.state_dim + self.action_dim])
        br = torch.FloatTensor(bm[:, -self.state_dim - 1: -self.state_dim])
        bs_ = torch.FloatTensor(bm[:, -self.state_dim:])

        # 训练Actor
        a = self.actor(bs)
        q = self.critic(bs, a)
        a_loss = -torch.mean(q)
        self.aopt.zero_grad()
        a_loss.backward(retain_graph=True)
        self.aopt.step()

        # 训练critic
        a_ = self.actor_target(bs_)
        q_ = self.critic_target(bs_, a_)
        q_target = br + self.gamma * q_
        q_eval = self.critic(bs, ba)
        td_error = self.mse_loss(q_target, q_eval)
        self.copt.zero_grad()
        td_error.backward()
        self.copt.step()

    def store_transition(self, s, a, r, s_):
        # a = np.ndarray(a)
        s = s.view(-1)
        s_ = s_.view(-1)
        a = a.reshape([3])

        transition = np.hstack((s, a, [r], s_))
        index = self.pointer % self.memory_capacity
        self.memory[index, :] = transition
        self.pointer += 1


if __name__ == '__main__':

    # hyper parameters
    VAR = 3  # control exploration
    MAX_EPISODES = 500
    MAX_EP_STEPS = 200
    MEMORY_CAPACITY = 10000
    REPLACEMENT = [
        dict(name='soft', tau=0.005),
        dict(name='hard', rep_iter=600)
    ][0]  # you can try different target replacement strategies

    RENDER = False

    # train
    env = enviroment.CircuitEnv(0, 0, 100)


    s_dim = env.state.shape[0]*env.state.shape[1]*env.state.shape[2]
    a_dim = 3
    # a_bound = env.action_space.high
    ddpg = DDPG(state_dim=s_dim,
                action_dim=a_dim,
                replacement=REPLACEMENT,
                memory_capacity=MEMORY_CAPACITY)

    t1 = time.time()
    for i in range(MAX_EPISODES):
        s = env.reset()
        ep_reward = 0
        for j in range(MAX_EP_STEPS):
            if RENDER:
                env.render()

            # Add exploration noise
            a = ddpg.choose_action(s)
            a = np.clip(np.random.normal(a, VAR), -2, 2)  # 在动作选择上添加随机噪声

            s_, r, done = env.step(a)

            ddpg.store_transition(s, a, r / 10, s_)

            if ddpg.pointer > MEMORY_CAPACITY:
                VAR *= .9995  # decay the action randomness
                ddpg.learn()

            s = s_
            ep_reward += r
            if j == MAX_EP_STEPS - 1:
                print('Episode:', i, ' Reward: %i' % int(ep_reward), 'Explore: %.2f' % VAR, )
                if ep_reward > -300: RENDER = True
                break

    print('Running time: ', time.time() - t1)
