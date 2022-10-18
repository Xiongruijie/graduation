#!/usr/bin/env python
# coding: utf-8

# In[1]:


import parl
import paddle
import paddle.nn as nn
import paddle.nn.functional as F
import numpy as np
import enviroment
from paddle.distribution import Normal
from parl.utils import logger
import argparse
from itertools import count
from parl.utils import logger, summary, ReplayMemory
from collections import deque
import random
from visualdl import LogWriter


# In[2]:




class Actor(parl.Model):
    def __init__(self, action_dim=3, is_train=True):
        super(Actor, self).__init__()
        self.action_dim = action_dim
        self.is_train = is_train
        self.conv = nn.Sequential(
            nn.Conv2D(3, 96, 11, 4),
            nn.ReLU(),
            nn.MaxPool2D(3, 3),
            nn.Conv2D(96, 256, 3, 1, 1),
            nn.MaxPool2D(3, 3),
            nn.ReLU(),
            nn.Conv2D(256, 64, 3, 1, 1),
            nn.ReLU(),
            nn.MaxPool2D(3, 2),
        )
        self.fc = nn.Sequential(
            # nn.Linear(20800, 4096),#64 * 9 * 18
            nn.Linear(10368,4096),
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
        # obs = obs.unsqueeze(0)
        feature = self.conv(obs)
        feature = feature.reshape([obs.shape[0], -1])
        output = self.fc(feature)
        return output

    def select_action(self, epsilon, state):
        state = paddle.to_tensor(state, dtype="float32")
        with paddle.no_grad():
            action = self.forward(state)
            action = action + self.is_train * epsilon * self.noisy.sample(action.shape)
        return action.numpy()
        # return action.shape


# In[3]:



class Critic(parl.Model):
    def __init__(self, action_dim=3, ):
        super(Critic, self).__init__()
        self.action_dim = action_dim
        self.conv = nn.Sequential(
            nn.Conv2D(3, 96, 11, 4),
            nn.MaxPool2D(3, 3),
            nn.Conv2D(96, 256, 3, 1, 1),
            nn.MaxPool2D(3, 3),
            nn.Tanh(),
            nn.Conv2D(256, 64, 3, 1, 1),
            nn.Tanh(),
            nn.MaxPool2D(3, 2),
        )
        # 因为中间要加动作判定所以需要fc分层
        self.fc = nn.Sequential(
            nn.Linear(64 * 9 * 18, 4096),
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
        feature = self.conv(obs)
        feature = feature.reshape([obs.shape[0], -1])
        output = self.fc(feature)

        output = paddle.concat((output, action), axis=1)
        output = self.fc2(output)

        return output


# In[4]:


# class MujocoAgent(parl.Agent):
#     def __init__(self, algorithm, act_dim, expl_noise=0.1):
#         assert isinstance(act_dim, int)
#         super(MujocoAgent, self).__init__(algorithm)
#         self.act_dim = act_dim
#         self.expl_noise = expl_noise
#         self.alg.sync_target(decay=0)
#
#     def sample(self, obs):
#         action_numpy = self.predict(obs)
#         action_noise = np.random.normal(0, self.expl_noise, size=self.act_dim)
#         action = (action_numpy + action_noise).clip(-1, 1)
#         return action
#
#     def predict(self, obs):
#         obs = paddle.to_tensor(obs.reshape(1, -1), dtype='float32')
#         action = self.alg.predict(obs)
#         action_numpy = action.cpu().numpy()[0]
#         return action_numpy
#
#     def learn(self, obs, action, reward, next_obs, terminal):
#         terminal = np.expand_dims(terminal, -1)
#         reward = np.expand_dims(reward, -1)
#
#         obs = paddle.to_tensor(obs, dtype='float32')
#         action = paddle.to_tensor(action, dtype='float32')
#         reward = paddle.to_tensor(reward, dtype='float32')
#         next_obs = paddle.to_tensor(next_obs, dtype='float32')
#         terminal = paddle.to_tensor(terminal, dtype='float32')
#         critic_loss, actor_loss = self.alg.learn(obs, action, reward, next_obs,
#                                                  terminal)
#         return critic_loss, actor_loss


# class MujocoModel(parl.Model):
#     def __init__(self, action_dim):
#         super(MujocoModel, self).__init__()
#         self.actor_model = Actor(action_dim)
#         self.critic_model = Critic(action_dim)
#
#     def policy(self, obs):
#         return self.actor_model(obs)
#
#     def value(self, obs, action):
#         return self.critic_model(obs, action)
#
#     def get_actor_params(self):
#         return self.actor_model.parameters()
#
#     def get_critic_params(self):
#         return self.critic_model.parameters()


# In[5]:


class Memory(object):
    def __init__(self, memory_size) -> None:
        self.memory_size = memory_size
        self.buffer = deque(maxlen=self.memory_size)

    def add(self,experience)->None:
        self.buffer.append(experience)

    def size(self):
        return len(self.buffer)

    def sample(self, batch_size: int, continuous: bool = False):
        if batch_size > len(self.buffer):
            batch_size = len(self.buffer)
        if continuous:
            rand = random.randint(0, len(self.buffer) - batch_size)
            return [self.buffer[i] for i in range(rand, rand + batch_size)]
        else:
            indexes = np.random.choice(np.arange(len(self.buffer)), size=batch_size, replace=False)
            return [self.buffer[i] for i in indexes]

    def clear(self):
        self.buffer.clear()


# In[6]:


def soft_update(target, source, tau):
    for target_param, param in zip(target.parameters(), source.parameters()):
        target_param.set_value( target_param * (1.0 - tau) + param * tau)


# In[7]:







def main():
    logger.info("------------------ DDPG ---------------------")
    logger.info('Env: {}, Seed: {}'.format(args.env, args.seed))
    logger.info("---------------------------------------------")
    logger.set_dir('./{}_{}'.format(args.env, args.seed))
    paddle.device.set_device('cpu')
    print(paddle.device.get_device())

    # 定义超参赛
    explore = 50000
    epsilon = 1
    gamma = 0.99
    tau = 0.001

    memory_replay = Memory(50000)
    begin_train = False
    batch_size = 1

    epochs = 250
    learn_steps = 0
    writer = LogWriter('logs')

    actor = Actor()
    critic = Critic()
    actor_target = Actor()
    critic_target = Critic()
    env = enviroment.CircuitEnv(0, 0, 99)

    # 优化器
    critic_optim = paddle.optimizer.Adam(parameters=critic.parameters(), learning_rate=3e-5)
    actor_optim = paddle.optimizer.Adam(parameters=actor.parameters(), learning_rate=1e-5)





    for epoch in range(0,epochs): #count():
        state = env.reset()
        episode_reward = 0
        for time_step in range(10):#99
            action = actor.select_action(epsilon, state)
            next_state, reward, done = env.step(action)
            episode_reward = reward
            # reward = (reward + 8.1) / 8.1
            memory_replay.add((state, action, reward))
            if memory_replay.size() > 10:#1280
                learn_steps += 1
                if not begin_train:
                    print('train begin!')
                    begin_train = True
                experiences = memory_replay.sample(batch_size, False)
                batch_state,  batch_action, batch_reward = zip(*experiences)
                batch_state = paddle.to_tensor(batch_state, dtype="float32").squeeze(0)
                batch_action = paddle.to_tensor(batch_action, dtype="float32").squeeze(1)
                batch_reward = paddle.to_tensor(batch_reward, dtype="float32").unsqueeze(1)

                # 均方误差 y - Q(s, a) ， y是目标网络所看到的预期收益， 而 Q(s, a)是Critic网络预测的操作值。
                # y是一个移动的目标，评论者模型试图实现的目标；这个目标通过缓慢的更新目标模型来保持稳定。
                with paddle.no_grad():
                    Q_next = critic_target(batch_state, batch_action)
                    Q_target = batch_reward + gamma * Q_next

                critic_loss = F.mse_loss(critic(batch_state, batch_action), Q_target)

                critic_optim.clear_grad()
                critic_loss.backward()
                critic_optim.step()

                writer.add_scalar('critic loss', critic_loss.numpy(), learn_steps)
                # 使用Critic网络给定值的平均值来评价Actor网络采取的行动。 我们力求使这一数值最大化。
                # 因此，我们更新了Actor网络，对于一个给定状态，它产生的动作尽量让Critic网络给出高的评分。
                critic.eval()
                actor_loss = - critic(batch_state, actor(batch_state))
                actor_loss = actor_loss.mean()
                actor_optim.clear_grad()
                actor_loss.backward()
                actor_optim.step()
                critic.train()
                writer.add_scalar('actor loss', actor_loss.numpy(), learn_steps)
                soft_update(actor_target, actor, tau)
                soft_update(critic_target, critic, tau)


            if epsilon > 0:
                epsilon -= 1 / explore
            state = next_state

        writer.add_scalar('episode reward', episode_reward, epoch)
        if epoch % 10 == 0:
            print('Epoch:{}, episode reward is {}'.format(epoch, episode_reward))

        # if epoch % 200 == 0:
        if epoch % 20 == 0:
            paddle.save(actor.state_dict(), 'model/ddpg-actor' + str(epoch) + '.para')
            paddle.save(critic.state_dict(), 'model/ddpg-critic' + str(epoch) + '.para')
            print('model saved!')





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--env", default="CircuitEnv", help='CircuitEnv environment name')
    parser.add_argument("--seed", default=0, type=int, help='Sets Gym seed')
    # parser.add_argument(
    #     "--train_total_steps",
    #     default=5e6,
    #     type=int,
    #     help='Max time steps to run environment')
    # parser.add_argument(
    #     '--test_every_steps',
    #     type=int,
    #     default=int(5e3),
    #     help='The step interval between two consecutive evaluations')
    args = parser.parse_args(args=[])
    main()
    print('end')


# In[ ]:




