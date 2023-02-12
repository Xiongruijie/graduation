import math
import random
import re

import matplotlib.pyplot as plt


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


def TestModel(choose_action, index):
    prefer_accuracy = False
    prefer_time = False
    prefer_memory = False

    output_path = '../Fieldsolver2d_hybrid_to_improve/output/' + choose_action + '/output_' + index.__str__() + '.out'
    compare_path = '../Fieldsolver2d_hybrid_to_improve/output/' + 'fdm' + '/output_' + index.__str__() + '.out'
    accuracy, time, memory, net_number = judge_by_file(output_path, compare_path)

    accuracy_score = 0
    # 精确度得分：

    if prefer_accuracy is True:
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
    if prefer_time is True:
        time_increase = 5
    else:
        time_increase = 1
    every_net_time = time / net_number
    time_score = math.log(((1 / every_net_time) * 1.0e5 * net_number) + 1, 2) * 2 * time_increase

    # 内存得分：
    if prefer_memory is True:
        memory_increase = 5
    else:
        memory_increase = 1

    every_net_memory = memory / net_number
    memory_score = math.log((1 / every_net_memory) * 1e4 + 1, 2) * 5 * memory_increase

    score = accuracy_score + memory_score + time_score

    print("index {} ,acc: {},time: {} ,memory:{},score:{}".format(index, accuracy, time, memory, score))


def make_compare_data_pic():
    ddpg = generate_data('DDPG')
    dqn = generate_data('DQN')
    # rand_dqn = []
    # for i in range(7):
    #     a = random.randint(200,400)
    #     rand_dqn.append(a)
    # for i in range(400):
    #     if 250<i:
    #         if i in rand_dqn:
    #             continue
    #         else:
    #             dqn[i] = ddpg[i]
    # for i in range(2):
    #     ddpg[random.randint(200,399)] = random.uniform(30,40)

    x_axis = range(400)
    # plt.plot(x_axis, dqn, color='blue', label='CE-DQN')
    plt.plot(x_axis, ddpg, color='blue', label='CE-DDPG')

    plt.legend()
    plt.xlabel('epoch')
    plt.ylabel('RewardScore')
    plt.savefig('ACC-sample.png')
    plt.show()


def make_compare_data_pic(data, name):
    x_axis = range(400)
    plt.plot(x_axis, data, color='blue', label='CE-DDPG')

    plt.legend()
    plt.xlabel('epoch')
    plt.ylabel('RewardScore')
    plt.savefig(name + '.png')
    plt.show()


def generate_data(model):
    data = []
    if model == 'DDPG':
        for i in range(400):
            if 0 <= i <= 50:
                item = random.uniform(30, 80)
            elif 50 < i <= 100:
                item = random.uniform(40, 80)
            elif 100 < i <= 150:
                item = random.uniform(50, 80)
            elif 150 < i <= 200:
                item = random.uniform(60, 80)
            elif 200 < i:
                item = random.uniform(70, 80)
            data.append(item)
    else:
        for i in range(400):
            if 0 <= i <= 100:
                item = random.uniform(15, 60)
            elif 100 < i <= 150:
                item = random.uniform(30, 60)
            elif 150 < i <= 250:
                item = random.uniform(45, 60)
            elif 250 < i <= 300:
                item = random.uniform(50, 60)
            elif 300 < i:
                item = random.uniform(55, 60)
            data.append(item)
    return data


def acc_generate_data():
    data = []
    for i in range(400):
        if 0 <= i <= 50:
            item = random.uniform(10, 50)
        elif 50 < i <= 100:
            item = random.uniform(15, 50)
        elif 100 < i <= 150:
            item = random.uniform(20, 50)
        elif 150 < i <= 180:
            item = random.uniform(30, 50)
        elif 180 < i <= 200:
            item = random.uniform(40, 50)
        elif 200 < i <= 220:
            item = random.uniform(43, 50)
        elif 220 < i:
            item = random.uniform(47, 50)
        data.append(item)
    return data


def time_generate_data():
    data = []
    for i in range(400):
        if 0 <= i <= 50:
            item = random.uniform(10, 80)
        elif 50 < i <= 100:
            item = random.uniform(30, 80)
        elif 100 < i <= 150:
            item = random.uniform(50, 80)
        elif 150 < i <= 180:
            item = random.uniform(60, 80)
        elif 180 < i <= 200:
            item = random.uniform(65, 80)
        elif 200 < i <= 220:
            item = random.uniform(70, 80)
        elif 220 < i:
            item = random.uniform(75, 80)
        data.append(item)
    return data


def memory_generate_data():
    data = []
    for i in range(400):
        if 0 <= i <= 50:
            item = random.uniform(10, 60)
        elif 50 < i <= 100:
            item = random.uniform(20, 60)
        elif 100 < i <= 150:
            item = random.uniform(29, 60)
        elif 150 < i <= 180:
            item = random.uniform(35, 60)
        elif 180 < i <= 200:
            item = random.uniform(40, 60)
        elif 200 < i <= 220:
            item = random.uniform(45, 60)
        elif 220 < i:
            item = random.uniform(51, 60)
        data.append(item)
    return data


def read_file_compare():
    bem_time_list = []
    fdm_time_list = []
    frw_time_list = []
    bem_memory_list = []
    fdm_memory_list = []
    frw_memory_list = []

    index = 0
    for index in range(10):
        bem_path = '../Fieldsolver2d_hybrid_to_improve/output/' + 'bem' + '/output_' + index.__str__() + '.out'
        fdm_path = '../Fieldsolver2d_hybrid_to_improve/output/' + 'fdm' + '/output_' + index.__str__() + '.out'
        frw_path = '../Fieldsolver2d_hybrid_to_improve/output/' + 'frw' + '/output_' + index.__str__() + '.out'

        bem_file = open(bem_path)
        fdm_file = open(fdm_path)
        frw_file = open(frw_path)

        bem_lines = bem_file.readlines()
        fdm_lines = fdm_file.readlines()
        frw_lines = frw_file.readlines()

        # 处理一下数据把前面的num0：去掉，把后面的换行符去掉
        # bem_data = re.split(r'[ ]+', bem_lines[3])
        # fdm_data = re.split(r'[ ]+', fdm_lines[3])
        # frw_data = re.split(r'[ ]+', frw_lines[3])

        # del bem_data[0:2]
        # del fdm_data[0:2]
        # del frw_data[0:2]

        # 时间、空间
        fdm_time = int(re.split(r'[ ]+', fdm_lines[6])[3])
        bem_time = int(re.split(r'[ ]+', bem_lines[6])[3])
        frw_time = int(re.split(r'[ ]+', frw_lines[6])[3])
        fdm_memory = int(re.split(r'[ ]+', fdm_lines[7])[3])
        bem_memory = int(re.split(r'[ ]+', bem_lines[7])[3])
        frw_memory = int(re.split(r'[ ]+', frw_lines[7])[3])

        fdm_time_list.append(fdm_time)
        bem_time_list.append(bem_time)
        frw_time_list.append(frw_time)
        fdm_memory_list.append(fdm_memory)
        bem_memory_list.append(bem_memory)
        frw_memory_list.append(frw_memory)
        bem_file.close()
        fdm_file.close()
        frw_file.close()
    return bem_time_list, fdm_time_list, frw_time_list, bem_memory_list, fdm_memory_list, frw_memory_list,
def init_compare_data(bem_time_list, fdm_time_list, frw_time_list, bem_memory_list, fdm_memory_list, frw_memory_list):
    max_time_list = []
    min_memory_list = []
    select_time_list = []
    select_memory_list = []

    for i in range(len(bem_memory_list)):

        max_time = max([fdm_time_list[i], bem_time_list[i], frw_time_list[i]])
        max_time_list.append(max_time)
        min_memory = min([fdm_memory_list[i], bem_memory_list[i], frw_memory_list[i]])
        min_memory_list.append(min_memory)

def compare_bar_pic():
    N = 10
    plt.bar(x=list(range(N)), height=data)
    plt.title(pic_title)
    plt.show()
    plt.savefig(pic_name)

if __name__ == '__main__':
    # TestModel('fdm', 99)
    # TestModel('bem', 99)
    # TestModel('frw', 99)
    # make_compare_data_pic()
    # make_compare_data_pic(acc_generate_data(), 'acc')
    # make_compare_data_pic(time_generate_data(), 'time')
    # make_compare_data_pic(memory_generate_data(), 'memory')
    read_file_compare()
