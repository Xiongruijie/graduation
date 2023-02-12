import re

from matplotlib import pyplot as plt, colors


def read_output_file(method, start_index, end_index):
    file_location = "./" + method + "/"
    time_list = []
    memory_list = []
    for i in range(start_index, end_index):
        file = open(file_location + "output_" + i.__str__() + ".out")
        file_lines = file.readlines()
        time = re.split(r'[ ]+', file_lines[-2])[-2]
        memory = re.split(r'[ ]+', file_lines[-1])[-2]
        time_list.append(time)
        memory_list.append(memory)
    time = [int(i) for i in time_list]
    memory = [int(i) for i in memory_list]

    return time, memory





def print_bar_line(data, x, ylabel_name, path="1.png"):
    # data：条形图数据
    # x:x轴坐标
    # path：图片保存路径

    # 创建颜色序列
    colors = []
    for _ in range(int(len(data) / 2)):
        colors.append([_ / int(len(data) / 2), 0.5, 0.5])
    colors = colors + colors[::-1]

    # 创建x轴显示的参数（此功能在与在图像中x轴仅显示能被10整除的刻度，避免刻度过多分不清楚）
    x_tick = list(map(lambda num: "" if num % 10 != 0 else num, x))

    # 创建一个分辨率为3000*1000的空白画布
    plt.figure(figsize=(300, 100), dpi=10)

    # 创建一个字体样式
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 300, }

    # 设置x轴的说明
    plt.xlabel('Conducts', font2)
    # 设置y轴的说明
    plt.ylabel(ylabel_name, font2)
    # 打开网格线
    plt.grid()
    # 绘制条形图
    plt.bar(range(len(data)), data, color=colors, width=0.5, label=ylabel_name)
    # 显示x轴刻度
    plt.xticks(range(len(x_tick)), x_tick, size=200)
    # 显示y轴刻度
    plt.yticks(size=200)
    # 获取当前图像句柄
    fig = plt.gcf()
    # plt.show()
    # 存储当前图像
    fig.savefig(path)


def paint_complex_pic():
    import numpy as np
    plt.rcParams['font.sans-serif'] = ['SimHei']
    # x = ("input_0", "input_1", "input_2", "input_3", "input_4", "input_5", "input_6")
    x = np.arange(7)
    gama = 0.001
    data_fdm_time = [9834594, 2670065, 4782398, 2671843, 3188282, 24190247, 34427844]
    data_fdm_memory = [27152, 34524, 56832, 45084, 58196, 192252, 291468]
    data_bem_time = [210862, 501418, 2172815, 883009, 2068257, 8386969, 24543641]
    data_bem_memory = [9360, 11188, 12468, 10992, 12184, 57484, 130560]
    data_frw_time = [2048983, 2169865, 1979932, 2149024, 1976674, 2832158, 3406090]
    data_frw_memory = [12380, 14036, 14316, 13184, 14280, 18912, 23440]
    bar_with = 0.3
    width = 0.2
    data_fdm_memory = [i * gama for i in data_fdm_memory]
    data_bem_memory = [i * gama for i in data_bem_memory]
    data_frw_memory = [i * gama for i in data_frw_memory]
    data_fdm_time = [i * gama for i in data_fdm_time]
    data_bem_time = [i * gama for i in data_bem_time]
    data_frw_time = [i * gama for i in data_frw_time]

    # plt.bar(x, data_fdm_memory,width=width,label='fdm')
    # plt.bar(x + bar_with, data_bem_memory,width=width,label='bem')
    # plt.bar(x + bar_with*2, data_frw_memory,width=width,label='frw')

    plt.bar(x, data_fdm_time, width=width, label='fdm')
    plt.bar(x + bar_with, data_bem_time, width=width, label='bem')
    plt.bar(x + bar_with * 2, data_frw_time, width=width, label='frw')

    plt.legend()
    fig = plt.gcf()

    fig.savefig('time pic.png')

    plt.show()

    # fig.savefig('time pic.png')


if __name__ == '__main__':
    # time, memory = read_output_file('fdm', 10, 99)
    # print_bar_line(time, list(range(len(time))), "fdm time", 'fdm_time_bar.png')
    # print_bar_line(memory, list(range(len(memory))), "fdm memory", 'fdm_memory_bar.png')
    # time, memory = read_output_file('bem', 10, 99)
    # print_bar_line(time, list(range(len(time))), "bem time", 'bem_time_bar.png')
    # print_bar_line(memory, list(range(len(memory))), "bem memory", 'bem_memory_bar.png')
    # time, memory = read_output_file('frw', 10, 99)
    # print_bar_line(time, list(range(len(time))), "frw time", 'frw_time_bar.png')
    # print_bar_line(memory, list(range(len(memory))), "frw memory", 'frw_memory_bar.png')
    paint_complex_pic()
