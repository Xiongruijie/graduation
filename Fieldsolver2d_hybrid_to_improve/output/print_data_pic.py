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
    memory=[int(i) for i in memory_list]
    return time, memory


def print_broken_line(data):
    pass


def print_bar_line(data, x, path="1.png"):
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
    plt.ylabel('Time Consume', font2)
    # 打开网格线
    plt.grid()
    # 绘制条形图
    plt.bar(range(len(data)), data, color=colors, width=0.5)
    # 显示x轴刻度
    plt.xticks(range(len(x_tick)), x_tick, size=200)
    # 显示y轴刻度
    plt.yticks(size=200)
    # 获取当前图像句柄
    fig = plt.gcf()
    # plt.show()
    # 存储当前图像
    fig.savefig(path)



if __name__ == '__main__':
    time, memory = read_output_file('fdm', 10, 99)
    print_bar_line(time, list(range(len(time))),'fdm_time_bar.png')
    print_bar_line(memory, list(range(len(memory))), 'fdm_memory_bar.png')
    time, memory = read_output_file('bem', 10, 99)
    print_bar_line(time, list(range(len(time))), 'bem_time_bar.png')
    print_bar_line(memory, list(range(len(memory))), 'bem_memory_bar.png')
    time, memory = read_output_file('frw', 10, 99)
    print_bar_line(time, list(range(len(time))), 'frw_time_bar.png')
    print_bar_line(memory, list(range(len(memory))), 'frw_memory_bar.png')
    # time

