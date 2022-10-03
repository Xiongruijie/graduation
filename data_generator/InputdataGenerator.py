import os
import random


# 判定是否存在文件
# 返回顺位第一个不存在的数据索引
def move_index(file_index=0):
    while (os.path.exists("../Fieldsolver2d_hybrid_to_improve/input/input_" + file_index.__str__() + ".data")):
        print("input_" + file_index.__str__() + ".data")
        file_index += 1

    return file_index


# 从初始的index 开始生成数据
def generate_files(start_index, generate_number):
    while (generate_number > 0):
        write_data(start_index)
        generate_number -= 1


class Square(object):
    def __init__(self, width=None, hight=None, left_low=None, net_number=0):
        self.net_number = net_number
        self.width = width
        self.hight = hight
        self.left_low = left_low
        self.right_high = [self.left_low[0] + self.width, self.left_low[1] + self.hight]

    def overlap(self, rectangle):
        min_x = max(self.left_low[0], rectangle.left_low[0])
        min_y = max(self.left_low[1], rectangle.left_low[1])
        max_x = min(self.right_high[0], rectangle.right_high[0])
        max_y = min(self.right_high[1], rectangle.right_high[1])

        if (min_x > max_x or min_y > max_y):
            return False
        else:
            return True

    def inspect_data_line(self, clip):
        if self.right_high[0] > clip.right_high[0] :

            return False
        else:
            return True
    def inspect_data_column(self, clip):
        if self.right_high[1] > clip.right_high[1]:
            return False
        else:
            return True

    def output(self):
        rec0 = ('%.4f' % self.left_low[0]).__str__() + "  " + (
                '%.4f' % self.left_low[1]).__str__() + "    " + \
               ('%.4f' % self.right_high[0]).__str__() + "  " + (
                       '%.4f' % self.right_high[1]).__str__()

        return rec0


class Trapezoid(Square):
    def __init__(self, width=None, hight=None, left_low=None, net_number=None):
        if left_low is None:
            left_low = [random.uniform(-9.9, 9.8), random.uniform(0, 9.8)]
            print("net" + net_number.__str__() + "初始值为空，net" + net_number.__str__() + "的位置被随机赋值")

        super().__init__(width=width, hight=hight, left_low=left_low)
        width_diff = self.width * 0.1
        self.width_diff = width_diff
        self.hight_diff = self.hight / 4
        self.rectangle_0_left_low = self.left_low
        self.rectangle_0_right_high = [self.right_high[0], self.left_low[1] + self.hight_diff]
        self.rectangle_1_left_low = [self.rectangle_0_left_low[0] + self.width_diff, self.rectangle_0_right_high[1]]
        self.rectangle_1_right_high = [self.rectangle_0_right_high[0] - self.width_diff,
                                       self.rectangle_1_left_low[1] + self.hight_diff]
        self.rectangle_2_left_low = [self.rectangle_1_left_low[0] + self.width_diff, self.rectangle_1_right_high[1]]
        self.rectangle_2_right_high = [self.rectangle_1_right_high[0] - self.width_diff,
                                       self.rectangle_2_left_low[1] + self.hight_diff]
        self.rectangle_3_left_low = [self.rectangle_2_left_low[0] + self.width_diff, self.rectangle_2_right_high[1]]
        self.rectangle_3_right_high = [self.rectangle_2_right_high[0] - self.width_diff,
                                       self.rectangle_3_left_low[1] + self.hight_diff]

    def output(self):
        rec0 = ('%.4f' % self.rectangle_0_left_low[0]).__str__() + "  " + (
                '%.4f' % self.rectangle_0_left_low[1]).__str__() + "    " + \
               ('%.4f' % self.rectangle_0_right_high[0]).__str__() + "  " + (
                       '%.4f' % self.rectangle_0_right_high[1]).__str__()
        rec1 = ('%.4f' % self.rectangle_1_left_low[0]).__str__() + "  " + (
                '%.4f' % self.rectangle_1_left_low[1]).__str__() + "    " + \
               ('%.4f' % self.rectangle_1_right_high[0]).__str__() + "  " + (
                       '%.4f' % self.rectangle_1_right_high[1]).__str__()
        rec2 = ('%.4f' % self.rectangle_2_left_low[0]).__str__() + "  " + (
                '%.4f' % self.rectangle_2_left_low[1]).__str__() + "    " + \
               ('%.4f' % self.rectangle_2_right_high[0]).__str__() + "  " + (
                       '%.4f' % self.rectangle_2_right_high[1]).__str__()
        rec3 = ('%.4f' % self.rectangle_3_left_low[0]).__str__() + "  " + (
                '%.4f' % self.rectangle_3_left_low[1]).__str__() + "    " + \
               ('%.4f' % self.rectangle_3_right_high[0]).__str__() + "  " + (
                       '%.4f' % self.rectangle_3_right_high[1]).__str__()
        return [rec0, rec1, rec2, rec3]


class Clip(object):
    def __init__(self, x1, y1, x2, y2, ClipList):
        self.left_low = [x1, y1]
        self.right_high = [x2, y2]
        self.flag = True
        for i in range(len(ClipList)):
            if self.overlap(ClipList[i]) == False:
                self.flag = False
                break
            else:
                self.flag = True
        self.fall_flag = False

    def overlap(self, clip):
        min_x = max(self.left_low[0], clip.left_low[0])
        min_y = max(self.left_low[1], clip.left_low[1])
        max_x = min(self.right_high[0], clip.right_high[0])
        max_y = min(self.right_high[1], clip.right_high[1])
        if (min_x > max_x or min_y > max_y):
            return True
        else:
            return False


def print_Squares_Trapezoids(net_count):
    # 随机设置初始点

    interval_x = random.uniform(0.12, 0.15).__round__(4)
    interval_y = random.uniform(0.12, 0.15).__round__(4)
    square_width = random.uniform(0.03, 0.05).__round__(4)
    square_hight = random.uniform(0.05, 0.1).__round__(4)
    master_square_width = random.uniform(0.5, 2).__round__(4)
    master_square_hight = random.uniform(0.5, 2).__round__(4)

    ConductorList = []
    ClipList = []
    # 创建master导体
    point = [random.uniform(-10, 7.9).__round__(4), random.uniform(0, 7.9).__round__(4)]
    clip = Clip(point[0], point[1], point[0] + 2, point[1] + 2, ClipList)
    ClipList.append(clip)
    is_trapezoid = random.randint(0, 1)
    if is_trapezoid == 0:
        # 构造正方形导体
        conductor = Square(master_square_width, master_square_hight, point, len(ConductorList))
    else:
        # 构造梯形导体
        conductor = Trapezoid(master_square_width, master_square_hight, point, len(ConductorList))
    ConductorList.append(conductor)


    # 创建net_count个矩形或者梯形
    while (len(ConductorList) < net_count - 1):
        # 生成3*3的线簇片
        while True:
            point = [random.uniform(-10, 7.9).__round__(4), random.uniform(0, 7.9).__round__(4)]
            clip = Clip(point[0], point[1], point[0] + 2, point[1] + 2, ClipList)
            if clip.flag is True:
                break
        ClipList.append(clip)

        line_flag = True
        column_flag = True
        while (line_flag):
            is_trapezoid = random.randint(0, 1)
            if is_trapezoid == 0:
                # 构造正方形导体
                conductor = Square(square_width, square_hight, point, len(ConductorList))
                line_flag = conductor.inspect_data_line(clip)
                column_flag = conductor.inspect_data_column(clip)
            else:
                # 构造梯形导体
                conductor = Trapezoid(square_width, square_hight, point, len(ConductorList))
                line_flag = conductor.inspect_data_line(clip)
                column_flag = conductor.inspect_data_column(clip)

            if line_flag == False:
                point = [clip.left_low[0], point[1] + square_hight + interval_y]
                line_flag = True
            elif column_flag == False:
                break
            else:
                point = [point[0] + interval_x, point[1]]

            ConductorList.append(conductor)
            if len(ConductorList) >= net_count:
                break



    return ConductorList


def write_data(file_index, maxNumber):
    # 判定是否该文件是否存在 调试时记得去除
    # if (os.path.exists("../Fieldsolver2d_hybrid_to_improve/input/input_" + file_index.__str__() + ".data")):
    #     print(("input_" + file_index.__str__() + ".data" + "文件已经存在"))
    #     return
    file_path = "../Fieldsolver2d_hybrid_to_improve/input/input_" + file_index.__str__() + ".data"
    file_name = "input_" + file_index.__str__() + ".data"
    net_count = random.randint(10, maxNumber)
    data_file = open(file_path, mode='w')
    print("创建文件" + file_name)
    data_file.write("boundary  -10 0  10 9.9\n")
    dielectric = random.uniform(3.0, 9.9).__round__(1)
    data_file.write("dielectric  " + dielectric.__str__() + "\n")
    # todo 创建data并写进文件
    ConductorList = print_Squares_Trapezoids(net_count)

    for i in range(len(ConductorList)):
        if isinstance(ConductorList[i], Trapezoid):
            data_file.write("net net" + i.__str__() + "    " + ConductorList[i].output()[0] + "\n")
            data_file.write("net net" + i.__str__() + "    " + ConductorList[i].output()[1] + "\n")
            data_file.write("net net" + i.__str__() + "    " + ConductorList[i].output()[2] + "\n")
            data_file.write("net net" + i.__str__() + "    " + ConductorList[i].output()[3] + "\n")
        else:
            data_file.write("net net" + i.__str__() + "    " + ConductorList[i].output() + "\n")
    data_file.close()


if __name__ == '__main__':
    for i in range(9,3000):
        write_data(i,100)
