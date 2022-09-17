from unittest import TestCase
from InputdataGenerator import Trapezoid, Square, print_Squares_Trapezoids


class TestTrapezoid(TestCase):
    def test_case_1(self):
        net7 = Trapezoid(0.0318, 0.0666, [0.1761, 0.0723])
        print(net7.rectangle_0_left_low.__str__() + "  " + net7.rectangle_0_right_high.__str__())
        print(net7.rectangle_1_left_low.__str__() + "  " + net7.rectangle_1_right_high.__str__())
        print(net7.rectangle_2_left_low.__str__() + "  " + net7.rectangle_2_right_high.__str__())
        print(net7.rectangle_3_left_low.__str__() + "  " + net7.rectangle_3_right_high.__str__())

        net71 = Square(0.044, 0.027, [1.3730, 0.2100])
        print(net71.left_low[0].__round__(4).__str__() + " " + net71.left_low[1].__round__(
            4).__str__() + "  " + net71.right_high.__str__())

    def test_case_2(self):
        a = 0.1
        a = ('%.4f' % a)
        print(a)


class Test(TestCase):
    pass


class Test(TestCase):
    def test_print_squares_trapezoids(self):

        ConductorList = print_Squares_Trapezoids(110)
        for i in range(len(ConductorList)):
            if (isinstance(ConductorList[i], Trapezoid)):
                print("net net" + i.__str__() + "   " + ConductorList[i].output()[0])
                print("net net" + i.__str__() + "   " + ConductorList[i].output()[1])
                print("net net" + i.__str__() + "   " + ConductorList[i].output()[2])
                print("net net" + i.__str__() + "   " + ConductorList[i].output()[3])

            else:
                print("net net" + i.__str__() + "   " + ConductorList[i].output())
