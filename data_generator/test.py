import random


def overlap(rec1, rec2):
    min_x = max(rec1.left_low_x, rec2.left_low_x)
    min_y = max(rec1.left_low_y, rec2.left_low_y)
    max_x = min(rec1.right_high_x, rec2.right_high_x)
    max_y = min(rec1.right_high_y, rec2.right_high_y)

    if (min_x > max_x or min_y > max_y):
        return False;
    else:
        return True;


class Rectangle(object):
    def __init__(self, left_low_x, left_low_y, right_high_x, right_high_y):
        self.left_low_x = left_low_x
        self.left_low_y = left_low_y
        self.right_high_x = right_high_x
        self.right_high_y = right_high_y


if __name__ == '__main__':
    # rec_1 = Rectangle(0, 0, 5, 4)
    # rec_2 = Rectangle(3, -1, 6, 2)
    #
    # rec_3 = Rectangle(-3, 2, 1, 4)
    # print(overlap(rec_1, rec_2))
    # print(overlap(rec_2, rec_3))
    for i in range(10):
        print(i)
