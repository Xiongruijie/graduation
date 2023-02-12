import csv
import matplotlib.pyplot as plt
import numpy as np
def pain(filename,mode):
    f = open(filename,'r')
    data = f.readline()
    data = data.split(',')
    data.remove('')
    for i in range(len(data)):
        data[i] = float(data[i])
    print(data)

    reward = np.array(data)
    plt.plot(reward)
    plt.show()
    plt.savefig('./pic/'+mode+'brokenline.jpg')



if __name__ == '__main__':
    mode = 'acc'
    pain('rewardListFile.csv',mode)