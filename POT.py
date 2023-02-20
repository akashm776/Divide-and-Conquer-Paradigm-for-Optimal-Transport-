import numpy as np
import pandas as pd
import time
import ot

from decimal import Decimal
import numpy as np
with open('costMatrix') as file:
    cost = [[float(digit) for digit in line.split()] for line in file]
    


def getOT(n, useRandomPoints, filePath = ''):
    supply = np.full((1, n), 1).tolist()[0]
    demand = np.full((1, n), 1).tolist()[0]
    if useRandomPoints == False:
        df = pd.read_csv(filePath, delimiter=' ', names=['x', 'y'])
        p1 = []
        p2 = []
        for rows in df.itertuples():
            my_list = [rows.x, rows.y]
            if (len(p1) == n):
                if (len(p2) == n):
                    break
                p2.append(my_list)
            else:
                p1.append(my_list)
        p1 = np.array(p1)
        p2 = np.array(p2)
        costMatrix = ot.dist(p1, p2, 'euclidean')
    else:
        p1 = np.random.rand(n, 2)
        p2 = np.random.rand(n, 2)
        costMatrix = ot.dist(p1, p2, 'euclidean')
    if saveCostMatrix:
        np.savetxt('test.txt', costMatrix)
    start_time = time.time()
    print(ot.emd2(supply, demand, costMatrix, 1, 2000000000))
    end_time = time.time()
    print("--- %s seconds ---" % (end_time - start_time))

saveCostMatrix = False
N = 10000
months = [
    "jan_01_2014", "jan_24_2014", "mar_24_2014", "apr_13_2014", "jun_09_2014",
    "jul_17_2014", "aug_11_2014"
]
for month in months:
    print("--------------------Month: " + month + "--------------------")
    print()
    print("--------------------Pickup--------------------")
    for _ in range(5):
        getOT(N, False, "Datasets/10000/" + month + "_pickup.txt")
        print()
    print("--------------------Dropoff--------------------")
    for _ in range(5):
        getOT(N, False, "Datasets/10000/" + month + "_dropoff.txt")
        print()
    print("--------------------Pickup and Dropoff--------------------")
    for _ in range(5):
        getOT(N, False, "Datasets/10000/" + month + "_pick_and_drop.txt")
        print()