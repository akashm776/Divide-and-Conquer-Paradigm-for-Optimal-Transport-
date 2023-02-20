import ot
import pandas as pd

from decimal import Decimal
import numpy as np


with open("out.txt", "r") as file:
    costLine = file.readline().strip('\n')
    cost = [float(digit) for digit in costLine.split()]
    cost = np.array(cost)
    supplyLine = file.readline().strip('\n')
    supply = [float(digit) for digit in supplyLine.split()]
    supply = np.array(supply)
    #print(supply)
    demandLine = file.readline().strip('\n')
    demand = [float(digit) for digit in demandLine.split()]
    demand = np.array(demand)
    cost.resize((supply.size, demand.size))
    #print(cost)
    print(ot.emd2(supply, demand, cost, 1, 2000000000))

