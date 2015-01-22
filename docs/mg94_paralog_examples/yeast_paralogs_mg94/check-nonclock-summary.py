import pandas as pd
import numpy as np

def main():
    data_filename = 'MG94_nonclock_summary.txt'
    with open(data_filename) as fin:
        lines = fin.readlines()
    column_names = lines[0][1:].split()
    row_names = lines[-1][1:].split()
    data = np.array([[float(x) for x in line.split()] for line in lines[1:-1]])
    #print(len(row_names))
    #print(len(column_names))
    #print(data.shape)
    print('row names:')
    print(row_names)
    #print('column names:')
    #print(column_names)
    #print('data:')
    #print(data)
    F = pd.DataFrame(data, columns=column_names) #  rows=row_names)
    print(F)

main()
