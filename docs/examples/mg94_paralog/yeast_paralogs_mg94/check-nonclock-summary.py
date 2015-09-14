import pandas as pd

def main():
    data_filename = 'MG94_nonclock_summary.txt'
    with open(data_filename) as fin:
        lines = fin.readlines()
    # Column names are on the first line, and the first character is #.
    column_names = lines[0][1:].split()
    # Row names are on the last line, and the first character is #.
    row_names = lines[-1][1:].split()
    data = [[float(x) for x in line.split()] for line in lines[1:-1]]
    F = pd.DataFrame(data, columns=column_names, index=row_names)
    print(F)

main()
