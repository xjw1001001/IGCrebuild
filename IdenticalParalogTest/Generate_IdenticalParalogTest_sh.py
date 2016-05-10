if __name__ == '__main__':

    pairs = []
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    
    model = 'MG94_'
    sh_line = 'sbatch -o IPT-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/'

    #pairs = [pairs[0]]
    with open('./Run_IPT.sh', 'w+') as f:
        f.write('#!/bin/bash' + '\n')
        for paralog in pairs:
            f.write(sh_line + model + '_'.join(paralog) + '_IPT.sh\n')
            with open('./ShFiles/' + model + '_'.join(paralog) + '_IPT.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                printscreen = ' > ' + '_'.join(paralog) + '_IPT_PrintScreen.txt\n'
                g.write('python Run_IdenticalParalogTest.py ' + ' --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --no-force --model MG94 --no-clock ' + printscreen)


