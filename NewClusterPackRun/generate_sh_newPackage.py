if __name__ == '__main__':
    pairs = []
    with open('../All_Pairs.txt', 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    pairs.remove(['YLR028C', 'YMR120C'])
    Force_str = '{5:0.0}'

    sh_line = 'sbatch -o cd-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/'
    
    with open('./CdNPack.sh', 'w+') as f:
        f.write('#!/bin/bash' + '\n')
        for pair in pairs:
            save_txt = ' > Cd_' + '_'.join(pair) + '_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + 'cd.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + 'cd.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python CodonGeneconv.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + save_txt + '\n')
            
    with open('./CdNPack_ForceTau.sh', 'w+') as f:
        f.write('#!/bin/bash' + '\n')
        for pair in pairs:
            save_txt = ' > ForceCd_' + '_'.join(pair) + '_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + 'cd_ForceTau.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + 'cd_ForceTau.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python CodonGeneconv.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --Force '+ Force_str + save_txt + '\n')
