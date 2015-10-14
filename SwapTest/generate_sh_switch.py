if __name__ == '__main__':
    
    pairs = []
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    if ['YLR028C', 'YMR120C'] in pairs:
        pairs.remove(['YLR028C', 'YMR120C'])

    sh_line = 'sbatch -o cd-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/'
    ctl_lines = [' --model HKY --clock --force --no-dir --no-gBGC --switch',
                 ' --model HKY --clock --no-force --no-dir --no-gBGC --switch']
                 #' --model HKY --clock --no-force --no-dir --gBGC --switch',
                 #' --model HKY --clock --no-force --dir --no-gBGC --switch',
                 #' --model HKY --clock --no-force --dir --gBGC --switch']
    with open('./Switch.sh', 'w+') as f:
        f.write('#!/bin/bash' + '\n')
        for pair in pairs:
            #save_txt = ' > Hessian_' + '_'.join(pair) + '_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + '_switch.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + '_switch.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                for ctl_line in ctl_lines:
                    ctl_line = ctl_line.replace('--clock', '--no-clock')
                    #g.write('python Run_unfinished.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ctl_line + '\n')
                    
                    g.write('python Run_unfinished.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ctl_line.replace('HKY', 'MG94') + '\n')
                    
 
