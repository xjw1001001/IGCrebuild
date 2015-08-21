if __name__ == '__main__':
    
    pairs = []
    all_pairs = '../All_Pairs.txt'
    jeff_pairs = './Jeff_pairs_list.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    if ['YLR028C', 'YMR120C'] in pairs:
        pairs.remove(['YLR028C', 'YMR120C'])

    sh_line = 'sbatch -o cd-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/'
    ctl_lines = [' --model HKY --clock --force --no-dir --no-gBGC',
                 ' --model HKY --clock --no-force --no-dir --no-gBGC',
                 ' --model HKY --clock --no-force --no-dir --gBGC',
                 ' --model HKY --clock --no-force --dir --no-gBGC',
                 ' --model HKY --clock --no-force --dir --gBGC',
                 ' --model HKY --no-clock --force --no-dir --no-gBGC',
                 ' --model HKY --no-clock --no-force --no-dir --no-gBGC',
                 ' --model HKY --no-clock --no-force --no-dir --gBGC',
                 ' --model HKY --no-clock --no-force --dir --no-gBGC',
                 ' --model HKY --no-clock --no-force --dir --gBGC']
    with open('./SitePosterior.sh', 'w+') as f:
        f.write('#!/bin/bash' + '\n')
        for pair in pairs:
            #save_txt = ' > Hessian_' + '_'.join(pair) + '_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + '_SitePosterior.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + '_SitePosterior.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                for ctl_line in ctl_lines:
                    g.write('python get_site_posterior.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ctl_line + '\n')
                    
                for ctl_line in ctl_lines:
                    g.write('python get_site_posterior.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ctl_line.replace('HKY', 'MG94') + '\n')
                    
 
