if __name__ == '__main__':
    pairs = []
    with open('../All_Pairs.txt', 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    pairs.remove(['YLR028C', 'YMR120C'])
    Force_str = '{5:0.0}'

##    sh_line = 'sbatch -o cd-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/'
##    
##    with open('./DirNPack.sh', 'w+') as f:
##        f.write('#!/bin/bash' + '\n')
##        for pair in pairs:
##            save_txt = ' > Dir_' + '_'.join(pair) + '_PrintScreen.txt'
##            f.write(sh_line + '_'.join(pair) + 'dir.sh' + '\n')
##            with open('./NewRun/' + '_'.join(pair) + 'dir.sh', 'w+') as g:
##                g.write('#!/bin/bash' + '\n')
##                g.write('python DirGeneconv.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + save_txt + '\n')
##            
##    with open('./DirNPack_ForceTau.sh', 'w+') as f:
##        f.write('#!/bin/bash' + '\n')
##        for pair in pairs:
##            save_txt = ' > ForceDir_' + '_'.join(pair) + '_PrintScreen.txt'
##            f.write(sh_line + '_'.join(pair) + 'dir_ForceTau.sh' + '\n')
##            with open('./NewRun/' + '_'.join(pair) + 'dir_ForceTau.sh', 'w+') as g:
##                g.write('#!/bin/bash' + '\n')
##                g.write('python DirGeneconv.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --Force '+ Force_str + save_txt + '\n')


    sh_line = 'sbatch -o summary-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/'
    clock_suffix = ' --clock True --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/'
    nonclock_suffix = ' --clock False --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/'
    
    with open('./CdNPack_summary.sh', 'w+') as f:
        f.write('#!/bin/bash' + '\n')
        for pair in pairs:
            #save_txt = ' > Dir_' + '_'.join(pair) + '_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + 'summary.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + 'summary.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python GenerateIndvidualSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --force False' + clock_suffix + '\n')
                g.write('python GenerateIndvidualSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --force False' + nonclock_suffix + '\n')
                g.write('python GenerateIndvidualSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --force True' + clock_suffix + '\n')
                g.write('python GenerateIndvidualSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --force True' + nonclock_suffix + '\n')
            
