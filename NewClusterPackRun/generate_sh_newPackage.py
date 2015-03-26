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


##    sh_line = 'sbatch -o summary-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/'
##    clock_suffix = ' --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/'
##    nonclock_suffix = ' --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/'
##    
##    with open('./CdNPack_summary.sh', 'w+') as f:
##        f.write('#!/bin/bash' + '\n')
##        for pair in pairs:
##            #save_txt = ' > Dir_' + '_'.join(pair) + '_PrintScreen.txt'
##            f.write(sh_line + '_'.join(pair) + 'summary.sh' + '\n')
##            with open('./NewRun/' + '_'.join(pair) + 'summary.sh', 'w+') as g:
##                g.write('#!/bin/bash' + '\n')
##                g.write('python GenerateIndvidualSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --no-force' + clock_suffix + '\n')
##                g.write('python GenerateIndvidualSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --no-force' + nonclock_suffix + '\n')
##                g.write('python GenerateIndvidualSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --force' + clock_suffix + '\n')
##                g.write('python GenerateIndvidualSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --force' + nonclock_suffix + '\n')



##    sh_line = 'sbatch -o summary-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/'
##    clock_suffix = ' --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/'
##    nonclock_suffix = ' --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/'
##    
##    with open('./DirNPack_summary.sh', 'w+') as f:
##        f.write('#!/bin/bash' + '\n')
##        for pair in pairs:
##            #save_txt = ' > Dir_' + '_'.join(pair) + '_PrintScreen.txt'
##            f.write(sh_line + '_'.join(pair) + 'dir_summary.sh' + '\n')
##            with open('./NewRun/' + '_'.join(pair) + 'dir_summary.sh', 'w+') as g:
##                g.write('#!/bin/bash' + '\n')
##                g.write('python GenerateIndividualDirSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --no-force' + clock_suffix + '\n')
##                g.write('python GenerateIndividualDirSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --no-force' + nonclock_suffix + '\n')
##                g.write('python GenerateIndividualDirSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --force' + clock_suffix + '\n')
##                g.write('python GenerateIndividualDirSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --force' + nonclock_suffix + '\n')



##    sh_line = 'sbatch -o summary-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/'
##    clock_suffix = ' --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/'
##    nonclock_suffix = ' --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/'
##    
##    with open('./gBGCNPack_summary.sh', 'w+') as f:
##        f.write('#!/bin/bash' + '\n')
##        for pair in pairs:
##            #save_txt = ' > Dir_' + '_'.join(pair) + '_PrintScreen.txt'
##            f.write(sh_line + '_'.join(pair) + 'gBGC_summary.sh' + '\n')
##            with open('./NewRun/' + '_'.join(pair) + 'gBGC_summary.sh', 'w+') as g:
##                g.write('#!/bin/bash' + '\n')
##                g.write('python GenerateIndividualgBGCSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --dir --no-force' + clock_suffix + '\n')
##                g.write('python GenerateIndividualgBGCSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --dir --no-force' + nonclock_suffix + '\n')
##                g.write('python GenerateIndividualgBGCSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --no-dir --no-force' + clock_suffix + '\n')
##                g.write('python GenerateIndividualgBGCSummary.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --no-dir --no-force' + nonclock_suffix + '\n')
              

##    sh_line = 'sbatch -o gBGC_Dir-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/'
##    
##    with open('./gBGCDirNPack.sh', 'w+') as f:
##        f.write('#!/bin/bash' + '\n')
##        for pair in pairs:
##            save_txt = ' > gBGC_Dir_' + '_'.join(pair) + '_PrintScreen.txt'
##            f.write(sh_line + '_'.join(pair) + 'dir_gBGC.sh' + '\n')
##            with open('./NewRun/' + '_'.join(pair) + 'dir_gBGC.sh', 'w+') as g:
##                g.write('#!/bin/bash' + '\n')
##                g.write('python gBGCDirGeneconv.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + save_txt + '\n')

##    sh_line = 'sbatch -o gBGC-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/'
##    
##    with open('./gBGCCdNPack.sh', 'w+') as f:
##        f.write('#!/bin/bash' + '\n')
##        for pair in pairs:
##            save_txt = ' > gBGC_' + '_'.join(pair) + '_PrintScreen.txt'
##            f.write(sh_line + '_'.join(pair) + '_gBGC.sh' + '\n')
##            with open('./NewRun/' + '_'.join(pair) + '_gBGC.sh', 'w+') as g:
##                g.write('#!/bin/bash' + '\n')
##                g.write('python gBGCCodonGeneconv.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + save_txt + '\n')
             

# Generate .sh files for unfinished pairs


    sh_line = 'sbatch -o summary-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/'
    clock_suffix = ' --clock --model MG94'
    nonclock_suffix = ' --no-clock --model MG94'

    pairs_nonclock = []
    pairs_clock = []
    pairs_nonclock_force = []
    pairs_clock_force = []
    pairs_nonclock_dir = []
    pairs_clock_dir = []
    pairs_nonclock_dir_gBGC = []
    pairs_clock_dir_gBGC = []
    pairs_nonclock_gBGC = []
    pairs_clock_gBGC = []
    summary_path = '/Users/xji3/FromCluster03212015/'
    with open(summary_path + 'NewPackageNewRun/MG94_nonclock_unfinished.txt', 'r') as f:
        for line in f.readlines():
            pairs_nonclock.append(line.replace('\n','').split('_'))
    
    with open(summary_path + 'NewPackageNewRun/MG94_clock_unfinished.txt', 'r') as f:
        for line in f.readlines():
            pairs_clock.append(line.replace('\n','').split('_'))
    
    with open(summary_path + 'NewPackageNewRun/Force_MG94_nonclock_unfinished.txt', 'r') as f:
        for line in f.readlines():
            pairs_nonclock_force.append(line.replace('\n','').split('_'))
    
    with open(summary_path + 'NewPackageNewRun/Force_MG94_clock_unfinished.txt', 'r') as f:
        for line in f.readlines():
            pairs_clock_force.append(line.replace('\n','').split('_'))

    with open(summary_path + 'NewPackageNewRun/Dir_MG94_nonclock_unfinished.txt', 'r') as f:
        for line in f.readlines():
            pairs_nonclock_dir.append(line.replace('\n','').split('_'))
    
    with open(summary_path + 'NewPackageNewRun/Dir_MG94_clock_unfinished.txt', 'r') as f:
        for line in f.readlines():
            pairs_clock_dir.append(line.replace('\n','').split('_'))

    with open(summary_path + 'NewPackageNewRun/Dir_gBGC_MG94_nonclock_unfinished.txt', 'r') as f:
        for line in f.readlines():
            pairs_nonclock_dir_gBGC.append(line.replace('\n','').split('_'))
    
    with open(summary_path + 'NewPackageNewRun/Dir_gBGC_MG94_clock_unfinished.txt', 'r') as f:
        for line in f.readlines():
            pairs_clock_dir_gBGC.append(line.replace('\n','').split('_'))

    with open(summary_path + 'NewPackageNewRun/gBGC_MG94_nonclock_unfinished.txt', 'r') as f:
        for line in f.readlines():
            pairs_nonclock_gBGC.append(line.replace('\n','').split('_'))
    
    with open(summary_path + 'NewPackageNewRun/gBGC_MG94_clock_unfinished.txt', 'r') as f:
        for line in f.readlines():
            pairs_clock_gBGC.append(line.replace('\n','').split('_'))

            
    with open('./CdNPack_unfinished.sh', 'w+') as f:
        f.write('#!/bin/bash' + '\n')
        for pair in pairs_nonclock:
            save_txt = ' > ' + '_'.join(pair) + '_nonclock_unfinished_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + '_nonclock_unfinished.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + '_nonclock_unfinished.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python Run_unfinished.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --no-force' + ' --no-dir' + nonclock_suffix + save_txt + '\n')
        for pair in pairs_clock:
            save_txt = ' > ' + '_'.join(pair) + '_clock_unfinished_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + '_clock_unfinished.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + '_clock_unfinished.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python Run_unfinished.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --no-force' + ' --no-dir' + clock_suffix + save_txt + '\n')
        for pair in pairs_nonclock_force:
            save_txt = ' > ' + '_'.join(pair) + '_nonclock_force_unfinished_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + '_nonclock_force_unfinished.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + '_nonclock_force_unfinished.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python Run_unfinished.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --force' + ' --no-dir' + nonclock_suffix + save_txt + '\n')
        for pair in pairs_clock_force:
            save_txt = ' > ' + '_'.join(pair) + '_clock_force_unfinished_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + '_clock_force_unfinished.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + '_clock_force_unfinished.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python Run_unfinished.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --force' + ' --no-dir' + clock_suffix + save_txt + '\n')
        for pair in pairs_nonclock_dir:
            save_txt = ' > ' + '_'.join(pair) + '_nonclock_dir_unfinished_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + '_nonclock_dir_unfinished.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + '_nonclock_dir_unfinished.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python Run_unfinished.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --no-force' + ' --dir' + nonclock_suffix + save_txt + '\n')
        for pair in pairs_clock_dir:
            save_txt = ' > ' + '_'.join(pair) + '_clock_dir_unfinished_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + '_clock_dir_unfinished.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + '_clock_dir_unfinished.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python Run_unfinished.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --no-force' + ' --dir' + clock_suffix + save_txt + '\n')

        for pair in pairs_nonclock_dir_gBGC:
            save_txt = ' > ' + '_'.join(pair) + '_nonclock_dir_gBGC_unfinished_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + '_nonclock_dir_gBGC_unfinished.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + '_nonclock_dir_gBGC_unfinished.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python Run_unfinished.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --no-force --gBGC' + ' --dir' + nonclock_suffix + save_txt + '\n')
        for pair in pairs_clock_dir_gBGC:
            save_txt = ' > ' + '_'.join(pair) + '_clock_dir_gBGC_unfinished_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + '_clock_dir_gBGC_unfinished.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + '_clock_dir_gBGC_unfinished.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python Run_unfinished.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --no-force --gBGC' + ' --dir' + clock_suffix + save_txt + '\n')

        for pair in pairs_nonclock_gBGC:
            save_txt = ' > ' + '_'.join(pair) + '_nonclock_gBGC_unfinished_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + '_nonclock_gBGC_unfinished.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + '_nonclock_gBGC_unfinished.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python Run_unfinished.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --no-force --gBGC' + ' --no-dir' + nonclock_suffix + save_txt + '\n')
        for pair in pairs_clock_gBGC:
            save_txt = ' > ' + '_'.join(pair) + '_clock_gBGC_unfinished_PrintScreen.txt'
            f.write(sh_line + '_'.join(pair) + '_clock_gBGC_unfinished.sh' + '\n')
            with open('./NewRun/' + '_'.join(pair) + '_clock_gBGC_unfinished.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python Run_unfinished.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' --no-force --gBGC' + ' --no-dir' + clock_suffix + save_txt + '\n')
