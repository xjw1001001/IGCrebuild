from GenerateSummary import *

def get_individual_summary(pair, pair_path, model, summary_path, clock, force):
    if force:
        prefix = pair_path + 'Force_' + model + '_'
        prefix_summary = summary_path + 'Force_' + model + '_'
    else:
        prefix = pair_path + model + '_'
        prefix_summary = summary_path + model + '_'

    if clock:
        suffix = '_clock.p'
        suffix_summary = '_clock_summary.txt'
    else:
        suffix = '_nonclock.p'
        suffix_summary = '_nonclock_summary.txt'    

    p_file = prefix + '_'.join(pair) + suffix
    summary_file = prefix_summary + '_'.join(pair) + suffix_summary
    if os.path.isfile(p_file):
        res = get_summary(p_file, True)
        summary = np.matrix(res[0])
        label = res[1]
    else:
        summary = np.matrix(-1.234)
        label = ['unfinished']
        
    footer = ' '.join(label)  # row labels
    np.savetxt(open(summary_file, 'w+'), summary.T, delimiter = ' ', footer = footer)

def main(args):
    pair = [args.paralog1, args.paralog2]
    clock = args.clock
    force = args.force
    summary_path = args.sump
    model = args.model
    pair_path = args.pairp

    get_individual_summary(pair, pair_path, model, summary_path, clock, force)

if __name__ == '__main__':
##    pair = ['YAL056W', 'YOR371C']
##    clock = True
##    force = False
##    summary_path = './'
##    model = 'HKY'
##    pair_path = './NewPackageNewRun/'
##
##    get_individual_summary(pair, pair_path, model, summary_path, clock, force)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
    parser.add_argument('--force', required = True, help = 'Tau parameter control')
    parser.add_argument('--clock', required = True, help = 'clock parameter control')
    parser.add_argument('--model', required = True, help = 'model control')
    parser.add_argument('--sump', required = True, help = 'summary file path')
    parser.add_argument('--pairp', required = True, help = 'pair file path')
    main(parser.parse_args())
