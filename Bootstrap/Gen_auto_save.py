from Pdiff import read_txt
from IGCexpansion.CodonGeneconv import ReCodonGeneconv


if __name__ == '__main__':
    pairs = []
    all_pairs = './Filtered_pairs.txt'
    summary_path = '/Users/xji3/MixedFromCluster/NewPackageNewRun/'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    for paralog in pairs:
        model = 'MG94'
        alignment_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
        newicktree = './YeastTree.newick'

        clock = False
        Force = None

        test = ReCodonGeneconv(newicktree, alignment_file, paralog, Model = model, Force = Force, clock = clock)
        x = read_txt(summary_path, paralog, model, Force, clock, False, False)
        test.update_by_x(x)
        test.save_x()
