      seqfile = for.codeml.ECP_EDN.fasta
     treefile = for.paml.newick.tree
      outfile = mlc

        noisy = 9
      verbose = 1
      runmode = 0


      seqtype = 1
    * CodonFreq = 0
      * estFreq = 1

    CodonFreq = 4 * use the MG94 codon model
      estFreq = 1

        ndata = 1
        clock = 0
       aaDist = 0



        model = 0






      NSsites = 0 




        icode = 0
        Mgene = 0



    fix_kappa = 0 * allow kappa to be estimated freely
        kappa = 5.86954
    fix_omega = 0 * allow omega to be estimated freely
        omega = 0.087136

    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 5

        getSE = 1
 RateAncestor = 0

   Small_Diff = 5e-7
    cleandata = 0
  fix_blength = 0 * allow branch lengths to be estimated freely
       method = 0
