seqfile = aln.phy
treefile = tree.tre
outfile = mlc_link
   
        noisy = 3              * How much rubbish on the screen
      verbose = 1              * More or less detailed report

      seqtype = 1              * Data type
        ndata = NDAT           * Number of data sets or loci
        icode = 0              * Genetic code 
    cleandata = 0              * Remove sites with ambiguity data?
		
model = 0
NSsites = 2
    CodonFreq = CODFREQ        * Codon frequencies
	  estFreq = ESTFREQ        * Use observed freqs or estimate freqs by ML
        clock = CLOCK          * Clock model
fix_omega = 0
omega = 1
