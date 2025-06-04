#!/usr/bin/env python

"""
Create HMMs.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import argparse

class CreateHMMs(object):
  def __init__(self):
    self.prokSeqs = '/srv/whitlam/bio/db/mothur/gg_12_10/aligned_seqs/94_otu.align.fna'
    self.taxonomyFile = '/srv/whitlam/bio/db/mothur/gg_12_10/taxonomy/94_otu_taxonomy.txt'
    self.eukSeqs = '/srv/whitlam/bio/db/mothur/silva/silva.eukarya.fasta'

    self.nuc = {'A':'T','T':'A','G':'C','C':'G','K':'M','M':'K','R':'Y','Y':'R','S':'W','W':'W','B':'V','V':'B','H':'G','D':'C','X':'N','N':'N','*':'N','-':'-', '.':'-'}
    pass

  def readTaxonomy(self, taxonomyFile):
    taxonomy = {}
    for line in open(taxonomyFile):
      lineSplit = line.split('\t')
      ggId = lineSplit[0]
      taxa = [x.strip() for x in lineSplit[1].split(';')]

      taxonomy[ggId] = taxa

    return taxonomy

  def createSeqFile(self, taxonomy, outputFile, taxa):
    fout = open(outputFile, 'w')

    bKeep = False
    for line in open(self.prokSeqs):
      if line[0] == '>':
        ggId = line[1:].rstrip()
        if taxa == taxonomy[ggId][0]:
          fout.write(line)
          bKeep = True
        else:
          bKeep = False
      elif bKeep:
        fout.write(line)
    fout.close()

    self.revComplementFastaFile(outputFile)

  def run(self, output_dir, threads):
    taxonomy = self.readTaxonomy(self.taxonomyFile)

    self.createSeqFile(taxonomy, 'bacteria.fasta', 'k__Bacteria')
    os.system('hmmbuild --cpu ' + str(threads) + ' --informat afa ' + output_dir + '/SSU_bacteria.hmm bacteria.fasta')

    self.createSeqFile(taxonomy, 'archaea.fasta', 'k__Archaea')
    os.system('hmmbuild --cpu ' + str(threads) + ' --informat afa ' + output_dir + '/SSU_archaea.hmm archaea.fasta')

    self.revComplementFastaFile(self.eukSeqs)
    os.system('hmmbuild --cpu ' + str(threads) + ' --informat afa ' + output_dir + '/SSU_euk.hmm ' + self.eukSeqs)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Create HMMs.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('output_dir', help='output directory')
  parser.add_argument('-t', '--threads', help='number of threads', type=int, default = 1)

  args = parser.parse_args()

  createHMMs = CreateHMMs()
  createHMMs.run(args.output_dir, args.threads)
