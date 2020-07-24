#file for getting peptides from other enzymes


##import libraries
import matplotlib.pyplot as plt
import numpy as np
from pyteomics import fasta, parser

print('hello world')

databasePath = 'C:/Users/chris/OneDrive/Documents/bccrc/projectsRepository/sorensenLab/relatedToYbx1/design20200718_multipleProteaseSp3/uniprotHumanJul2020.fasta'

print('Cleaving the proteins with trypsin...')
uniquePeptides = set()
for description, sequence in fasta.read(databasePath):
    newPeptides = parser.cleave(sequence, 'K', min_length = 6)
    uniquePeptides.update(newPeptides)
print('Done, {0} sequences obtained!'.format(len(uniquePeptides)))
