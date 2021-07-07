'''
Created on 06.03.2017

@author: marisakoe

The files in the input folder are samples.
Please add your own files to the folder in order to compute the SDM.
'''


##create matrices combined with geo distances for nelex
from nelex import create_dm

if __name__ == '__main__':
    
    ##nelex
    ##inputFolder: input/nelex/phylip
    ##method of the concept tree reconstruction
    methods = ["cog", "pmi", "sigmoid"]
    
    for method in methods:
        create_dm(method)