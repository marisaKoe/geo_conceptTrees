'''
Created on 01.08.2017

@author: marisakoe
'''
from numpy import *
import re
import glob


def read_pmi_matrices(method):
    '''
    read all the distance matrices from a folder
    '''
    list_matrices = glob.glob("input/nelex/"+method+"/phylip/*.phy")
    return list_matrices

def transform_distance_matrices(matrix):
    """
    !!!When language names are stored in an extra file, delete the comments in the corresponding line
    Read interleaved phylip distance matrix into a numpy matrix.
    Return taxon names, matrix
    """
    #only needed when the language names are stored in an extra file
    #names, longnames = get_language_sample() 
    #reads the matrix into a list
    elements = open(matrix,'rU').read().split()
    #removes the first row (the number of taxa in the matrix)
    N = int(elements.pop(0))
    #creates a list for the taxa
    taxa = []
    #fills the default matrix with -inf
    S = asmatrix(-Inf*ones([N,N]))
    #got through the list and save the distances in the matrix
    for row, i in enumerate(range(0, len(elements), N+1)):
        #append the first element of the list to the list of taxa
        taxa.append(elements[i])
        #fill the matrix
        S[row,:] = [float(x) for x in elements[i+1:i+N+1]]
    
    #ONLY needed when the language names are in an extra file
    #taxa = longnames.tolist()
    #returns the taxa and the matrix
    return taxa, S

def phylip_output(m,names=[],file='data.phy'):
    '''
    auxiliary function for writing the distance matrices to a file
    :param m: matrix
    :type m:an array
    :param names:languages
    :type names:list
    :param file:a file
    :type file:file
    '''
    if len(names) != len(m):
        names = xrange(len(m))
        mx = 10
    else:
        mx = max([len(x) for x in names])
    f = open(file,'w')
    f.write(str(len(names))+'\n')
    for nm, row in zip(names,m):
        f.write(str(nm).ljust(mx))
        for cell in row:
            f.write(' '+str(cell))
        f.write('\n')
    f.close()


if __name__ == '__main__':
    pass