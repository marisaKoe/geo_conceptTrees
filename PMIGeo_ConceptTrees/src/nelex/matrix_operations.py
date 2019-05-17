'''
Created on 01.08.2017

@author: marisakoe


+ read both matrices in -> getting them as numpy matrices
+ make sure both matrices contain the same set of languages -> edit the geo dist matrix
+ if both matrices are in numpy sum them up

+ save & take the resulting matrix to reconstruct trees with ape
'''

from collections import defaultdict

from itertools import combinations
from rpy2.robjects.packages import importr

from helper_functions import *

#import the basis of R
base = importr("base")
utils = importr("utils")
stats = importr("stats")
#imports ape, required for phangorn
ape = importr("ape")
#imports phangorn
phangorn = importr("phangorn")


def create_dm(method):
    #get all the pmi matrices
    all_pmi_matrices = read_pmi_matrices(method)
    #print all_pmi_matrices
    #get the geo distance matrix
    geo_dm = "input/nelex/distances_langPairs.phy"
    #print geo_dm
    geo_taxa, geo_matrix = transform_distance_matrices(geo_dm)
    #for each matrix in the pmi matrices
    #transform the matrix into numpy
    #write a helper function to edit the matrix
    #sum the matrices up and write them into a file
    all_geo_matrices = defaultdict()
    for mtx in all_pmi_matrices:
        #get the concept of the matrix
        concept = mtx.split("/")
        concept = concept[-1]
        #transform the pmi matrix into a numpy array
        pmi_taxa, pmi_matrix = transform_distance_matrices(mtx)
        #if the language sample is different
        if not len(pmi_taxa) == len(geo_taxa):
            #helper funtion to edit the matrix
            new_geo_taxa, new_geo_matrix = edit_matrix(geo_taxa, geo_matrix, pmi_taxa)
            #print new_geo_matrix
            #writes the new geo_matrix to a file (can also be outcommended)
            matrixArray = squeeze(asarray(new_geo_matrix))
            folder = "output/nelex/"+method+"/geomatrices/geo+"+concept
            phylip_output(matrixArray, pmi_taxa, folder)
            all_geo_matrices[folder] = pmi_taxa

        else:
            #writes the geo matrix for the concept into a file (can also be outcommended)
            matrixArray = squeeze(asarray(geo_matrix))
            folder = "output/nelex/"+method+"/geomatrices/geo+"+concept
            phylip_output(matrixArray, pmi_taxa, folder)
            all_geo_matrices[folder] = pmi_taxa
    
    #use the method for combining the matrices, write the supermatrix, create the trees (NJ & FastME), write the trees in files
    combine_matrices(all_pmi_matrices, all_geo_matrices, method)
    
def edit_matrix(geo_taxa,geo_matrix, pmi_taxa):
    '''
    edit the geo matrix to the same length and with the same samples as in the pmi matrix
    needed for combining both matrices
    :param geo_taxa: the taxa sample of the geo dist
    :param geo_matrix:the geo dist matrix
    :param pmi_taxa:the pmi taxa
    '''
    #get the difference between the two taxa lists
    diff = list(set(geo_taxa)-set(pmi_taxa))
    #print diff
    #get the intersection of the two lists to form a new taxa list
    newTaxa = list(set(geo_taxa) & set(pmi_taxa))
    #for each taxa in the difference set, get the index and save them to a list
    diffIndex = []
    for t in diff:
        i = geo_taxa.index(t)
        diffIndex.append(i)
    #print diffIndex
    #delete all rows and colums at the indices
    matrixNew = delete(geo_matrix,diffIndex,0)
    matrixNew = delete(matrixNew,diffIndex,1)
    
    #print len(matrixNew)         
    #return the new taxa sample and the new matrix
    return newTaxa, matrixNew     
    
def combine_matrices(all_pmi_matrices, all_geo_matrices, method):
    '''
    Combines the matrices of two same concepts with the SDM method.
    The method involves brining the matrices closer (in the least-square-sense) relative to the others.
    "SDM deforms the source matrices, without modifying their topological message, to bring them as close as possible to each other;
    these deformed matrices are then averaged to obtain the distance supermatrix." (Criscoulo et al. SDM, 2006)
    The new matrix is saved.
    Trees are build using the NJ and fastme algorithm in ape and are saved accordingly.
    :param all_pmi_matrices:a list of all pmi matrices
    :param all_geo_matrices:a dictionary of all geo matrices and corresponding taxa list key=file name of matrix value=list of taxa in the sample
    '''
    #count for checking
    count = 0
    #for each matrix file and taxa list in the geo matrix dictionary
    for mtx,taxa in all_geo_matrices.items():
        #get the concept for comparison (includes .phy ending)
        concept = mtx.split("/")[-1].split("+")[-1]
        #print concept
        #for each matrix in the pmi list
        for mtx2 in all_pmi_matrices:
            #get the concept (includes .phy ending)
            concept2 = mtx2.split("/")[-1]
            #check if the concepts are the same
            if concept == concept2:
                #count for checking
                count += 1
                #get the concept name without the .phy ending 8for writing the trees)
                conceptName = concept.split(".")[0]
                ###use all the functions in R
                #read it in geo matrix
                t = utils.read_table(mtx,header=False, skip=1, row_names=1)
                #make a matrix
                mx = base.as_matrix(t)
                #make a distance matrix
                dm = stats.as_dist(mx)
                 
                #read it in the pmi matrix
                t2 = utils.read_table(mtx2,header=False, skip=1, row_names=1)
                #make a matrix
                mx2 = base.as_matrix(t2)
                #make a distance matrix
                dm2 = stats.as_dist(mx2)
                 
                #combine the matrices with the super distance matrix method from ape
                dd = ape.SDM(dm, dm2,len(taxa),len(taxa))
                #create a numpy array out of it
                matrixArray = squeeze(asarray(dd[0]))
                #write the matrix into a file
                phylip_output(matrixArray, taxa, "output/nelex/"+method+"/combinedMatrices/"+concept)
                #create a NJ tree out of the supermatrix
                treeNJ = ape.nj(dd[0])
                #write the tree in a file
                ape.write_tree(treeNJ, file="output/nelex/"+method+"/conceptTrees/nj/"+conceptName+"+njTree.nwk")
                
                #create a fastme tree out of the supermatrix
                treeFastme = ape.fastme_bal(dd[0],nni=True, spr=True, tbr=False)
                #write the tree in a file
                ape.write_tree(treeFastme, file="output/nelex/"+method+"/conceptTrees/fastme/"+conceptName+"+fastmeTree.nwk")
    print count          

if __name__ == '__main__':
    pass