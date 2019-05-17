'''
Created on 06.03.2017

@author: marisakoe
'''

##create matrices combined with geo distances for the northeuralex+ielex data
#from IELexNLex import create_dm

##create matrices combined with geo distances for nelex
from nelex import create_dm

if __name__ == '__main__':
    #NLex+IELex data
    #create_dm()
    
    ##nelex
    ##inputFolder: input/nelex/phylip
    methods = ["cog"]
    
    for method in methods:
        create_dm(method)