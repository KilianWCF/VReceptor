#!/usr/bin/env python
# coding: utf-8

# Virtual Receptor
# The vrecptor class below simulates a receptor based on a reference sequence and a function defining the pharmacophores
# Written by Kilian Conde-Frieboes, version below is from 2025-07-01


import numpy as np
import scipy
from scipy.spatial import distance
import pandas as pd
from sklearn.preprocessing import OneHotEncoder

import random
from packaging.version import Version
import sklearn
from sklearn.base import BaseEstimator, TransformerMixin


cosine = scipy.spatial.distance.cosine





aminoacids = {
    'A':['C3H5NO',{'C':3,'H':5,'N':1,'O':1,'S':0},[0.24,-2.32,0.6,-0.14,1.3],        [2.34,	9.69]],
    'C':['C3H5NOS',{'C':3,'H':5,'N':1,'O':1,'S':1},[0.84,-1.67,3.71,0.18,-2.65],     [1.96,	10.28,	-8.18]],
    'D':['C4H5NO3',{'C':4,'H':5,'N':1,'O':3,'S':0},[3.98,0.93,1.93,-2.46,0.75],      [1.88,	9.60,	-3.65]],
    'E':['C5H7NO3',{'C':5,'H':7,'N':1,'O':3,'S':0},[3.11,0.26,-0.11,-3.04,-0.25],    [2.19,	9.67,	-4.25]],
    'F':['C9H9NO',{'C':9,'H':9,'N':1,'O':1,'S':0},[-4.22,1.94,1.06,0.54,-0.62],      [1.83,	9.13]],
    'G':['C2H3NO',{'C':2,'H':3,'N':1,'O':1,'S':0},[2.05,-4.06,0.36,-0.82,-0.38],     [2.34,	9.60]],
    'H':['C6H7N3O',{'C':6,'H':7,'N':3,'O':1,'S':0},[2.47,1.95,0.26,3.9,0.09],        [1.82,	9.17,	6.00]],
    'I':['C6H11NO',{'C':6,'H':11,'N':1,'O':1,'S':0},[-3.89,-1.73,-1.71,-0.84,0.26],  [2.36,	9.60]],
    'K':['C6H12N2O',{'C':6,'H':12,'N':2,'O':1,'S':0},[2.29,0.89,-2.49,1.49,0.31],    [2.18,	8.95,	10.53]],
    'L':['C6H11NO',{'C':6,'H':11,'N':1,'O':1,'S':0},[-4.28,-1.3,-1.49,-0.72, 0.84],  [2.36,	9.60]],
    'M':['C5H9NOS',{'C':5,'H':9,'N':1,'O':1,'S':1},[-2.85,-0.22,0.47,1.94,-0.98],    [2.28,	9.21]],
    'N':['C4H6N2O2',{'C':4,'H':6,'N':2,'O':2,'S':0},[3.05,1.62,1.04,-1.15,1.61],     [2.02,	8.80]],
    'P':['C5H7NO',{'C':5,'H':7,'N':1,'O':1,'S':0},[-1.66,0.27,1.84,0.7,2.0],         [1.99,	10.60]],
    'Q':['C5H8N2O2',{'C':5,'H':8,'N':2,'O':2,'S':0},[1.75,0.5,-1.44,-1.34,0.66],     [2.17,	9.13]],
    'R':['C6H12N4O',{'C':6,'H':12,'N':4,'O':1,'S':0},[3.52,2.5,-3.5,1.99,-0.17],     [2.17,	9.04,	12.48]],
    'S':['C3H5NO2',{'C':3,'H':5,'N':1,'O':2,'S':0},[2.39,-1.07,1.15,-1.39,0.67],     [2.21,	9.15]],
    'T':['C4H7NO2',{'C':4,'H':7,'N':1,'O':2,'S':0},[0.75,-2.18,-1.12,-1.46,-0.4],    [2.09,	9.10]],
    'V':['C5H9NO',{'C':5,'H':9,'N':1,'O':1,'S':0},[-2.59,-2.64,-1.54,-0.85,-0.02],   [2.32,	9.62]],
    'W':['C11H10N2O',{'C':11,'H':10,'N':2,'O':1,'S':0},[-4.36,3.94,0.59,3.44,-1.59], [2.83,	9.39]],
    'Y':['C9H9NO2',{'C':9,'H':9,'N':1,'O':2,'S':0},[-2.54,2.44,0.43,0.04,-1.47],     [2.20,	9.11,	-10.07]],
    ' ':[''       ,{'C':0,'H':0,'N':0,'O':0,'S':0},[0.0, 0.0, 0.0, 0.0, 0.0],     []],
}

def func(x, *args):
    ''' Flexible function for a sum of gaussians

    Parameters
    ----------
    x : numpy array
        x values for the function
    *args : list
        First argument is the number of gaussians, second is the constant c, then for
        each gaussian the parameters mu, sig and var are given in triplets.
        Example: func(x, 2, 0, 45, 6.4, 89.3, 40.5, 8.9, 88) will return a sum of two gaussians
        with c=0, first gaussian with mu=45, sig=6.4 and var=89.3, second gaussian with mu=40.5, sig=8.9
    Returns
    -------
    numpy array
        The result of the function evaluated at x
    '''
    if len(args) < 2:
        print ('Parameter missing!')
        return None
    number_of_gaussian = int(args[0])
    c = args[1]
    if number_of_gaussian < 1:
        return c
    elif (not (len(args)-2)%3 == 0) or (len(args) == 0):
        print ('Gauss Parameter missing!')
        return None
    mu = []
    sig = []
    var = []
    for i in range (0, number_of_gaussian*3,3):
        mu.append(args[i+2])
        sig.append(args[i+3])
        var.append(args[i+4])
    z = c
    for i in range(number_of_gaussian):
        z = z + gauss(x, mu=mu[i], sig = sig[i]) * var[i] 
    return z


def gauss(x, mu=0, sig=1):
    ''' Gaussian function
    Parameters
    ----------
    x : numpy array
        x values for the function
    mu : float, optional
        mean of the gaussian, by default 0
    sig : float, optional
        standard deviation of the gaussian, by default 1
    Returns
    -------
    numpy array
        The result of the gaussian function evaluated at x
    '''
    return 1/(sig * np.sqrt(2 * np.pi)) * np.exp(-(x - mu)**2 / (2 * sig**2))


aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
def cosine_dist(aa, aminoacids = aminoacids):
    ''' Cosine distance between amino acids based on their properties
    Parameters
    ----------
    aa : list of characters
        amino acids to calculate the distance for
    aminoacids : dict, optional
        dictionary of amino acids with their properties, by default aminoacids
    Returns
    -------
    dict
        dictionary with cosine distance between amino acids, keys are amino acids, 
        values are dictionaries with other amino acids as keys and cosine distance as values
    '''
    return {j: {i: cosine(aminoacids[i][2],aminoacids[j][2]) for i in aa} for j in aa}



def chebyshev_dist(aa, aminoacids = aminoacids):
    ''' Chebyshev distance between amino acids based on their properties
    Parameters
    ----------
    aa : list of characters
        amino acids to calculate the distance for
    aminoacids : dict, optional
        dictionary of amino acids with their properties, by default aminoacids
    Returns
    -------
    dict
        dictionary with chebyshev distance between amino acids, keys are amino acids, 
        values are dictionaries with other amino acids as keys and chebyshev distance as values
    '''
    z=[0]*5
    for i in aa:
        for j in aa:
            diff=np.array(aminoacids[i][2])-np.array(aminoacids[j][2])
            z=[max(z[i],np.abs(diff[i])) for i in range (len(z))]
    sim = {j: {i: max(np.abs(np.array(aminoacids[i][2])-np.array(aminoacids[j][2])))/z[np.argmax(np.abs(np.array(aminoacids[i][2])-np.array(aminoacids[j][2])))]*2 
               for i in aa} for j in aa}
    return sim


def L2_dist(aa, aminoacids = aminoacids):
    ''' L2 distance between amino acids based on their properties
    Parameters
    ----------
    aa : list of characters
        amino acids to calculate the distance for
    aminoacids : dict, optional
        dictionary of amino acids with their properties, by default aminoacids
    Returns
    -------
    dict
        dictionary with L2 distance between amino acids, keys are amino acids, 
        values are dictionaries with other amino acids as keys and L2 distance as values
    '''
    zscale_v = {}
    for a in aa:
        v = aminoacids[a][2]    
        zscale_v[a] = v/np.linalg.norm(v)
        
    return {j: {i: np.linalg.norm(zscale_v[i]-zscale_v[j]) for i in aa} for j in aa}
    


class vreceptor:
    ''' Virtual Receptor class, simulates a receptor based on a reference sequence and a function defining the pharmacophores
    Parameters
    ----------
    ref : str
        reference sequence of the ligand
    *args : list
        First argument is the number of gaussians, second is the constant c, then for
        each gaussian the parameters mu, sig and var are given in triplets.
        Example: vreceptor(ref, 2, 0, 45, 6.4, 89.3, 40.5, 8.9, 88) will create a receptor with two gaussians
        with c=0, first gaussian with mu=45, sig=6.4 and var=89.3, second gaussian with mu=40.5, sig=8.9
    func : function, optional
        function to use for the pharmacophores, by default func
    aa : list of characters, optional
        amino acids to use for the ligand, by default aa
    aa_norm : function, optional
        function to calculate distance between the amino acids, by default cosine_sim
    Attributes
    ----------
    reference : str
        reference sequence of the ligand
    args : list
        arguments for the function defining the pharmacophores
    weights : numpy array
        weights for the pharmacophores, calculated based on the reference sequence and the function func
    aa : list of characters
        amino acids used for the ligand
    '''

    def __init__(self, ref, *args, func=func, aa=aa, aa_norm=cosine_dist):
        # initialization of the receptor, checks the input and calculates the weights for the pharmacophores
        if not type(ref) == str:
            print ('Please enter a peptide as string!')
            return None
        aa_norm = aa_norm(aa)
        self.aa = aa
        self.reference = ref
        self.args = args
        x = np.array([i for i in range(0, len(ref))])
        y = func(x, *args)
        self.weights = np.array([[aa_norm[j][i] for i in aa] for j in self.reference])
        for i in range(len(self.weights)):
            self.weights[i] = self.weights[i] * y[i]

    def virtual_assay(self, pep_list, csv_out=None, noise=0, seed=42):
        ''' Virtual assay for the receptor, calculates the logIC50 values for a list of peptides
        Parameters
        ----------
        pep_list : list of str
            list of peptides to calculate the logIC50 values for
        csv_out : str, optional
            path to the output csv file, if None, the output is not saved to a file, by default None
        noise : float, optional
            noise to add to the logIC50 values, by default 0
        seed : int, optional
            seed for the random number generator, by default 42
        Returns
        -------
        pandas DataFrame
            DataFrame with the peptides and their logIC50 values, if csv_out is not None, the DataFrame is saved to a csv file
        '''
        if not type(pep_list) == list:
            print ('Please enter a list of peptides as string!')
            return None
        np.random.seed(seed)
        aa_ = [self.aa for i in range(len(self.reference))]
        # latest version of OneHotEncoder is sklearn uses sparse_output instead of sparse
        if Version(sklearn.__version__) > Version('1.1.2'):
            enco = OneHotEncoder(sparse_output = False, categories = aa_)
        else:
            enco = OneHotEncoder(sparse = False, categories = aa_)
        X1 = [list(i) for i in pep_list]
        X1 = enco.fit_transform(X1)
        X1 = X1.reshape(X1.shape[0], int(X1.shape[1] / len(aa)), len(aa))
        results = np.multiply(X1, self.weights)
        df = pd.DataFrame({'Sequence': pep_list, 'logIC50' : np.array([np.sum(i)+ np.random.normal(0, noise) for i in results])})
        if csv_out != None:
            df.to_csv(csv_out, index=False, sep=',')
        return df





Receptor_1=vreceptor('SRTHRHSMEIRTPDINPAWYASRGIRPVGRF', *[2, 0, 45, 6.4, 89.3, 40.5, 8.9, 88], aa_norm=chebyshev_dist)
Receptor_2=vreceptor('SRTHRHSMEIRTPDINPAWYASRGIRPVGRF', *[2, 0, 130, 34.8, -2547, 29.5, 0.12, 6249], aa_norm=chebyshev_dist)


def peptide_generator(pep, number_of_mut=1, start = 0, stop=None, aa=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'], exclude_C=True):    
    '''
    Peptide generator, enumerates all mutants with exact number_of_mut mutations in order given the aa list provided
    Parameters
    ----------
    pep : string
        peptide to mutate        
    number_of_mut : int
        number of mutations
    start : int
        first position to mutate (zero based)
    stop : int
        last position to mutate, if None stop=len(pep) (dfault)
    aa : list of characters
        amino acids to mutate with
    exclude_C : boolean
        if True, Cysteine positions are not changed (default)
    '''
    if stop == None:
        stop = len(pep)
    if number_of_mut == 0:
        yield pep
    else:
        for aminoacid in aa:
            for n in range(start, stop):
                if aminoacid == pep[n]:
                    continue
                if (pep[n] == 'C') and (exclude_C == True):
                    continue
                dummy = pep[:n]+aminoacid+pep[n+1:]
                yield from peptide_generator(dummy, number_of_mut=number_of_mut-1, aa=aa, start=n+1, stop=stop)
