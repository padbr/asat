#! /usr/bin/python

__date__ = "02/2021"

def padjust(pvalues):
    """
    padjust(pvals)
    padjust([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1])
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    n = float(pvalues.shape[0])
    new_pvalues = empty(int(n))
    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in range(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues

def padjNaN(pvals):
    """
    Does Benjamini-Hochberg FDR on samples with NaN values
    """
    from numpy import isnan
    enum_pvals = list(enumerate(pvals))
    nanIndexes = [elem[0] for elem in enum_pvals if isnan(elem[1])]
    pvals_noNaN = [elem[1] for elem in enum_pvals if not isnan(elem[1])]
    padj_noNaN = list(padjust(pvals_noNaN))
    padj = list(enumerate(padj_noNaN))
    padj = [list(elem) for elem in padj]
    nanList = [[i, float('nan')] for i in nanIndexes]
    for i in nanIndexes:
        for elem in padj:
            if elem[0] >= i:
                elem[0] += 1
    padj = padj + nanList
    padj.sort(key=lambda x:x[0])
    padj = [elem[1] for elem in padj]
    return padj
