import numpy

from itertools import combinations, combinations_with_replacement

from _discover import maxent, fdr


def count_unique(keys):
    # From http://stackoverflow.com/questions/10741346/numpy-frequency-counts-for-unique-values-in-an-array
    uniq_keys = numpy.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, numpy.bincount(bins)


def estimateBackground(events):
    rowSums = events.sum(1)
    colSums = events.sum(0)

    rowValues, rowWeights = count_unique(rowSums)
    colValues, colWeights = count_unique(colSums)

    mu = maxent.fit(rowValues, rowWeights, colValues, colWeights)
    nRows = len(rowValues)
    eA = numpy.dot(numpy.exp((mu[:nRows] / rowWeights)[:, numpy.newaxis]),
                   numpy.exp((mu[nRows:] / colWeights)[numpy.newaxis]))
    P = 1.0 / (eA + 1)

    if P.max() > 1:
        import warnings
        warnings.warn("Some background estimates are greater than 1")
        #P[P > 1] = 1

    return P[rowValues.searchsorted(rowSums)[:, numpy.newaxis],
             colValues.searchsorted(colSums)[numpy.newaxis]]


def estimateBackgroundStratified(events, strata):
    result = numpy.empty(events.shape)
    for x in numpy.unique(strata):
        result[:, strata == x] = estimateBackground(events[:, strata == x])
    return result


def analyse(events, subset=None, strata=None):
    if strata is None:
        bg = estimateBackground(events.data)
    else:
        bg = estimateBackgroundStratified(events.data, strata)

    return fdr.mutex(events[subset], bg[subset])


def labelResultMatrix(result, data1, data2):
    return biotk.la.LabelledArray(
        result,
        [data1.featureNames, data2.featureNames],
        ["genes1", "genes2"])


def analyseBlockStructure(data1, bg1, data2, bg2, lowerTail, symmetric=True):
    assert len(data1) == len(data2)

    if symmetric:
        pairs = combinations
    else:
        pairs = combinations_with_replacement

    numTests = sum(data1[i].shape[0] * data2[j].shape[0] for i, j in pairs(xrange(len(data1)), 2))
    pValues = numpy.empty(numTests)

    offset = 0
    for i, j in pairs(xrange(len(data1)), 2):
        if i != j:
            size = data1[i].shape[0] * data2[j].shape[0]
            if size > 0:
                view = pValues[offset:offset+size].reshape((data1[i].shape[0], data2[j].shape[0]))
                view[:] = fdr.computep(data1[i], bg1[i], data2[j], bg2[j], lowerTail)
                offset += size

    uniqueP, inv = numpy.unique(pValues, return_inverse=True)
    uniqueQ = numpy.zeros_like(uniqueP)
    expP0 = numpy.zeros(1)

    for i, j in pairs(xrange(len(data1)), 2):
        if i != j:
            fdr.updatemultiq(bg1[i], bg2[j], lowerTail, uniqueP, uniqueQ, expP0)

    counts = numpy.bincount(inv).cumsum()
    uniqueQ = numpy.minimum.accumulate(numpy.minimum(uniqueQ.cumsum() / counts, 1)[::-1])[::-1]

    qValues = uniqueQ[inv]

    qValueList = []
    pValueList = []
    offset = 0
    for i, j in pairs(xrange(len(data1)), 2):
        if i != j:
            size = data1[i].shape[0] * data2[j].shape[0]
            if size > 0:
                qValueList.append(labelResultMatrix(
                    qValues[offset:offset+size].reshape((data1[i].shape[0], data2[j].shape[0])),
                    data1[i], data2[j]))
                pValueList.append(labelResultMatrix(
                    pValues[offset:offset+size].reshape((data1[i].shape[0], data2[j].shape[0])),
                    data1[i], data2[j]))
                offset += size

    pi0 = pValues.sum() / expP0[0]

    return pValueList, qValueList, pi0
