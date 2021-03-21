import numpy


def coocParams(events):
    assert events.shape[0] == 2
    
    colSums = events.sum(0)
    rowSums = events.sum(1)
    covered = numpy.sum(colSums > 0)
    impurity = 1.0 * numpy.sum(colSums == 1) / covered
    balance = 1.0 * (covered - rowSums[1]) / (covered - rowSums).sum()
    
    return (1.0 * covered / events.shape[1],
            impurity, balance)


def mutexParams(events):
    assert events.shape[0] == 2
    
    colSums = events.sum(0)
    rowSums = events.sum(1)
    covered = numpy.sum(colSums > 0)
    impure = numpy.sum(colSums == 2)
    balance = 1.0 * (rowSums[1] - impure) / (rowSums - impure).sum()
    
    return (1.0 * covered / events.shape[1],
            1.0 * impure / covered,
            balance)
