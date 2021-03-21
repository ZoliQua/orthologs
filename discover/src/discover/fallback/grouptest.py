import numpy
import scipy.stats

from . import poisbinom


def firstBitSet(x):
    assert x > 0
    i = 0
    while not x & 1 == 1:
        x >>= 1
        i += 1
    return i


def isPowerOfTwo(x):
    return (x != 0) and not (x & (x - 1))


m1  = 0x5555555555555555
m2  = 0x3333333333333333
m4  = 0x0f0f0f0f0f0f0f0f
h01 = 0x0101010101010101

def popcount(x):
    x -= (x >> 1) & m1             # put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2) # put count of each 4 bits into those 4 bits
    x = (x + (x >> 4)) & m4        # put count of each 8 bits into those 8 bits
    return (x * h01)>>56           # returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...


def groupTest(events, bg):
    features = range(events.shape[0])

    lbg = numpy.log(bg)
    numSubsets = 2**events.shape[0]
    observed = numpy.empty((numSubsets, events.shape[1]), dtype=int)
    expected = numpy.empty((numSubsets, events.shape[1]))
    for subset in xrange(1, numSubsets):
        first = firstBitSet(subset)
        if isPowerOfTwo(subset):
            observed[subset] = events[first]
            expected[subset] = lbg[first]

        else:
            rest = subset ^ (1 << first)
            assert 0 < rest < subset
            observed[subset] = events[first] * observed[rest]
            expected[subset] = lbg[first] + expected[rest]

    p = numpy.array([
            poisbinom.cdf(numpy.exp(expected[i]), observed[i].sum())
            for i in xrange(1, numSubsets) if popcount(i) >= 2])

    expected[0] = 0

    cov = numpy.exp(numpy.array([
        [numpy.logaddexp.reduce(expected[i | j] + numpy.log1p(-numpy.exp(expected[i & j])))
         for i in xrange(1, numSubsets) if popcount(i) >= 2]
        for j in xrange(1, numSubsets) if popcount(j) >= 2]))

    var = numpy.diag(cov)

    corr = cov / numpy.sqrt(var)[:, numpy.newaxis] / numpy.sqrt(var)

    return scipy.stats.norm.cdf(scipy.stats.norm.ppf(p).sum() / numpy.sqrt(corr.sum()))
