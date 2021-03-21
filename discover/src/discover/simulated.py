import networkx
import numpy
import scipy.stats


def generateNullData(geneMarginals, sampleMarginals):
    graph = networkx.bipartite_configuration_model(
        geneMarginals, sampleMarginals, create_using=networkx.Graph())
    events = numpy.asarray(networkx.algorithms.bipartite.biadjacency_matrix(
        graph, [n for n in graph if graph.node[n]["bipartite"] == 0]))
    return events


def generateCooc(n, coverage, impurity, impurityBalance=0.5):
    genes = numpy.zeros((2, n))
    subset = scipy.stats.bernoulli.rvs(coverage, size=n).astype(bool)
    genes[:, subset] = 1
    impure = scipy.stats.bernoulli.rvs(impurity, size=subset.sum()).astype(bool)
    impureGenes = scipy.stats.bernoulli.rvs(impurityBalance, size=impure.sum())
    genes[impureGenes, subset.nonzero()[0][impure]] = 0
    return genes


def generateMutex(n, coverage, impurity, coverageBalance=0.5):
    genes = numpy.zeros((2, n))
    subset = scipy.stats.bernoulli.rvs(coverage, size=n).astype(bool)
    affectedGenes = scipy.stats.bernoulli.rvs(coverageBalance, size=subset.sum())
    genes[affectedGenes, subset.nonzero()[0]] = 1
    impure = scipy.stats.bernoulli.rvs(impurity, size=n).astype(bool)
    genes[:, subset & impure] = 1
    return genes


def generateMutexGroup(numGenes, sampleMarginals, coverage=0.5, impurity=0.05):
    numSamples = len(sampleMarginals)
    numSamplesCov = int(coverage * numSamples)
    events = numpy.zeros((numGenes, numSamples))

    p = 1.0 * sampleMarginals / sampleMarginals.sum()
    alteredSamples = numpy.random.choice(numSamples, numSamplesCov, False, p)

    events[numpy.random.choice(numGenes, numSamplesCov), alteredSamples] = 1

    for i in xrange(events.shape[0]):
        altered = events[i].nonzero()[0]
        numSamplesImp = int(impurity * (numSamplesCov - len(altered)))
        availableSamples = numpy.setdiff1d(alteredSamples, altered)
        alteredSamplesImp = numpy.random.choice(availableSamples, numSamplesImp, False,
                                                p[availableSamples] / p[availableSamples].sum())
        events[i, alteredSamplesImp] = 1

    return events
