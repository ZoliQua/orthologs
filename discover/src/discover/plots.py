import matplotlib
import matplotlib.pyplot as plt
import numpy


def eventPlot(events, tickDist=None):
    geneOrder = events.sum(1).argsort()
    sampleOrder = numpy.lexsort(events[geneOrder])[::-1]

    plt.imshow(
        events[geneOrder][:, sampleOrder], aspect="auto", interpolation="none",
        cmap=matplotlib.colors.ListedColormap(["#dddddd", "#377eb8"]) , origin="lower")

    if tickDist is not None:
        plt.xticks(numpy.arange(0, events.shape[1], tickDist), [])
    
    plt.grid(ls="-", c="white", axis="x")
    plt.setp(plt.gca().spines.values(), color="#888888", alpha=0)

    for tic in plt.gca().xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    for tic in plt.gca().yaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    try:
        plt.yticks(range(len(events)), events.featureNames[geneOrder])
    except AttributeError:
        pass

    plt.gca().yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(
        numpy.linspace(-0.5, len(events) - 0.5, len(events) + 1)))
    plt.grid(ls="-", c="white", axis="y", lw=4, alpha=1, which="minor")
    plt.grid(ls="None", axis="y", which="major")
    for tic in plt.gca().yaxis.get_minor_ticks():
        tic.tick1On = tic.tick2On = False

    plt.gcf().patch.set_facecolor("white")
