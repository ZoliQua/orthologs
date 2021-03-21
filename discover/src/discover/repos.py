import numpy
import os
import tables

import biotk.la


# From http://blog.genforma.com/wp-content/uploads/2011/05/baypiggies_talk_20110526.pdf
#def require_methods(*method_args):
#    """Class decorator to require methods on a subclass.
#    
#    Example usage
#    ------------
#    @require_methods('m1', 'm2')
#    class C(object):
#        'This class cannot be instantiated unless the subclass defines m1() and m2().'
#        def __init__(self):
#            pass
#    """
#    def fn(cls):
#        orig_init = cls.__init__
#        def init_wrapper(self, *args, **kwargs):
#            for method in method_args:
#                if (not (method in dir(self))) or \
#                   (not callable(getattr(self, method))):
#                    raise Exception("Required method %s not implemented" % method)
#            orig_init(self, *args, **kwargs)
#        cls.__init__ = init_wrapper
#        return cls
#    return fn


def require_methods(*methods):
    def wrapped(cls):
        cls.__require_methods__ = methods
        return cls
    return wrapped


class MixinContainerMetaClass(type):

    def __init__(self, name, bases, dict):
        attributes = reduce(set.union, map(dir, bases), set(dict))
        for base in bases:
            if hasattr(base, "__require_methods__"):
                for method in getattr(base, "__require_methods__"):
                    if method not in attributes:
                        raise Exception("Class %s does not implement %s" % (name, method))
        super(MixinContainerMetaClass, self).__init__(name, bases, dict)


class DataRepository(object):

    __metaclass__ = MixinContainerMetaClass

    def __init__(self, dataDir, filename="repos.h5"):
        self.data = tables.openFile(os.path.join(dataDir, filename), "a")

    def close(self):
        self.data.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()


class LabelledArrayMixin(object):

    def _storeLabelledArray(self, data, groupName, identifier, attrs={}):
        group = self.data.createGroup(groupName, identifier, createparents=True)
        for key, value in attrs.iteritems():
            group._v_attrs[key] = value
        self.data.createArray(group, "data", data)
        self.data.createArray(group, "featureNames", data.featureNames)
        self.data.createArray(group, "sampleNames", data.sampleNames)

    def _loadLabelledArray(self, groupName, identifier):
        group = self.data.getNode(groupName, identifier)

        if "subsetOf" in group._v_attrs:
            data = self._loadLabelledArray(groupName, group._v_attrs["subsetOf"])
            from biotk.util.join import match

            for node in group._f_iterNodes("Leaf"):
                try:
                    if node._v_attrs["subsetFilter"]:
                        axis = data.dimNames.index(node._v_name)
                        indices = match(node[:], data.labels[axis])
                        data = data.take(indices, axis)
                except KeyError:
                    pass
            
            return data
        else:
            return biotk.la.LabelledArray(
                group.data[:],
                [group.featureNames[:], group.sampleNames[:]],
                ["featureNames", "sampleNames"])

    def _defineSubset(self, groupName, identifier, dataset, **filters):
        group = self.data.createGroup(groupName, identifier, createparents=True)
        group._v_attrs["subsetOf"] = dataset
        for name, values in filters.iteritems():
            filter = self.data.createArray(group, name, numpy.asarray(values))
            filter._v_attrs["subsetFilter"] = True


@require_methods("_storeLabelledArray", "_loadLabelledArray", "_defineSubset")
class ExpressionDataMixin(object):

    def importExpressionData(self, exprs, identifier, chip):
        self._storeLabelledArray(exprs, "/data/exprs", identifier, {"chip": chip})

    def loadExpressionData(self, identifier):
        exprs = self._loadLabelledArray("/data/exprs", identifier)
        return exprs

    def defineExpressionSubset(self, identifier, dataset, **filters):
        self._defineSubset("/data/exprs", identifier, dataset, **filters)

    def listExpressionDataSets(self):
        return self.data.root.data.exprs._v_children.keys()

    #def loadExpressionData(self, identifier, entrezIds=False):
    #    group = getattr(self.data.root.data.exprs, identifier)
    #
    #    if "subsetOf" in group._v_attrs:
    #        exprs = self.loadExpressionData(group._v_attrs["subsetOf"], entrezIds)
    #        return exprs[:, match(group.sampleNames[:], exprs.sampleNames)]
    #
    #    else:
    #        featureNames = group.probeNames[:]
    #
    #        if entrezIds:
    #            probe2entrez = self.getChipAnnotationTable(group._v_attrs["chip"])
    #            i = probe2entrez["probe"].searchsorted(featureNames)
    #            probe2entrez = probe2entrez[i]
    #            featureNames = numpy.where(featureNames == probe2entrez["probe"],
    #                                       probe2entrez["entrezgene"], "NA")
    #            assert not numpy.any(featureNames == "NA")
    #
    #        return biotk.la.LabelledArray(
    #            group.exprs[:],
    #            [featureNames, group.sampleNames[:]],
    #            ["featureNames", "sampleNames"])


@require_methods("_storeLabelledArray", "_loadLabelledArray", "_defineSubset")
class CopyNumberDataMixin(object):

    def importCopyNumberData(self, cn, identifier, chip):
        self._storeLabelledArray(cn, "/data/cn", identifier, {"chip": chip})

    def loadCopyNumberData(self, identifier):
        exprs = self._loadLabelledArray("/data/cn", identifier)
        return exprs

    def defineCopyNumberSubset(self, identifier, dataset, **filters):
        self._defineSubset("/data/cn", identifier, dataset, **filters)

    def listCopyNumberDataSets(self):
        return self.data.root.data.cn._v_children.keys()
