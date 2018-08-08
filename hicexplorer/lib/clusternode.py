class ClusterNode():

    def __init__(self, pChromosome, pStart, pEnd, pValueLeft, pValueRight, pParent=None, pChildLeft=None, pChildRight=None,
                 pParentId=None, pId=None, pChildLeftId=None, pChildRightId=None):
        self.chromosome = pChromosome
        self.start = pStart
        self.end = pEnd
        self.parent = pParent
        self.childLeft = pChildLeft
        self.childRight = pChildRight
        self.valueRight = pValueRight
        self.valueLeft = pValueLeft
        self.parentId = pParentId
        self.id = pId
        self.childLeftId = pChildLeftId
        self.childRightId = pChildRightId
