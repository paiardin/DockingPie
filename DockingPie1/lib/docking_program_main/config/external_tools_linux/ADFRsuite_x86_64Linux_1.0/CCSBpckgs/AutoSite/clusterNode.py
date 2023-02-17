import numpy,os

class clusterNode(object):
    def __init__(self, id, size, gen=0, rg=100.0, buriedness=0.0, score=0.0, totalE=0.0, children=None):
        self.id=id
        self.size = size
        self.children = []
        self.rg=rg
        self.score=score
        self.totalE=totalE
        self.buriedness=buriedness
        self.generation=gen
        if children is not None:
            for child in children:
                self.add_child(child)

    def __repr__(self):
        return self.size

    def add_child(self, node):
        assert isinstance(node, clusterNode)
        self.children.append(node)

    # get leaf level node index
    def getNodeIndex(self):
        resultlist=""
        if not self.children:
            return " "+str(self.id)
        for child in self.children:
            resultlist=resultlist+child.getNodeIndex()
        return resultlist

    
    def getNodebySize(self,sizecutoff):
        resultlist=[]
        if self.size<sizecutoff or len(self.children)==0:
            resultlist.append(self)
            return resultlist
        if self.size<1.2*sizecutoff:
            for node in self.children:
                if node.size<0.8*sizecutoff:
                    resultlist.append(self)
                    return resultlist
        for node in self.children:
            #import pdb; pdb.set_trace()
            resultlist=resultlist+node.getNodebySize(sizecutoff)
            
        return resultlist
    
    def getAllNodes(self):
        resultlist=[]
        if len(self.children)==0:
            resultlist.append(self)
            #print 1
            return resultlist
        for node in self.children:
            resultlist=resultlist+node.getAllNodes()
        if self.size!=999999:
            resultlist.append(self)
            #import pdb;pdb.set_trace()
        #print len(resultlist)
        return resultlist


    def updateAllNodes(self,clProp):
        #resultlist=[]
        if len(self.children)==0:
            self.buriedness = clProp[self.id-1][4]
            self.rg = clProp[self.id-1][3]
            self.size = clProp[self.id-1][1]
            self.score = clProp[self.id-1][5]
            self.totalE=clProp[self.id-1][0]
            #resultlist.append(self)
            #print 1
            return 
        for node in self.children:
           node.updateAllNodes(clProp)
        if self.size!=999999:
            #resultlist.append(self)
            self.buriedness = clProp[self.id-1][4]
            self.rg = clProp[self.id-1][3]
            self.size = clProp[self.id-1][1]
            self.score = clProp[self.id-1][5]
            self.totalE=clProp[self.id-1][0]
            #import pdb;pdb.set_trace()
        #print len(resultlist)
        #return resultlist


    def writeJSON(self):
        d={'id':self.id,'size':self.size,'score':self.score,'totalE':round(self.totalE,3),'rg':self.rg,'buriedness':self.buriedness}
        #d={'id':self.id,'size':self.size,'score':self.score}
        d['children']=[child.writeJSON() for child in self.children]
        return d


def buildJSON(obj):
    node = clusterNode(id=obj['id'], size=obj['size'],score=obj['score'],totalE=obj['totalE'],rg=obj['rg'],buriedness=obj['buriedness'])
    for child in obj.get('children',[]):
        node.add_child(buildJSON(child))
    return node
