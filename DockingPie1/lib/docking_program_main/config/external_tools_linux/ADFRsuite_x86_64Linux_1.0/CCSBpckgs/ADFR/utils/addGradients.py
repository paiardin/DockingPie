import os, sys, numpy, shutil
from ADFR.utils.maps import MapsFile
from Volume.IO.AutoGridWriter import WriteAutoGrid
from Volume.Grid3D import Grid3DF
from ADFRcc.adfrcc import AddGradients 
from ADFRcc.adfr import GridMap
from time import time
from threading import Thread

def addGradientToMaps(mapFiles, mapTypes, spacing, cutoffvalue, errorCut=0.01, neighborPts=1, logFileName=None):
    _map = GridMap()
    mapsFolder = os.path.split(mapFiles[0])[0]
    _map.loadFromMapFile("e", mapsFolder, 'rigidReceptor.%s.map'%"e")
    _emapdata = _map.getGridDataPy()
    _emin = min(_emapdata.flatten())
    _emax = max(_emapdata.flatten())
    _gradLow = int(-_emin)+1
    _writer = WriteAutoGrid()
    #print 'processing',
    maps = []
    for atype in mapTypes:
        if atype in ['e', 'd', 'sd']: continue
        #print atype,
        #sys.stdout.flush()
        map_ = GridMap()
        mapname = 'rigidReceptor.%s.map'%atype
        map_.loadFromMapFile(atype, mapsFolder, mapname)
        maps.append(map_)
    # C++ class
    addgrad = AddGradients(_gradLow)
    t0= time()
    print "processing maps ...",
    
    thread = Thread(target=addgrad.processMaps,
                    args = (maps, (spacing, spacing, spacing), 
                            neighborPts, cutoffvalue, errorCut, logFileName))
    thread.start()
    #we wait for the thread to finish
    #while (True):
    #    if not thread.isAlive():
    #        break
    #    thread.join(0.1)
    thread.join()
    #addgrad.processMaps(maps, (spacing, spacing, spacing), neighborPts, cutoffvalue, errorCut)
    print "done", time()-t0
    print "writing maps ...",
    t1 = time()
    for map_ in maps:
        _mapdata = map_.getGridDataPy()
        origin = map_.getOriginPy().astype('f')
        atype = map_.getMapType()
        mapname = 'rigidReceptor.%s.map'%atype
        grid = Grid3DF(_mapdata.astype("f"), origin, (spacing, spacing, spacing),
                       {'GRID_PARAMETER_FILE':'None',
                        'GRID_DATA_FILE':'None',
                        'MACROMOLECULE':'None'})
        _writer.write(grid, os.path.join(mapsFolder, mapname))
    print "done", time()-t1

def addGradientToMapsOld(mapFiles, mapTypes, spacing, cutoffvalue, errorCut=0.01):
    _map = GridMap()
    #print "mapFiles:", mapFiles
    mapsFolder = os.path.split(mapFiles[0])[0]
    _map.loadFromMapFile("e", mapsFolder, 'rigidReceptor.%s.map'%"e")
    _emapdata = _map.getGridDataPy()
    _emin = min(_emapdata.flatten())
    _emax = max(_emapdata.flatten())
    _gradLow = int(-_emin)+1
    _writer = WriteAutoGrid()
    #print 'processing', 
    for atype in mapTypes:
        if atype in ['e', 'd', 'sd']: continue
        #print atype,
        #sys.stdout.flush()
        map_ = GridMap()
        mapname = 'rigidReceptor.%s.map'%atype
        map_.loadFromMapFile(atype, mapsFolder, mapname)
        origin = map_.getOriginPy().astype('f')
        _mapdata = map_.getGridDataPy()

        ## find all grid points with negative values
        i,j,k = numpy.where(_mapdata<=_gradLow)
        indices = numpy.zeros((len(i),3),'i')
        indices[:,0] =  i
        indices[:,1] =  j
        indices[:,2] =  k
        neighborPts = 1
        # C++ class
        addgrad = AddGradients(_gradLow)
        addgrad.processMap(_mapdata, indices, (spacing, spacing, spacing), neighborPts, cutoffvalue, errorCut)
        #map_.setGridData(_mapdata)
        #map_.saveToMapFile(mapsFolder, mapname)
        grid = Grid3DF(_mapdata.astype("f"), origin, (spacing, spacing, spacing),
                       {'GRID_PARAMETER_FILE':'None',
                        'GRID_DATA_FILE':'None',
                        'MACROMOLECULE':'None'})
        _writer.write(grid, os.path.join(mapsFolder, mapname))
