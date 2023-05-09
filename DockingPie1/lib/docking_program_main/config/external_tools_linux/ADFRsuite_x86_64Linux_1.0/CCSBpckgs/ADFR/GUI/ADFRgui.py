################################################################################
##
## This library is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## 
## This library is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
##
## (C) Copyrights Dr. Michel F. Sanner and TSRI 2016
##
################################################################################

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner and TSRI 2016
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/GUI/ADFRgui.py,v 1.5 2016/12/08 18:39:50 sanner Exp $
#
# $Id: ADFRgui.py,v 1.5 2016/12/08 18:39:50 sanner Exp $
#
import os, thread, numpy
from math import ceil, sqrt

from ADFRcc.adfr import RigidReceptorScorer, GridMap, Parameters
parameters = Parameters.getParameters()
feCoeffVdw = parameters.feCoeffVdw
feCoeffHbond = parameters.feCoeffHbond
feCoeffEstat = parameters.feCoeffEstat
feCoeffDesolv = parameters.feCoeffDesolv
feCoeffTors = parameters.feCoeffTors

from PySide import QtCore, QtGui
from ADFR import ADFR
from ADFR.GUI import ICONPATH
from ADFR.utils.runADFR import runADFR
from ADFR.utils.cluster import clusterPoses, oneCluster
from MolKit2 import Read
from MolKit2.molecule import getAtomIndicesPerType
from mglutil.math.rmsd import HungarianMatchingRMSD_prody
from mglutil.util.callback import CallbackFunction

class TableWidgetItem(QtGui.QTableWidgetItem):

    def __lt__(self, other):
        v1 = float(self.data(QtCore.Qt.DisplayRole))
        v2 = float(other.data(QtCore.Qt.DisplayRole))
        return (v1 < v2)

class GARunsMap(QtGui.QWidget):

    def __init__(self, dockButton):
        super(GARunsMap, self).__init__()
        self.dockButton = dockButton
        self.jobsStatus = None
        self.dx = 0
        self.dy = 0
        self.nrow = 0
        self.ncol = 0
        self.minSize = 8 # smallest width of a GA run patch
        self.setStyleSheet("background-color: white")
        #self.setMinimumSize(QtCore.QSize(202, self.minSize+2))
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Ignored, QtGui.QSizePolicy.Fixed)
        self.setSizePolicy(sizePolicy)

    def _getSize(self):
        # calculate ncol, nrow, dx, dy
        size = self.size()
        w = size.width()-2
        h = size.height()-2
        dy = self.minSize
        maxi = len(self.jobsStatus)
        if maxi*self.minSize < w: # single line
            nrow = 1
            ncol = maxi
            dx = (w)/ncol
            #dy = h
        else:
            dx = self.minSize
            ncol = (w)/self.minSize
            nrow = 1+maxi/ncol
            #dy = (h)/nrow
        self.dx = dx
        self.dy = dy
        self.nrow = nrow
        self.ncol = ncol
        
    def setJobs(self, jobsStatus):
        self.jobsStatus = jobsStatus
        self._getSize()
        self.setMinimumSize(QtCore.QSize(min(200, self.dockButton.size().width()), self.dy*self.nrow+2))
        self.resize(QtCore.QSize(self.minSize*self.ncol+2, self.dy*self.nrow+2))
        
    def paintEvent(self, e):

        qp = QtGui.QPainter()
        qp.begin(self)
        self.drawRectangles(qp)
        qp.end()

    def drawRectangles(self, qp):
        size = self.size()
        w = size.width()-2
        h = size.height()-2
        qp.eraseRect(0, 0, w, h)
        if self.jobsStatus is None: return
        
        maxi = len(self.jobsStatus)
        nb = 0
        for i in range(self.nrow):
            for j in range(self.ncol):
                color = QtGui.QColor(0, 0, 0) # outline color
                color.setNamedColor('#d4d4d4')
                qp.setPen(color)
                if self.jobsStatus[nb]==0: # not sarted
                    qp.setBrush(QtGui.QColor(200, 200, 200)) # rectangle color
                elif self.jobsStatus[nb]==1: # running
                    qp.setBrush(QtGui.QColor(200, 200, 0)) # rectangle color
                elif self.jobsStatus[nb]==2: # done
                    qp.setBrush(QtGui.QColor(0, 200, 0)) # rectangle color
                elif self.jobsStatus[nb]==3: # Failed
                    qp.setBrush(QtGui.QColor(200, 0, 0)) # rectangle color
                
                qp.drawRect(1+(j*self.dx), 1+(i*self.dy), self.dx, self.dy-1)
                nb += 1
                if nb == maxi:
                    return

class DockingClustersStackedHistograms(QtGui.QWidget):
    """
    draw a histograms with a bar for each cluster. The basr is broken into sub-bars
    for energy buckets.
    """
    def __init__(self):
        super(DockingClustersStackedHistograms, self).__init__()
        self.eBinWidth = 0.5
        self.colorsRGB = [ [0, 196,0, 128], [244, 232,0,128], [244,138,0,128], [242, 0, 0, 128]]
        #self.colorsRGB = [ [51,132,96,128], [86,199,62,128], [157,188,57,128],
        #                   [217,166,51,128], [223,92,40,128], [210,79,89,128] ]
        #[51,132,96],[216,105,213],[86,199,62],[119,124,223],[157,188,57],
        #                   [150,97,175],[217,166,51],[140,151,215],[223,92,40],[82,176,216],
        #                   [183,102,55],[57,186,178],[214,65,135],[91,198,124],[210,79,89],
        #                   [77,136,54],[224,135,190],[133,123,38],[75,112,148],[159,91,128] ]
        #self.colorsHEX = ["#338460", "#D869D5","#56C73E","#777CDF","#9DBC39",
        #               "#9661AF", "#D9A633","#8C97D7","#DF5C28","#52B0D8",
        #               "#B76637","#39BAB2","#D64187","#5BC67C","#D24F59",
        #               "#4D8836", "#E087BE","#857B26","#4B7094","#9F5B80"]
        self.maxBarHeight = 0
        self.barWidth = 8
        self.bestE = 0
        self.eRange = 0
        self.eBinWidth = 1.0
        self.clusters = None
        self.scores = None
        self.setMinimumSize(QtCore.QSize(10, 100))
        self.setStyleSheet("background-color: white")
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
        self.setSizePolicy(sizePolicy)

    def setMaxBarHeight(self, value):
        self.maxBarHeight = value
        
    def setClusters(self, clusters, scores):
        """
        len(clusters) is the number of clusters and clusters[0] is a list of
        solution indices with scrore within the first bucket
        """
        if clusters is None:
            self.clusters = None
            self.scores = None
            return
        
        # find best energy and biggest cluster
        eBest = 999999
        clbestE = []
        for cl in clusters:
            clbestE.append(scores[cl[0]])
            if scores[cl[0]] < eBest:
                eBest = scores[cl[0]]
            if len(cl) > self.maxBarHeight:
                self.maxBarHeight = len(cl)

        #compute frequencies of 1Kcal energy bins
        nbEbins = int(max(scores)-eBest)+1
        freq = []
        for i in range(len(clusters)):
            freq.append([0]*nbEbins) # set freq to 0 for each cluster and each ebin

        #print 'Clusters', clusters
        #print 'Scores', scores
        for cnum, cl in enumerate(clusters):
            #print cnum, cl
            for sol in cl:
                #print sol, int(scores[sol] - eBest)
                freq[cnum][int(scores[sol] - eBest)] += 1

        self.clusters = clusters
        self.scores = scores
        self.freq = freq
        #print 'Freq', self.freq
        self.bestE = eBest
        self.eRange = max(clbestE)-eBest

    def paintEvent(self, e):
        qp = QtGui.QPainter()
        qp.begin(self)
        self.drawRectangles(qp)
        qp.end()

    def drawRectangles(self, qp):
        size = self.size()
        w = size.width() - 20 # X labels
        h = size.height() - 20 # 20 used for labels
        qp.eraseRect(20, 0, w,h)
        if self.clusters is None:
            return

        # draw Y axis labels
        if self.maxBarHeight > 20:
            step = 5
        else:
            step = 1
        for i in range(0, self.maxBarHeight, step):
            qp.drawText(1, h-int(i*h/self.maxBarHeight), str(i))

        # draw X axis labels
        # first label is best energy
        qp.drawText(21, h+13,'%.3f'%self.bestE)
        if self.eRange == 0.0:
            pixelsPerKcal = float(w-2)/2
            tickSpacing = 1.
        elif self.eRange < 1.0:
            tickSpacing = .5 # one label every 1Kcal
            pixelsPerKcal = float(w-2)/ceil(self.eRange) # energy per pixel
        elif self.eRange < 4.0:
            tickSpacing = 1.0 # one label every 1Kcal
            pixelsPerKcal = float(w-2)/ceil(self.eRange) # energy per pixel
        else:
            tickSpacing = 3.
            pixelsPerKcal = float(w-2)/ceil(self.eRange) 

        for i in numpy.arange(tickSpacing, int(ceil(self.eRange/tickSpacing))+2, tickSpacing):
            qp.drawText(21+i*pixelsPerKcal, h+13,'+%.1f'%i)
            
        #barWidth = (w-2)/len(self.clusters)
        colqt = QtGui.QColor(0, 0, 0) # outline color
        #colqt.setNamedColor('#d4d4d4')
        qp.setPen(colqt)
        nbCol = len(self.colorsRGB)
        for cnum, frequencies in enumerate(self.freq): # loop over clusters
            if self.eRange == 0.0:
                x1 = 21 + (w-2-2*self.barWidth)
            else:
                x1 = 21 + (w-2-2*self.barWidth)*(self.scores[self.clusters[cnum][0]]-self.bestE)/self.eRange
            #print cl[0][0],self.scores[cl[0][0]], self.bestE, self.scores[cl[0][0]]-self.bestE, self.eRange, x1
            yoff = 0
            for j, freq in enumerate(frequencies):
                colorIndex = min(j, nbCol-1)
                color = self.colorsRGB[colorIndex]
                qp.setBrush(QtGui.QColor(*color)) # rectangle color
                barHeight = int(freq*h/self.maxBarHeight)
                #print x1, h-yoff-barHeight, self.barWidth, barHeight
                qp.drawRect(x1, h-yoff-barHeight, self.barWidth, barHeight)
                yoff += barHeight

class runGAThread(QtCore.QThread):
    endGASignal = QtCore.Signal(int, str, float, str, str)
    startGASignal = QtCore.Signal(int, str)

    def __init__(self, nbRuns, parent=None):
        super(runGAThread, self).__init__(parent)
        self.todo = nbRuns
        self.done = 0

    def gaStart_cb(self, jobNum, logFile):
        self.startGASignal.emit(jobNum, logFile)
        
    def gaDone_cb(self, jobNum, logFile, percent, status, error):
        #print 'in thread', jobNum, logFile, percent, status, error
        self.endGASignal.emit(jobNum, logFile, percent, status, error)

    def run(self, args):
        runADFR(args, cb_end=self.gaDone_cb, cb_start=self.gaStart_cb)

class MyQLineEdit(QtGui.QLineEdit):
    
    def __init__(self, title, filters, parent=None):
        super(MyQLineEdit, self).__init__(parent)
        self.title = title
        self.filters = filters
        
    def mouseDoubleClickEvent(self, event):
        self.setStyleSheet("background-color: None")
        currentText = self.text()
        filename = self.getLigandFilename(event)
        if filename:
            self.setText(filename)
        else:
            self.clear()
            self.setText(currentText)
            
    def getLigandFilename(self, event):
        filename, selfilter = QtGui.QFileDialog.getOpenFileName(
            self, self.tr(self.title), '', self.tr(self.filters))
        return filename.encode('ascii', 'replace')
        
class SingleDockingInputWidget(QtGui.QWidget):

    def __init__(self, pmvViewer=None, parent=None):
        super(SingleDockingInputWidget, self).__init__(parent)
        self.pmvViewer = pmvViewer
        self.ligandFilename = None
        self.mapsFilename = None
        self.outputFilename = None
        self.unzippedMapsFolder = None
        self.receptor = None
        self._scores = []
        self._genes = []
        self._energies = [] # list of energies fromt eh GAs
        self._rmsdsRef = [] # list of RMSDs to refernce structure if provided
        self.dockedLigand = None # prody molecule with multiple coordinate sets
        self.minimizedLigand = None # prody molecule with multiple coordinate sets
        if pmvViewer:
            from DejaVu2.IndexedPolylines import IndexedPolylines
            self.deltaLines = IndexedPolylines('deviations', visible = 0)
            pmvViewer.AddObject(self.deltaLines)
        self.buildUI()
            
    def checkReady(self):
        self.dockButton.setDisabled(False)
        return True
        ## if self.ligandFilename is None or \
        ##    self.mapsFilename is None or \
        ##    self.outputFilename is None or len(self.outputFilename)==0:
        ##     self.dockButton.setDisabled(True)
        ##     return False
        ## else:
        ##     self.dockButton.setDisabled(False)
        ##     return True
        # build a scorer to get energy breakdown

    def makeScorer(self, flexRes=None):
        from ADFR.utils.scorer import ADFRscorer
        folder = self.unzippedMapsFolder.encode('ascii', 'replace')
        adfr, ind, _score = ADFRscorer(self.dockedLigand, folder, 'receptor',
                                       flexRes=flexRes)
        self._adfr = adfr
        self._ind = ind
        self._pop = adfr.createPopulation(1)
        from ADFR.GA import GA
        self._ga = GA(self._pop, self.dockedLigand)

    def minimize(self):
        deltas = []
        v = self._ind.phenotype[1].tolist()
        nv = len(v)
        f = numpy.arange(nv)
        faces = []
        for i in range(10):
            self._ind.setGenes(self._ind.genomePy.getIdentityGenesPy())
            s0 = self._ind.score()
            sc, nb = self._ga.minimize(self._ind, nbSteps=100, noImproveStop=20,
                                       max_steps=3000, MAX_FAIL=100, searchRate=.1)
            if self.pmvViewer:
                #self.minimizedLigand.geomContainer.allCoords[:] = self._ind.phenotye[1]
                #self.pmvViewer.pmv.displayLines(self.minimizedLigand)
                v1 = self._ind.phenotype[1]
                delta = v1 - self.dockedLigand.geomContainer.allCoords[:]
                #self.dockedLigand.geomContainer.allCoords[:] = v1
                #self.pmvViewer.pmv.displayLines(self.dockedLigand)
                v.extend(v1)
                for a in range(nv):
                    faces.append( (a, a+(i+1)*nv) )
                print -s0, -self._ind._score
        self.deltaLines.Set(visible=True, vertices = v, faces=faces)
        
    #def getLigandFilename(self):
    #    filename, selfilter = QtGui.QFileDialog.getOpenFileName(
    #        self, self.tr("Read Ligand"), '',
    #        self.tr("PDBQT Files (*.pdbqt);; All files (*)"))
    #    return filename.encode('ascii', 'replace')

    #def getReceptorMapsFilename(self):
    #    filename, selfilter = QtGui.QFileDialog().getOpenFileName(
    #        self, self.tr("Read Receptor Maps"), '',
    #        self.tr("zip Files (*.zip);; All files (*)"))
    #    return filename.encode('ascii', 'replace')

    def getLigand(self, filename):
        if os.path.exists(filename):
            self.ligandEntryWidget.setStyleSheet("background-color: None")
            mol = Read(filename.encode('ascii', 'replace'))
            self.setLigand(mol)
            self.checkReady()
            if self.unzippedMapsFolder is not None:
                self.makeScorer()
        else:
            self.ligandEntryWidget.setStyleSheet("background-color: #F14D81")

    def setLigand(self, mol):
        if self.dockedLigand:
            self.pmvViewer.pmv.deleteMolecule(self.dockedLigand)
        self.dockedLigand = mol
        atoms = mol.select()
        d1 = getAtomIndicesPerType(atoms)
        self.rmsdCalc = HungarianMatchingRMSD_prody(atoms.getCoords(), d1, d1)
        
        if self.pmvViewer:
            pmv = self.pmvViewer.pmv
            pmv.addMolecule(mol)
            pmv.customColor(mol.select('element C'), [(0.,1.,1.)], geomsToColor=['lines'])
            #pmv.displaySticksAndBalls(mol)
            if len(pmv.Mols)==1:
                self.pmvViewer.Reset_cb()
                self.pmvViewer.Normalize_cb()
                self.pmvViewer.Center_cb()
            
    def setGridVisible(self, value):
        # value is 0 for unchecked and 2 for checked for checkbox
        # not(value==0) make it work for 0, 1, 2, False, True
        self.boxGeom.master.Set(visible = not(value==0))
        for c in self.boxGeom.master.children:
            if c.name=='faces':
                c.Set(visible = 0)
            else:
                c.Set(visible = not(value==0))

    def getMaps(self, filename):
        if os.path.exists(filename):
            from ADFR.utils.maps import MapsFile
            self.mf = mf = MapsFile(filename)
            mf.unzipMaps()
            self.unzippedMapsFolder = unzippedMapsFolder = mf.getMapsFolder()
            receptorFilename = os.path.join(mf.getMapsFolder(),
                                        mf.getReceptorFilename())
            flexRes = mf.getFlexResStr()
            #flexResStr = mf.getFlexResStr()
            #from ADFR.utils.maps import flexResStr2flexRes
            #flexRes = flexResStr2flexRes(flexResStr)
            covalentRec = mf.getCovalentBond()
            if covalentRec is not None:
                covalentRec.insert(
                    0, int(mf.getCovalentBondTorsionAtom().split()[1][1:-1]))

            self.mapsFilename = filename
            self.checkReady()
            if self.receptor and self.pmvViewer:
                self.pmvViewer.pmv.deleteMolecule([self.receptor])
            self.receptor = Read(receptorFilename)
            #if self.dockedLigand is not None:
            #    self.makeScorer(flexRes=flexRes)
            if self.pmvViewer:
                self.pmvViewer.pmv.addMolecule(self.receptor)

                from DejaVu2.Box import NiceBox
                b = self.boxGeom = NiceBox('gridOutline')
                b.setCenter(*mf.getBoxCenter())
                b.setSides(*mf.getBoxSize())
                self.boxGeom.addToViewer(self.pmvViewer)
                self.setGridVisible(True)

                ## from DejaVu2.Points import Points
                ## self.TPoints = Points(
                ##     'tpoints', visible=1, inheritMaterial=False,
                ##     materials=[(1,0,0)], inheritPointWidth=False,
                ##     pointWidth=4.)
                ## self.pmvViewer.AddObject(self.TPoints)

                #from DejaVu2.Spheres import Spheres
                #self.anchorAtomGeom = Spheres('rootAtom', visible=0, inheritMaterial=False,
                #                               materials=[(1,0,1)], inheritFrontPolyMode=False,
                #                              frontPolyMode='line', quality=2,
                #                              inheritLineWidth=0, lineWidth=1)
                #self.pmvViewer.AddObject(self.anchorAtomGeom)
            

    def setOutput(self, text):
        self.outputFilename = text.encode('ascii', 'replace')
        self.checkReady()
        
    def gaStart_cb(self, jobNum, logFile):
        #print 'in main Start', jobNum, logFile, percent
        self._jobStatus[jobNum] = 1
        self.gaRunsMap.setJobs(self._jobStatus)
        self.gaRunsMap.update()
        
    def getPoseData(self, logFile):
        f = open(logFile)
        lines = f.readlines()
        f.close()
        w1 = lines[-3].split()
        w2 = lines[-2].split()
        return float(w1[2]), float(w1[4]),{
                'RRL': float(w2[1][:-1]), 'FRFR': float(w2[3][:-1]),
                'RRFR': float(w2[5][:-1]), 'wRR': float(w2[7][:-1]),
                'LL': float(w2[9][:-1]), 'FRL': float(w2[11][:-1])}

    def updateBestLabels(self, jobNum, score, rmsdRef, energies):
        self.bestScoreLabel.setText('job: %d score: %.3f'%(jobNum+1, score))

        if energies['FRFR'] != 0.0:
            lab = "LL: %.3f, RL: %.3f, 'FRL: %.3f, FRFR: %.3f, RRFR: %.3f"%(energies['LL'], energies['RRL'], energies['FRL'], energies['FRFR'], energies['RRFR'])
        else:
            lab = "LL: %.3f, RL: %.3f"%(energies['LL'], energies['RRL'])
        self.bestScoreEnergyLabel.setText(lab)

        self.rmsdCalc.setRefCoords(self.dockedLigand._ag._coords[self.best_score_jobnum])
        rmsdBest = self.rmsdCalc.computeRMSD(self.dockedLigand._ag._coords[jobNum])

        self.rmsdsLabel.setText('ref: %.3f solution: %.3f'%(rmsdRef, rmsdBest))

    def gaDone_cb(self, jobNum, logFile, percent, status, error):
        #print 'in main end', jobNum, logFile, percent, status, error

        if status=='OK':
            self._jobStatus[jobNum] = 2
            self.gaRunsMap.setJobs(self._jobStatus)
            self.gaRunsMap.update()
            score, rmsdRef, energies = self.getPoseData(logFile)
            self._scores[jobNum] =  score
            self._rmsdsRef[jobNum] =  rmsdRef
            self._energies[jobNum] =  energies

            # get pose coordinates
            ligandFilename = '%s_%04d_lig.pdbqt'%(self.outputNameWidget.text(), jobNum)
            lig = Read(ligandFilename)
            ag = self.dockedLigand._ag
            ag.setACSIndex(jobNum)
            ag.setCoords(lig._ag.getCoords())

            # get ligand genes
            f = open(ligandFilename)
            lines = f.readlines()
            f.close()
            ln = 3
            words = lines[ln].split()
            if words[1]=='GENES':
                nbGenesLines = int(words[2])
                genes = []
                for i in range(nbGenesLines):
                    words = lines[ln+1+i].split('|==|')
                    genes.extend([float(x) for x in words[1].split()])
                self._genes[jobNum] = genes
            else:
                print 'ERROR: GENES not found', lines[0]

            if score < self.best_score:
                if self.pmvViewer:
                    self.dockedLigand.geomContainer.allCoords[:] = lig._ag.getCoords()
                    self.pmvViewer.pmv.displayLines(self.dockedLigand)
            
                self.best_score = score
                self.best_score_jobnum = jobNum
                self.best_score_rmsdRef = rmsdRef
                self.best_score_energies = energies
                self.updateBestLabels(jobNum, score, rmsdRef, energies)

        elif status=='FAILED':
            #b.setStyleSheet("background-color: red")
            self._jobStatus[jobNum] = 3
            self.gaRunsMap.setJobs(self._jobStatus)
            self.gaRunsMap.update()
            print 'ERROR', error
            
        ## cluster solutions
        #order = numpy.argsort([x for x in self._scores if x is not None])
        order = []
        scores = [] # list of scores from the self._scores for jobs that have completed
        #build scores list and list of indices of solutions to be clustered
        for i, sc in enumerate(self._scores):
            if sc is not None:
                order.append(i) # because solution coords start at self.dockedLigand._ag._coords[1]
                scores.append(sc)
        # make sure the 'order' list is sorted by score
        oorder = numpy.argsort(scores)
        order = numpy.array(order)[oorder]
        if len(order)>1:
            # cluster all solutions
            #print 'ORDER', order
            #print 'scores', self._scores
            self.clusters = clusterPoses(self.dockedLigand._ag._coords, order,
                                         self.rmsdCalc, cutOff=2.0)
            #print 'clusters', self.clusters
            #for i, c in enumerate(self.clusters):
            #    print i, c, [self._scores[j] for j in c]
                
            self.gaRunsMap.setJobs(self._jobStatus)
            # bin scores in each cluster
            ## eBinWidth = 0.5
            ## minE = min(scores)
            ## maxE = max(scores)
            ## nBins = int(ceil((maxE-minE)/eBinWidth))
            ## #print 'NBINS', nBins, maxE, minE, eBinWidth
            ## #print 'energies', min(self._scores), max(self._scores)
            ## histo = [None]* len(self.clusters)
            ## for cnum, cl in enumerate(self.clusters):
            ##     count = [0]*nBins
            ##     for solInd in cl:
            ##         count[int((self._scores[solInd]-minE)/eBinWidth)] += 1
            ##     histo[cnum] = count
            #print 'HISTO', histo
            self.clustersWidget.setClusters(self.clusters, self._scores)
            self.clustersWidget.update()
            
        if percent==1.0:
            self.dockButton.setText('dock')
            if len(order)>1:
                self.setNbCusterButtons(len(self.clusters))

            self.setPose(self.best_score_jobnum)
            
    def setPose(self, i):
        if self.dockButton.text() == 'stop': return
        if self.pmvViewer:
            self.dockedLigand.geomContainer.allCoords[:] = self.dockedLigand._ag._coords[i]
            self.pmvViewer.pmv.displayLines(self.dockedLigand)
        self.updateBestLabels(i, self._scores[i], self._rmsdsRef[i], self._energies[i])

        self._ind.setGenes(self._genes[i])
        _score = self._ind.score()
        #print 'POSE', i, _score,
        
        from ADFR.utils.analyze import getHBPairs, addHBlines
        atoms =  self._adfr.ligandFT.mol.select()
        #hbPairs, hbEne = getHBPairs(self._ind, atoms, cutOffEne=-0.001)
        #if len(hbPairs):
        #    geoms = addHBlines(self.pmvViewer, hbPairs, hbEne, atoms.getCoords())
        #import pdb; pdb.set_trace()
        self.detailsWidget.fillTable(self._ind, self._adfr, self.dockedLigand._ag._coords[i])

    def setNbGA(self, num):
        self._jobStatus = [0]*num
        self.gaRunsMap.setJobs(self._jobStatus)
        self.gaRunsMap.update()
        self.setNbCusterButtons(0)
        self.clustersWidget.setMaxBarHeight(num)
        self.clustersWidget.setClusters(None, None)
        self.clustersWidget.update()

    def setNbCusterButtons(self, num):
        for b in self.clButtons:
            self.clButtonsLayout.removeWidget(b)
            b.setParent(None)
            b.deleteLater()
        self.clButtons = []
        n = 0
        nbPerRow = 3
        for i in range(num):
            self.rmsdCalc.setRefCoords(self.dockedLigand._ag._coords[self.best_score_jobnum])
            rmsdBest = self.rmsdCalc.computeRMSD(self.dockedLigand._ag._coords[self.clusters[i][0]])
            w = QtGui.QPushButton("%d (%.2f)"%(n+1,rmsdBest))
            w.setFixedSize(QtCore.QSize(50, 15))
            self.clButtons.append(w)
            cb = CallbackFunction(self.setPose, self.clusters[i][0])
            w.clicked.connect(cb)
            self.clButtonsLayout.addWidget(w, n/nbPerRow, n-nbPerRow*(n/nbPerRow))
            n += 1
        self.clButtonsLayout.update()
        #import pdb; pdb.set_trace()
        
    def runDocking(self, inThread=True):
        # reset buttons to default color
        #for b in self.buttons:
        #    b.setStyleSheet("background-color: None")

        # delete the cluster buttons
        self.setNbCusterButtons(0)

        self.best_score = 9999999999.
        self.best_score_jobnum = -1
        self.best_score_rmsd = -1
        self.best_score_energies = {}
        nbGA = self.gaNbWidget.value()
        self._scores =  [None]*nbGA
        self._genes =  [None]*nbGA
        self._rmsdsRef =  [None]*nbGA
        self._energies =  [None]*nbGA

        # makes sure we have enough coord sets to store poses
        # first job is ni coordinate set 1 NOT 0
        ag = self.dockedLigand._ag
        coords = ag.getCoords()
        if ag.numCoordsets() < nbGA:
            for i in range(ag.numCoordsets(), nbGA):
                self.dockedLigand._ag.addCoordset(coords, 'pose %d'%(i))

        self.bestScoreLabel.setText('job: %d score: %.3f'%(-1, 0.))
        self.bestScoreEnergyLabel.setText("")
        self.rmsdsLabel.setText('ref: %.3f solution: %.3f'%(-1, -1))
        
        args = [None, self.ligandEntryWidget.text().encode('ascii', 'replace'),
                '--target', '"%s"'%self.mapsEntryWidget.text().encode('ascii', 'replace'),
                '--jobName', '"%s"'%self.outputNameWidget.text(),
                '--maxCores', str(self.coreNbWidget.value()),
                '-o', '"%s"'%self.outputFilename,
                '-O',
                '--nbRuns', str(nbGA),
                '--maxEvals',  str(self.maxEvalsWidget.value()),
                ] # first agument is ignored
        refLig = self.refLigWidget.text()

        if refLig:
            args.append('-r')
            args.append(refLig)
        
        #print args
        #print ' '.join(args[1:])
    
        self.dockButton.setText('stop')
        #self.dockButton.setDisabled(True)
        nga = self.gaNbWidget.value()
        self._jobStatus = [0]*nga
        self.gaRunsMap.setJobs(self._jobStatus)
        self.gaRunsMap.update()
        #self.clustersWidget.setMaxBarHeight(nga)
        self.clustersWidget.setMaxBarHeight(1)
        self.clustersWidget.setClusters(None, None)
        self.clustersWidget.update()
        
        gaThread = runGAThread(nga)
        gaThread.startGASignal.connect(self.gaStart_cb,
                                       QtCore.Qt.QueuedConnection)
        gaThread.endGASignal.connect(self.gaDone_cb,
                                     QtCore.Qt.QueuedConnection)
        if inThread:
            thread.start_new_thread( gaThread.run, (args,) )
        else:
            gaThread.run(args) 
        
    def buildUI(self):
        layout = QtGui.QVBoxLayout()

        grp1 = QtGui.QGroupBox("input")
        formLayout = QtGui.QFormLayout()
        w = self.ligandEntryWidget = MyQLineEdit("Read Ligand", "PDBQT Files (*.pdbqt);; All files (*)")
        #w.textChanged.connect(self.checkReady)
        w.textChanged.connect(self.getLigand)
        formLayout.addRow(self.tr("ligand:"), self.ligandEntryWidget)

        w = self.mapsEntryWidget = QtGui.QLineEdit()
        #w.textChanged.connect(self.checkReady)
        w.textChanged.connect(self.getMaps)
        formLayout.addRow(self.tr("target:"), self.mapsEntryWidget)

        w = self.refLigWidget = QtGui.QLineEdit()
        formLayout.addRow(self.tr("reference ligand:"), self.refLigWidget)

        grp1.setLayout(formLayout)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        grp1.setSizePolicy(sizePolicy)
        layout.addWidget(grp1)
        #ret = layout.setStretch(grp1, 1)
        #print 'FUGU', layout.stretch(0)
        #print 'FUGU1', layout.stretch(1)
        
        grp2 = QtGui.QGroupBox("parameters")
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        grp2.setSizePolicy(sizePolicy)
        formLayout = QtGui.QFormLayout()

        self.outputNameWidget = QtGui.QLineEdit()
        self.outputNameWidget.textChanged.connect(self.setOutput)
        formLayout.addRow(self.tr("output name:"), self.outputNameWidget)
        
        import multiprocessing
        ncpu = multiprocessing.cpu_count()
        w = self.coreNbWidget = QtGui.QSpinBox()
        w.setValue(ncpu-1)
        w.setRange(1, ncpu)        
        formLayout.addRow(self.tr("cores:"), self.coreNbWidget)

        w = self.gaNbWidget = QtGui.QSpinBox()
        w.setValue(50)
        w.setMinimum(1)
        w.setMaximum(999999)
        w.valueChanged.connect(self.setNbGA)
        formLayout.addRow(self.tr("GA runs:"), self.gaNbWidget)
        grp2.setLayout(formLayout)

        w = self.maxEvalsWidget = QtGui.QSpinBox()
        w.setMinimum(1)
        w.setMaximum(99999999)
        w.setValue(5000000)
        formLayout.addRow(self.tr("max. evals.:"), self.maxEvalsWidget)
        grp2.setLayout(formLayout)
        layout.addWidget(grp2)

        w = self.minimizeButton = QtGui.QPushButton('minimize')
        w.clicked.connect(self.minimize)
        layout.addWidget(w)

        grp3 = QtGui.QGroupBox("run")
        gLayout = QtGui.QVBoxLayout()
        w = self.dockButton = QtGui.QPushButton('dock')
        w.setDisabled(True)
        gLayout.addWidget(w)
        self.dockButton.clicked.connect(self.runDocking)

        w = self.gaRunsMap = GARunsMap(self.dockButton)
        gLayout.addWidget(w)

        bestForm = QtGui.QFormLayout()
        self.bestScoreLabel = QtGui.QLabel('None')
        bestForm.addRow(self.tr("solution:"), self.bestScoreLabel)
        self.bestScoreEnergyLabel = QtGui.QLabel('-1')
        bestForm.addRow(self.tr("energies:"), self.bestScoreEnergyLabel)
        self.rmsdsLabel = QtGui.QLabel('None')
        bestForm.addRow(self.tr("RMSD:"), self.rmsdsLabel)
        gLayout.addLayout(bestForm)

        # add cluster histogram
        self.clustersWidget = w = DockingClustersStackedHistograms()

        # add color legend
        hlayout = QtGui.QHBoxLayout()
        colors = self.clustersWidget.colorsRGB
        for i in range(len(colors)):
            if i==len(colors)-1:
                l1 = QtGui.QLabel(">%dKcal"%(i+1))
            else:
                l1 = QtGui.QLabel("%dKcal"%(i+1))
            l1.setAlignment(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignHCenter)
            l1.setFrameStyle(QtGui.QFrame.StyledPanel | QtGui.QFrame.Plain)
            qcol = QtGui.QColor(*colors[i])
            l1.setStyleSheet("background-color: %s"%qcol.name())
            l1.setMinimumSize(QtCore.QSize(30, 15))
            hlayout.addWidget(l1)
        gLayout.addLayout(hlayout)

        gLayout.addWidget(w)

        self.clButtons = []
        self.clButtonsLayout = QtGui.QGridLayout()
        gLayout.addLayout(self.clButtonsLayout)
        
        grp3.setLayout(gLayout)
        layout.addWidget(grp3)
        self.setLayout(layout)

class SingleDockingDetailsWidget(QtGui.QWidget):

    def __init__(self, PmvViewer, parent=None):
        super(SingleDockingDetailsWidget, self).__init__(parent)
        self.PmvViewer = PmvViewer
        self.buildUI(parent)
        from DejaVu2.Spheres import Spheres
        self.LRSpheres = Spheres('Ligand-Receptor grid interactions', visible=False,
                                 inheritMaterial=False, transparent=True, opacity=0.4)
        PmvViewer.AddObject(self.LRSpheres)
        
    def buildUI(self, parent):
        self.tabWidget = QtGui.QTabWidget(parent)
        w = self.interactionsTableWidget = QtGui.QTableWidget(parent)
        w.setColumnCount(6)

        w.setHorizontalHeaderLabels(
            ["name", "element", "energy", "x", "y", "z"])
        w.itemSelectionChanged.connect(self.onSelectLR)
        
        self.recTableWidget = QtGui.QTableWidget(parent)

        w = self.ligTableWidget = QtGui.QTableWidget(parent)
        w.setColumnCount(8)
        w.setHorizontalHeaderLabels(
            ["non-bond", "at1-at2", "distance", "total",
             "elec", "vdW", "hb", "desolv"] )
        w.itemSelectionChanged.connect(self.onSelectLL)

        self.tabWidget.addTab(self.interactionsTableWidget, 'ligand-receptor')
        self.tabWidget.addTab(self.ligTableWidget, 'ligand')
        self.tabWidget.addTab(self.recTableWidget, 'receptor')
        self.layout = QtGui.QVBoxLayout()
        self.layout.addWidget(self.tabWidget)
        self.setLayout(self.layout)

    def onSelectLR(self):
        w = self.interactionsTableWidget
        indices = w.selectionModel().selectedIndexes()
        #print 'OnSelect', indices
        coords = []
        radii = []
        
        for index in sorted(indices):
            row = index.row()
            x = float(w.item(row, 3).text())
            y = float(w.item(row, 4).text())
            z = float(w.item(row, 5).text())
            coords.append( (x,y,z) )
            radii.append( -float(w.item(row, 2).text()) )

        if len(coords):
            self.LRSpheres.Set(visible=True, vertices=coords, radii=numpy.array(radii)/self.RRL_range)
        else:
            self.LRSpheres.Set(visible=False)
            
    def onSelectLL(self):
        w = self.ligTableWidget
        indices = w.selectionModel().selectedIndexes()
        
    def fillTable(self, ind, adfr, coords):
        scorer = ind.genomePy.scorer.getLlPairwiseScorer()
        energiesLL = scorer.getTotalScore()
        distance = scorer.getDistanceMatrix()
        eArray = scorer.getEstatEnergyMatrix()
        vdwArray = scorer.getVdwEnergyMatrix()
        hArray = scorer.getHbEnergyMatrix()
        dsArray = scorer.getSolvEnergyMatrix()
        atomSet = scorer.getAtomSet1().atomSetStatic

        table = self.ligTableWidget
        table.clear()
        table.setRowCount(0)
        table.setHorizontalHeaderLabels(
            ["non-bond", "at1-at2", "distance", "total",
             "elec", "vdW", "hb", "desolv"] )
        table.setColumnCount(8)
        #for i in range(table.rowCount()):
        #    table.removeRow(0)

        ligAtoms = adfr.ligandFT.mol.select()
        table.setSortingEnabled(False)
        n = 0
        for i, a1 in enumerate(ligAtoms):
            n1 = a1.getName()
            for j, a2 in enumerate(ligAtoms):
                if j<=i: continue
                n2 = a2.getName()
                if atomSet.getPairScorable(i,j):
                    e = eArray[i][j]*feCoeffEstat
                    v = vdwArray[i][j]*feCoeffVdw
                    h = hArray[i][j]*feCoeffHbond
                    d = dsArray[i][j]*feCoeffDesolv
                    table.insertRow(n)
                    newItem = QtGui.QTableWidgetItem("%d-%d"%(i+1, j+1))
                    table.setItem(n, 0, newItem)
                    newItem = QtGui.QTableWidgetItem("%s-%s"%(n1, n2))
                    table.setItem(n, 1, newItem)
                    newItem = TableWidgetItem("%7.3f"% distance[i][j])
                    table.setItem(n, 2, newItem)
                    newItem = TableWidgetItem("%7.4f"%(e+v+h+d))
                    table.setItem(n, 3, newItem)
                    newItem = TableWidgetItem("%7.4f"%(e))
                    table.setItem(n, 4, newItem)
                    newItem = TableWidgetItem("%7.4f"%(v))
                    table.setItem(n, 5, newItem)
                    newItem = TableWidgetItem("%7.4f"%(h))
                    table.setItem(n, 6, newItem)
                    newItem = TableWidgetItem("%7.4f"%(d))
                    table.setItem(n, 7, newItem)
                    n += 1
        table.resizeColumnsToContents()
        table.setSortingEnabled(True)
            
        table = self.interactionsTableWidget
        table.clear()
        table.setRowCount(0)
        table.setColumnCount(6)
        table.setHorizontalHeaderLabels(
            ["name", "element", "energy", "x", "y", "z"])
        ADelem = adfr.ligandFT.mol._ag._data['AD_element']
        RRL = ind.genomePy.scorer.getLrrGridScorer().getScoreArray()
        self.RRL_range = abs(max(RRL) - min(RRL))
        table.setSortingEnabled(False)
        for i, a1 in enumerate(ligAtoms):
            table.insertRow(i)
            x, y ,z = coords[i]
            newItem = QtGui.QTableWidgetItem("%4s"%(a1.getName()))
            table.setItem(i, 0, newItem)
            newItem = QtGui.QTableWidgetItem("%2s"%(ADelem[a1.getIndex()]))
            table.setItem(i, 1, newItem)
            newItem = TableWidgetItem("%9.3f"%(RRL[i]))
            table.setItem(i, 2, newItem)
            newItem = TableWidgetItem("%9.3f"%(x))
            table.setItem(i, 3, newItem)
            newItem = TableWidgetItem("%9.3f"%(y))
            table.setItem(i, 4, newItem)
            newItem = TableWidgetItem("%9.3f"%(z))
            table.setItem(i, 5, newItem)
        table.resizeColumnsToContents()
        table.setSortingEnabled(True)

class SingleDockingWidget(QtGui.QWidget):

    def __init__(self, pmvViewer=None, parent=None):
        super(SingleDockingWidget, self).__init__(parent)
        self.buildUI(pmvViewer, parent)

    def buildUI(self, pmvViewer=None, parent=None):
        mainLayout = QtGui.QHBoxLayout()
        splitter = QtGui.QSplitter()
        self.tabWidget = QtGui.QTabWidget(parent)
        self.inputWidget = SingleDockingInputWidget(pmvViewer, parent)
        self.tabWidget.addTab(self.inputWidget, 'docking')

        self.detailsWidget = SingleDockingDetailsWidget(pmvViewer, parent)
        self.inputWidget.detailsWidget = self.detailsWidget
        
        self.tabWidget.addTab(self.detailsWidget, 'details')
        splitter.addWidget(self.tabWidget)
        if pmvViewer:
            splitter.addWidget(pmvViewer.cameras[0])
        mainLayout.addWidget(splitter)
        self.setLayout(mainLayout)

        #self.gridLayout = QtGui.QGridLayout(self)
        #self.tabWidget = QtGui.QTabWidget(self)
        #self.inputWidget = SingleDockingInputWidget(pmvViewer, parent)
        #self.resultsWidget = SingleDockingResultsWidget(parent)
        #self.tabWidget.addTab(self.inputWidget, 'docking')
        #self.tabWidget.addTab(self.resultsWidget, 'details')
        #self.gridLayout.addWidget(self.tabWidget)
        

        
if __name__=='__main__':
    import sys

    #from MolKit2 import Read
    #rec = Read('4EK3_rec.pdbqt')
    #lig = Read('4EK4_lig.pdbqt')
    app = QtGui.QApplication(sys.argv)

    from PmvApp import mkPmvApp
    pmv = mkPmvApp()
    
    from PmvApp.GUI.Qt.PmvGUI import PmvViewer
    viewer = PmvViewer(pmv)
    viewer.cameras[0].setMinimumWidth(450)

    widget = SingleDockingWidget(pmvViewer=viewer)    
    widget.inputWidget.pmvViewer=viewer
    widget.resize(150,300)
    #widget.inputWidget.ligandEntryWidget.setText('4EK4_lig.pdbqt')
    #widget.inputWidget.mapsEntryWidget.setText('4EK3_rec.zip')
    #widget.inputWidget.refLigWidget.setText('4EK4_lig.pdbqt')

    #widget.inputWidget.ligandEntryWidget.setText('Astex/ligands/1jje_lig.pdbqt')
    #widget.inputWidget.mapsEntryWidget.setText('1jje_rec.zip')
    #widget.inputWidget.refLigWidget.setText('Astex/ligands/1jje_lig.pdbqt')

    #widget.inputWidget.ligandEntryWidget.setText('HSP90/1uy6_lig.pdbqt')
    #widget.inputWidget.mapsEntryWidget.setText('2wi2_bestRes.zip')
    #widget.inputWidget.refLigWidget.setText('/mnt/ark/home/sanner/AD5/versions/flexLoop/2wi2_bestRes.pdbqt')

    #widget.inputWidget.ligandEntryWidget.setText('ZikaPro/CP01187.pdbqt')
    widget.inputWidget.ligandEntryWidget.setText('ZikaPro/CP01194.pdbqt')
    widget.inputWidget.mapsEntryWidget.setText('ZikaPro/5lc0_A_cl1and2.zip')
    #widget.inputWidget.mapsEntryWidget.setText('ZikaPro/5lc0_A_full_allsites.zip')

    widget.inputWidget.outputNameWidget.setText('NoName')        
    widget.inputWidget.maxEvalsWidget.setValue(1000)
    widget.inputWidget.gaNbWidget.setValue(10)
    
    widget.show()
    sys.exit(app.exec_())
