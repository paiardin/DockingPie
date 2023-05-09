## Adopted from c++ code wwWidgets <wwwidgets@wysota.eu.org>

from PySide import QtGui, QtCore
import sys
import numpy as np
import math

# This class implements a simple color wheel widget for picking colors.

class HueSatRadialPicker(QtGui.QWidget):
    #colorPicked = QtCore.Signal(int, int, int, int)
    colorPicked = QtCore.Signal(object)
    valueChanged = QtCore.Signal(int)
    
    def __init__(self, parent=None, size=150):
        super(HueSatRadialPicker, self).__init__()
        rwidth = rheight = size
        self.setMinimumWidth(rwidth)
        self.setMinimumHeight(rheight)
        self.setMouseTracking(True)
        #policy= QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        #policy.setHeightForWidth(True);
        #self.setSizePolicy(policy)
        self.conical = QtGui.QConicalGradient()
        self.m_value = 255
        self.m_color = QtGui.QColor(255,255, 255, 255)
        self.m_pt = None #QtCore.QPoint(self.rect().center())
        self.px = None
        for i in range(360):
            self.conical.setColorAt(i*1.0/360, QtGui.QColor.fromHsv(i, 255, self.m_value))
        self.conical.setColorAt(1, QtGui.QColor.fromHsv(359, 255, self.m_value))
        #print "size:", self.size(), "rect", self.rect()
    
    def paintEvent(self, event):
        if not self.px:
            self.buildPixmap()
        if not self.m_pt:
            self.m_pt = QtCore.QPoint(self.rect().center())
        p = QtGui.QPainter()
        p.begin(self)
        p.drawPixmap(0, 0, self.px)
        self.drawCrosshair(p, self.m_pt)
        p.end()

    def buildPixmap(self):
        #build huecircle
        self.conical.setCenter(self.rect().center())
        huecircle = QtGui.QImage(self.size(), QtGui.QImage.Format_ARGB32)
        huecircle.fill(QtCore.Qt.transparent)
        phcircle = QtGui.QPainter(huecircle)
        #phcircle.begin(self)
        phcircle.setRenderHint(QtGui.QPainter.Antialiasing, True)
        phcircle.setBrush(QtGui.QBrush(self.conical))
        phcircle.setPen(self.palette().color(QtGui.QPalette.Shadow))
        phcircle.drawEllipse(self.rect().adjusted(1, 1, -1, -1))
        phcircle.end()
        self.px = px = QtGui.QPixmap.fromImage(huecircle)

        ## ## #alpha gradient
        rg = QtGui.QRadialGradient(self.rect().center(), self.width()/2-2, self.rect().center())

        for i in np.arange(0.0, 1.0, 0.1):
            rg.setColorAt(i, QtGui.QColor.fromHsv(0, 0, int(256*i)))
        rg.setColorAt(1, QtCore.Qt.white)

        # alpha channel
        ac = QtGui.QImage(self.size(), QtGui.QImage.Format_RGB32)
        ac.fill(QtCore.Qt.transparent)
        acpainter = QtGui.QPainter(ac)
        acpainter.setPen(QtCore.Qt.NoPen)
        acpainter.setBrush(QtGui.QBrush(rg))
        acpainter.drawEllipse(self.rect().adjusted(1, 1, -1, -1))
        acpainter.end()
        px.setAlphaChannel(QtGui.QPixmap.fromImage(ac))
        # destination image
        dst = QtGui.QImage(self.size(), QtGui.QImage.Format_ARGB32)
        dst.fill(QtCore.Qt.transparent)
        dstp = QtGui.QPainter(dst)
        dstp.setBrush(QtGui.QColor.fromHsv(0, 0, self.m_value))
        dstp.setRenderHint(QtGui.QPainter.Antialiasing, True)
        dstp.setPen(self.palette().color(QtGui.QPalette.Shadow))
        dstp.drawEllipse(self.rect().adjusted(1, 1, -1, -1))
        dstp.setCompositionMode(QtGui.QPainter.CompositionMode_SourceOver)
        dstp.drawPixmap(0, 0, px)
        dstp.end()
        self.px = QtGui.QPixmap.fromImage(dst)
        #p.drawPixmap(0, 0, px)
        #p.end()

    def drawCrosshair(self, painter, pt):
        if pt == None: return
        painter.save()
        painter.setPen(QtCore.Qt.black)
        painter.drawLine(pt-QtCore.QPoint(0, -3), pt-QtCore.QPoint(0, -1))
        painter.drawLine(pt-QtCore.QPoint(0, 1), pt-QtCore.QPoint(0, 3))
        painter.drawLine(pt-QtCore.QPoint(-3, 0), pt-QtCore.QPoint(-1, 0))
        painter.drawLine(pt-QtCore.QPoint(1, 0), pt-QtCore.QPoint(3, 0))
        painter.restore()

    def mousePressEvent(self, event):
        if (event.button()==QtCore.Qt.LeftButton):
            pt = QtCore.QPoint(event.pos())
            r = self.radius(pt)
            if (r>=self.width()/2): return
            self.m_pt = pt
            h = self.hue(pt)
            s = int(r*255.0/(self.width()/2))
            col = QtGui.QColor.fromHsv(h, s, self.m_value) #.getRgbF()
            #print "mousePressEvent", col
            self.m_color = col
            self.colorPicked.emit(col)
            self.update()
        else : QtGui.QWidget.mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if (event.buttons() and  QtCore.Qt.LeftButton):
            pt = QtCore.QPoint(event.pos())
            r = self.radius(pt)
            if (r>=self.width()/2):
                r = self.width()/2
                pt = self.m_pt
            self.m_pt = pt
            h = self.hue(pt)
            s = int(r*255.0/(self.width()/2))
            col = QtGui.QColor.fromHsv(h, s, self.m_value)#.getRgbF()
            #print "in mouseMoveEvent", col
            self.m_color = col
            self.colorPicked.emit(col)
            self.update()
        #else:
        #    QtGui.QWidget.mouseMoveEvent(self)

    def hue(self, pt):
        c = QtCore.QPoint(self.rect().center())
        x = c.x()-pt.x()
        y = pt.y()-c.y()
        d = QtCore.qAbs(x)
        if d == 0: d = 1e-6
        a = QtCore.qAbs(y)*1.0 / d   # tangent
        rd = math.atan(a)
        deg = rd * 180 /math.pi
        h = QtCore.qRound(deg)
        if x>=0 and y>=0:
            h=180+h
        elif x<0 and y>=0:
            h=360-h
        elif x<0 and y<0 :
            h=h
        else: h=180-h
        #print "hue:", h % 360
        return h % 360

    def radius(self, pt):
        c = QtCore.QPoint(self.rect().center())
        x = pt.x() - c.x()
        y = pt.y() - c.y()
        r = QtCore.qRound(math.sqrt(x*x + y*y))
        return r

    def setValue(self, val):
        if (val<0 or val>255 or val==self.m_value): return
        self.m_value = val
        self.valueChanged.emit(val)
        #self.conical.setStops(QtGui.QGradientStops());
        #for i in range(360):
        #    self.conical.setColorAt(i*1.0/360, QtGui.QColor.fromHsv(i, 255, self.m_value));
        #col = QtGui.QColor.fromHsv(359, 255, self.m_value)
        #self.conical.setColorAt(1, col)
        h = self.m_color.hue()
        s = self.m_color.saturation()
        self.m_color = QtGui.QColor.fromHsv(h, s, val)
        #self.buildPixmap()
        #self.update()

    def value(self):
        return self.m_value

    def setColor(self, col):
        if (col==self.m_color): return
        h, s, v = col.getHsv()
        if (v!=self.m_value): self.setValue(v)
        r = QtCore.qRound(s*1.0*(self.width()/2-2)/255.0)
        x = QtCore.qRound(r*math.cos(h*math.pi/180.0))
        y = QtCore.qRound(r*math.sin(h*math.pi/180.0))
        ctr = QtCore.QPoint(self.rect().center())
        self.m_pt = QtCore.QPoint(x+ctr.x(), -y+ctr.y())
        self.m_color = col
        self.colorPicked.emit(QtGui.QColor(col))
        self.update()
        
    def color(self):
        return self.m_color

    def heightForWidth(self, w):
        return w

    def resizeEvent(self, event):
        #print "in resizeEvent: width", self.width(), "height:", self.height()
        if (self.width()!=self.height()):
            s = self.qMin(self.width(), self.height())
            self.resize(s,s);
            return
        self.buildPixmap();
        self.update()

    def qMin(self, a, b):
        return a if a < b else b


class CBox(QtGui.QFrame):
    def __init__(self, parent=None, color=QtGui.QColor("black"), size=None):
        super(CBox, self).__init__(parent)
        self.setFrameShape(QtGui.QFrame.Box)
        self.setFrameShadow(QtGui.QFrame.Plain)
        #self.setLineWidth(0)
        #self.setMidLineWidth(0)
        #self.setContentsMargins(0, 0, 0, 0)
        self.setColor(color)
        if size is not None:
            assert len(size)==2
            self.setMaximumSize(size[0], size[1])
        
    def setColor(self, color):
        pal = self.palette()
        pal.setColor(QtGui.QPalette.Window, color)
        pal.setColor(QtGui.QPalette.WindowText, color)
        self.setAutoFillBackground(True)
        self.setPalette(pal)   

class ColorPicker(QtGui.QWidget):
    colorChanged = QtCore.Signal(object)
    
    def __init__(self, parent=None, colorWheelSize=150):
        super(ColorPicker, self).__init__(parent)
        grid = QtGui.QVBoxLayout()
        groupBox = QtGui.QGroupBox("")
        self.colorWheel = cw = HueSatRadialPicker(parent=groupBox, size=colorWheelSize)
        cw.colorPicked.connect(self.colorChanged_cb)
        vbox = QtGui.QVBoxLayout()
        lay1 = QtGui.QVBoxLayout()
        lay1.addWidget(cw)
        self.colorBox = cb = CBox(color=QtGui.QColor("white"))
        self.valueSlider = vs = QtGui.QSlider(groupBox)
        
        vs.setOrientation(QtCore.Qt.Horizontal)
        #vs.setOrientation(QtCore.Qt.Vertical)
        vs.setRange(0, 255)
        vs.setValue(255)
        #vs.setMaximumWidth(colorWheelSize)
        lay1.addWidget(vs)
        cb.setMinimumHeight(20)
        cb.setMaximumHeight(20)
        #cb.setMaximumWidth(colorWheelSize)
        lay1.addWidget(cb)
        lay2 = QtGui.QHBoxLayout()
        self.opac = op = QtGui.QDoubleSpinBox()
        op.setSingleStep(0.1)
        op.setValue(1.0)
        op.setRange(0.0, 1.0)
        op.valueChanged.connect(self.opacity_cb)
        vs.valueChanged.connect(self.valueChanged_cb)
        lab = QtGui.QLabel("alpha:")
        lay2.addWidget(lab)
        lay2.addWidget(op)
        vbox.addLayout(lay1)
        vbox.addLayout(lay2)        
        groupBox.setLayout(vbox)
        #vbox.addStretch(1)

        grid.addWidget(groupBox)
        self.setLayout(grid)

    def colorChanged_cb(self, color):
        self.colorBox.setColor(color)
        self.colorChanged.emit(self.addAlpha(color))
        
    def valueChanged_cb(self, value):
        self.colorWheel.setValue(value)
        color = self.colorWheel.color()
        #print "valueChanged_cb", value, color
        self.colorBox.setColor(color)
        self.colorChanged.emit(self.addAlpha(color))
        
    def getRgba(self):
        # get colorwheel QColor obj (alpha is 255) 
        return self.addAlpha(self.colorWheel.color())

    def addAlpha(self, color, alpha=None):
        # modify alpha value of the input with the value from the opacity spin box
        r, g, b, a = color.getRgb()
        if alpha is not None:
            a = alpha
        else:
            a = round(255 * float(self.opac.value()))
        return QtGui.QColor(r,g,b,int(a))
        
    def opacity_cb(self, val):
        alpha = int(round(255 * val))
        color = self.colorWheel.color()         
        self.colorChanged.emit(self.addAlpha(color, alpha))

    def getOpacity(self):
        return float(self.opac.value())
    
        
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    cl = ColorPicker(colorWheelSize=150)
    cl.show()
    sys.exit(app.exec_())

