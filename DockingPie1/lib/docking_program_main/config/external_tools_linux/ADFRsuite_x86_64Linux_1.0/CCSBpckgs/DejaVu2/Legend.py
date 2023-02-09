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

########################################################################
#
# Date: Febuary 2006 Authors: Guillaume Vareille, Michel Sanner
#
#    vareille@scripps.edu
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Guillaume Vareille, Michel Sanner and TSRI
#
# Revision:
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/Legend.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: Legend.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

from copy import deepcopy

from opengltk.OpenGL import GL
from mglutil.util.colorUtil import ToHSV
from pyglf import glf

def drawSelfOrientedLegend(
                    fullWidth,
                    fullHeight,
                    ramp,
                    mini,
                    maxi,
                    name = '',
                    unit = '',
                    labelValues = [],
                    roomLeftToLegend = 0,
                    roomBelowLegend = 50,
                    legendShortSide = 10,
                    legendLongSide = 150,
                    significantDigits = 3,
                    backgroundColor = (0,0,0,.8),
                    interpolate = True,
                    frame=True,
                    selected=False,
                    numOfLabels=None,
                    resizeSpotRadius=5,
                    fontScale=8,
                    glfFontID=0,
                    tile=None,
                    ):   

    lRoomLeftToViewportCenter = fullWidth / 2
    lRoombelowViewportCenter = fullHeight / 2

    twiceFontHeight = 2*fontScale
    
    if roomLeftToLegend < lRoomLeftToViewportCenter:
        if roomBelowLegend < lRoombelowViewportCenter: # 1 and 8
            if roomBelowLegend > roomLeftToLegend: #1
                #print "1"
                lVerticalLegend=True
                lLeftOrBelowLabels=False
                if roomLeftToLegend < 0:
                    roomLeftToLegend = 0
                if roomBelowLegend < 0:
                    roomBelowLegend = 0   
            else: # 8
                #print "8"
                lVerticalLegend=False
                lLeftOrBelowLabels=False
                if roomLeftToLegend < 0:
                    roomLeftToLegend = 0
                if roomBelowLegend < 0:
                    roomBelowLegend = 0   
        else: # 2 and 3
            if fullHeight-roomBelowLegend > roomLeftToLegend + (fullHeight/8): #2
                #print "2"
                lVerticalLegend=True
                lLeftOrBelowLabels=False
                if roomLeftToLegend < 0:
                    roomLeftToLegend = 0
                if roomBelowLegend > fullHeight:
                    roomBelowLegend = fullHeight  
            else: #3
                #print "3"
                lVerticalLegend=False
                lLeftOrBelowLabels=True
                if roomLeftToLegend < 0:
                    roomLeftToLegend = 0
                if roomBelowLegend > fullHeight:
                    roomBelowLegend = fullHeight  
    else:
        if roomBelowLegend < lRoombelowViewportCenter: # 6 and 7
            if roomBelowLegend + (fullWidth/8) < fullWidth-roomLeftToLegend: #7
                #print "7"
                lVerticalLegend=False
                lLeftOrBelowLabels=False
                if roomLeftToLegend > fullWidth :
                    roomLeftToLegend = fullWidth
                if roomBelowLegend < 0:
                    roomBelowLegend = 0   
            else: # 6
                #print "6"
                lVerticalLegend=True
                lLeftOrBelowLabels=True
                if roomLeftToLegend > fullWidth :
                    roomLeftToLegend = fullWidth
                if roomBelowLegend < twiceFontHeight:
                    roomBelowLegend = twiceFontHeight   
        else: # 4 and 5
            if fullHeight-roomBelowLegend < fullWidth-roomLeftToLegend: #2
                #print "4"
                lVerticalLegend=False
                lLeftOrBelowLabels=True
                if roomLeftToLegend > fullWidth :
                    roomLeftToLegend = fullWidth
                if roomBelowLegend > fullHeight:
                    roomBelowLegend = fullHeight  
            else: #5
                #print "5"
                lVerticalLegend=True
                lLeftOrBelowLabels=True
                if roomLeftToLegend > fullWidth :
                    roomLeftToLegend = fullWidth
                if roomBelowLegend > fullHeight:
                    roomBelowLegend = fullHeight 

    if lVerticalLegend is True:
        if legendLongSide > fullHeight:
            legendLongSide = fullHeight
        if legendShortSide > fullWidth:
            legendShortSide = fullWidth
    else:
        if legendLongSide > fullWidth:
            legendLongSide = fullWidth
        if legendShortSide > fullHeight:
            legendShortSide = fullHeight

    bgHSV = ToHSV(backgroundColor[:3]) 
    if bgHSV[2] < .5:
        labelColor = (1, 1, 1)
    else:
        labelColor = (0, 0, 0)

    return drawLegendLabelName(
        fullWidth=fullWidth,
        fullHeight=fullHeight,
        ramp=ramp,
        mini=mini,
        maxi=maxi,
        name=name,
        unit=unit,
        labelValues=labelValues,
        verticalLegend=lVerticalLegend,
        leftOrBelowLabels=lLeftOrBelowLabels,
        roomLeftToLegend=roomLeftToLegend,
        roomBelowLegend=roomBelowLegend,
        legendShortSide=legendShortSide,
        legendLongSide=legendLongSide,
        significantDigits=significantDigits,
        backgroundColor=backgroundColor,
        labelColor=labelColor,
        interpolate=interpolate,
        frame=frame,
        selected=selected,
        numOfLabels=numOfLabels,
        resizeSpotRadius=resizeSpotRadius,
        fontScale=fontScale,
        glfFontID=glfFontID,
        tile=tile,
        )


def drawLegendLabelName(
                    fullWidth,
                    fullHeight,
                    ramp,
                    mini,
                    maxi,
                    name='',
                    unit='',
                    labelValues=[],
                    verticalLegend=True,
                    leftOrBelowLabels=False,
                    roomLeftToLegend=0,
                    roomBelowLegend=50,
                    legendShortSide=10,
                    legendLongSide=150,
                    significantDigits=3,
                    backgroundColor=(0,0,0,.8),
                    labelColor=(1,1,1),
                    interpolate=True,
                    frame=True,
                    selected=False,
                    numOfLabels=None,
                    resizeSpotRadius=5,
                    fontScale=8,
                    glfFontID=0,
                    tile=None,
                    ):
    
    if ramp is None or len(ramp) == 0:
        return

    if tile is not None and selected is True:
        selected = False

    # some glf font initialisation 
    glf.glfSetCurrentFont(glfFontID)
    glf.glfStringCentering(GL.GL_FALSE)
    glf.glfStringDirection(glf.GLF_LEFT)

    # calculate name and unit size
    lNameMinAndMax = glf.glfGetStringBounds(name)
    lNameMinAndMax = ( lNameMinAndMax[0]*fontScale,
                       lNameMinAndMax[1]*fontScale,
                       lNameMinAndMax[2]*fontScale,
                       lNameMinAndMax[3]*fontScale )
    lNameWidth = lNameMinAndMax[2]-lNameMinAndMax[0]
    lNameHeight = lNameMinAndMax[3]-lNameMinAndMax[1]
    if lNameWidth == 0:
        lNameWidth = fontScale
    if unit is not None and (len(unit) > 0):
        lUnitMinAndMax = glf.glfGetStringBounds(unit)
        lUnitMinAndMax = ( lUnitMinAndMax[0]*fontScale,
                           lUnitMinAndMax[1]*fontScale,
                           lUnitMinAndMax[2]*fontScale,
                           lUnitMinAndMax[3]*fontScale )
        lUnitWidth = lUnitMinAndMax[2]-lUnitMinAndMax[0]
        lUnitHeight = lUnitMinAndMax[3]-lUnitMinAndMax[1]
    else:
        lUnitWidth = fontScale
        lUnitHeight = 0
    lMaxNameUnitHeight = max(lNameHeight, lUnitHeight)
    lMaxNameUnitWidth = max(lNameWidth, lUnitWidth)

    # deducted values
    fontScaleHalf = fontScale/2
    if verticalLegend is True:
        lRoomToLegendCloseLongSide = roomLeftToLegend
        lRoomToLegendCloseShortSide = roomBelowLegend
        if legendShortSide < 1:
            legendShortSide = 1
        if legendLongSide < 1:
            legendLongSide = 1
    else:
        lRoomToLegendCloseLongSide = roomBelowLegend 
        lRoomToLegendCloseShortSide = roomLeftToLegend 
        if legendShortSide < lMaxNameUnitHeight + fontScaleHalf:
            legendShortSide = lMaxNameUnitHeight + fontScaleHalf
        if legendLongSide < 1:
            legendLongSide = 1

    lRoomToLegendFarLongSide = lRoomToLegendCloseLongSide + legendShortSide 
    lRoomToLegendFarShortSide = lRoomToLegendCloseShortSide + legendLongSide 

    # prepare the legend labels
    if len(labelValues) > 0:
        lOnScreenLabelValues = labelValues
    else:        
        if numOfLabels == None or numOfLabels < 0:
            if verticalLegend is True:
                numOfLabels = 6
            else:
                numOfLabels = 4
        
        if ( numOfLabels == 0 ) or maxi is None or mini is None:
            lOnScreenLabelValues = []
        elif numOfLabels == 1:
            lOnScreenLabelValues = [(mini+maxi)*0.5]
        else: 
            delta = (maxi - mini) / float(numOfLabels-1)
            lOnScreenLabelValues = []
            for i in range(numOfLabels):
                lOnScreenLabelValues.append(mini+i*delta)

    lMaxLabelDigits = 0
    for v in lOnScreenLabelValues:
        lLabel = "%.*g" % (significantDigits,v)
        if lMaxLabelDigits < len(lLabel):
            lMaxLabelDigits = len(lLabel)

    # calculate labels size
    lLabelsMinAndMax = []
    lLabelsWidth = []
    lLabelsHeight = []
    lMaxLabelsWidth = 0
    lMaxLabelsHeight = 0
    for v in lOnScreenLabelValues:
        lLabel = "%.*g" % (significantDigits,v)
        lLabelsMinAndMax.append( glf.glfGetStringBounds(lLabel) )
        lLabelsMinAndMax[-1] = ( lLabelsMinAndMax[-1][0]*fontScale,
                           lLabelsMinAndMax[-1][1]*fontScale,
                           lLabelsMinAndMax[-1][2]*fontScale,
                           lLabelsMinAndMax[-1][3]*fontScale )
        lLabelsWidth.append( lLabelsMinAndMax[-1][2]-lLabelsMinAndMax[-1][0] )
        lLabelsHeight.append( lLabelsMinAndMax[-1][3]-lLabelsMinAndMax[-1][1] )
        if lLabelsWidth[-1] > lMaxLabelsWidth:
            lMaxLabelsWidth = lLabelsWidth[-1]
        if lLabelsHeight[-1] > lMaxLabelsHeight:
            lMaxLabelsHeight = lLabelsHeight[-1]

    # calculate frame size
    lMaxWidth = max(lMaxLabelsWidth + legendShortSide, lMaxNameUnitWidth)
    lMaxNameUnitHeightHalf = lMaxNameUnitHeight/2
    legendShortSideHalf = legendShortSide/2
    if verticalLegend is True:
        if leftOrBelowLabels is False:
            lPt1 = (lRoomToLegendCloseLongSide,
                    -fontScale+lRoomToLegendCloseShortSide-lMaxLabelsHeight-lNameHeight,0)
            lPt2 = (fontScale+lRoomToLegendCloseLongSide+lMaxWidth, 
                    -fontScale+lRoomToLegendCloseShortSide-lMaxLabelsHeight-lNameHeight,0)
            lPt3 = (fontScale+lRoomToLegendCloseLongSide+lMaxWidth, 
                    fontScale+lRoomToLegendFarShortSide+lMaxLabelsHeight+lUnitHeight,0)
            lPt4 = (lRoomToLegendCloseLongSide, 
                    fontScale+lRoomToLegendFarShortSide+lMaxLabelsHeight+lUnitHeight,0)
        else:
            lPt1 = (-fontScale+lRoomToLegendFarLongSide-lMaxWidth,
                    -fontScale+lRoomToLegendCloseShortSide-lMaxLabelsHeight-lNameHeight,0)
            lPt2 = (lRoomToLegendCloseLongSide+legendShortSide, 
                    -fontScale+lRoomToLegendCloseShortSide-lMaxLabelsHeight-lNameHeight,0)
            lPt3 = (lRoomToLegendCloseLongSide+legendShortSide, 
                    fontScale+lRoomToLegendFarShortSide+lMaxLabelsHeight+lUnitHeight,0)
            lPt4 = (-fontScale+lRoomToLegendFarLongSide-lMaxWidth, 
                    fontScale+lRoomToLegendFarShortSide+lMaxLabelsHeight+lUnitHeight,0)
    else:
        if leftOrBelowLabels is False:
            lPt1 = (-fontScale+lRoomToLegendCloseShortSide-lNameWidth,
                    lRoomToLegendCloseLongSide,0)
            lPt2 = (fontScale+lRoomToLegendFarShortSide+lUnitWidth,
                    lRoomToLegendCloseLongSide,0)
            lPt3 = (fontScale+lRoomToLegendFarShortSide+lUnitWidth, 
                    fontScale+lRoomToLegendCloseLongSide \
                    + legendShortSide \
                    + lMaxLabelsHeight,
                    0)
            lPt4 = (-fontScale+lRoomToLegendCloseShortSide-lNameWidth, 
                    fontScale+lRoomToLegendCloseLongSide \
                    + legendShortSide \
                    + lMaxLabelsHeight,
                    0)
        else:
            lPt1 = (-fontScale+lRoomToLegendCloseShortSide-lNameWidth,
                    -fontScale+lRoomToLegendCloseLongSide-lMaxLabelsHeight,
                    0)
            lPt2 = (fontScale+lRoomToLegendFarShortSide+lUnitWidth,
                    -fontScale+lRoomToLegendCloseLongSide-lMaxLabelsHeight,
                    0)
            lPt3 = (fontScale+lRoomToLegendFarShortSide+lUnitWidth,
                    lRoomToLegendFarLongSide,0)
            lPt4 = (-fontScale+lRoomToLegendCloseShortSide-lNameWidth,
                    lRoomToLegendFarLongSide,0)

    if selected is True:
        labelColor2 = (.5, .5, .5)
    else:
        labelColor2 = labelColor

    GL.glMatrixMode(GL.GL_PROJECTION)
    GL.glPushMatrix()
    GL.glLoadIdentity()
    if tile is None:
        GL.glOrtho(0, float(fullWidth), 0, float(fullHeight), -1, 1)
    else:
        GL.glOrtho(float(tile[0]), float(tile[1]), float(tile[2]), float(tile[3]), -1, 1)
    GL.glMatrixMode(GL.GL_MODELVIEW)
    GL.glPushMatrix()
    GL.glLoadIdentity()
    GL.glDisable( GL.GL_LIGHTING )

    GL.glPolygonMode(GL.GL_FRONT, GL.GL_FILL )
    
    if len(backgroundColor) == 4:
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)

    #because of an unexplained bug on michel's laptop,
    #we don't draw the background when there is no frame 
    #(so michel has a way to manage the legend)
    if frame is True:
        #draw transparent background
        GL.glDepthMask(GL.GL_FALSE) 
        if len(backgroundColor) == 3:
            GL.glColor3fv(backgroundColor)
        else:
            GL.glColor4fv(backgroundColor)           
        GL.glBegin(GL.GL_QUADS)
        GL.glVertex3fv(lPt1)
        GL.glVertex3fv(lPt2)
        GL.glVertex3fv(lPt3)
        GL.glVertex3fv(lPt4)
        GL.glEnd()
        GL.glDepthMask(GL.GL_TRUE) 

    #draw frame
    if frame is True:
        GL.glPolygonMode(GL.GL_FRONT, GL.GL_LINE )
        GL.glLineWidth(1)
        GL.glColor3fv(labelColor2)    
        GL.glBegin(GL.GL_QUADS)
        GL.glVertex3fv(lPt1)
        GL.glVertex3fv(lPt2)
        GL.glVertex3fv(lPt3)
        GL.glVertex3fv(lPt4)
        GL.glEnd()
        GL.glPolygonMode(GL.GL_FRONT, GL.GL_FILL)
    
    if mini is not None and maxi is not None and maxi > mini:
        lUnitStep = legendLongSide/float(maxi-mini)
    else:
        lUnitStep = legendLongSide
    
    GL.glDisable( GL.GL_LIGHTING )

    if verticalLegend is True:
        if leftOrBelowLabels is True:
            lRoomLeftToLabel = -fontScaleHalf + lRoomToLegendCloseLongSide
            lRoomLeftToName = -fontScaleHalf + lRoomToLegendFarLongSide - lNameWidth
            lRoomLeftToUnit = -fontScaleHalf + lRoomToLegendFarLongSide - lUnitWidth
            lRoomBelowName = -fontScaleHalf + roomBelowLegend - lMaxLabelsHeight - lNameHeight
            lRoomBelowUnit = fontScaleHalf + roomBelowLegend + legendLongSide + lMaxLabelsHeight
        else:
            lRoomLeftToLabel = fontScaleHalf + lRoomToLegendFarLongSide
            lRoomLeftToName = fontScaleHalf + lRoomToLegendCloseLongSide    
            lRoomLeftToUnit = lRoomLeftToName
            lRoomBelowName = -fontScaleHalf + roomBelowLegend - lMaxLabelsHeight - lNameHeight
            lRoomBelowUnit = fontScaleHalf + roomBelowLegend + legendLongSide + lMaxLabelsHeight
    else:
        if leftOrBelowLabels is True:
            lRoomBelowLabel = -fontScaleHalf + lRoomToLegendFarLongSide \
                              - legendShortSide - lMaxLabelsHeight
            lRoomBelowName = lRoomToLegendCloseLongSide + legendShortSideHalf \
                             - lMaxNameUnitHeightHalf
        else:
            lRoomBelowLabel = fontScaleHalf + lRoomToLegendCloseLongSide \
                              + legendShortSide
            lRoomBelowName = lRoomToLegendCloseLongSide + legendShortSideHalf \
                             - lMaxNameUnitHeightHalf
        lRoomBelowUnit = lRoomBelowName
        lRoomLeftToName = -fontScaleHalf + lRoomToLegendCloseShortSide - lNameWidth
        lRoomLeftToUnit = fontScaleHalf + lRoomToLegendFarShortSide

    # set the color of the text
    GL.glColor3fv(labelColor2)

    # print the legend name
    GL.glPushMatrix()
    GL.glTranslatef(
            float(lRoomLeftToName+fontScale),
            float(lRoomBelowName-lNameMinAndMax[1]),
            0)
    GL.glScalef(float(fontScale), float(fontScale), 0);
    glf.glfDrawSolidString(name)
    GL.glPopMatrix()
    if unit is not None and (len(unit) > 0):
        GL.glPushMatrix()
        GL.glTranslatef(
            float(lRoomLeftToUnit+fontScale),
            float(lRoomBelowUnit-lUnitMinAndMax[1]),
            0)
        GL.glScalef(float(fontScale), float(fontScale), 1);
        glf.glfDrawSolidString(unit)
        GL.glPopMatrix()
       
    if mini is not None and maxi is not None:
        i = 0
        for v in lOnScreenLabelValues:
            #calculate label position
            lLabel = "%.*g" % (significantDigits,v)
            lStep = (v - mini) * lUnitStep
            GL.glPushMatrix()
            if verticalLegend:
                if leftOrBelowLabels:
                    GL.glTranslatef(
                        float(lRoomLeftToLabel+fontScale-lLabelsWidth[i]),
                        float(roomBelowLegend-(lLabelsHeight[i]/2+lLabelsMinAndMax[i][1])+lStep),
                        0)
                else:
                    GL.glTranslatef(
                        float(lRoomLeftToLabel+fontScale),
                        float(roomBelowLegend-(lLabelsHeight[i]/2+lLabelsMinAndMax[i][1])+lStep),
                        0)
            else:
                GL.glTranslatef(
                    float(roomLeftToLegend-(lLabelsWidth[i]/2)+fontScale+lStep),
                    float(lRoomBelowLabel-lLabelsMinAndMax[i][1]),
                    0)
            GL.glScalef(float(fontScale), float(fontScale), 1);
            glf.glfDrawSolidString("%s" % lLabel)
            GL.glPopMatrix()
            i += 1

    if len(backgroundColor) == 4:
        GL.glDisable(GL.GL_BLEND)

#    GL.glEnable(GL.GL_LIGHTING)
    GL.glPopMatrix()
    GL.glMatrixMode(GL.GL_PROJECTION)
    GL.glPopMatrix()
    GL.glMatrixMode(GL.GL_MODELVIEW)

    drawLegendOnly(
        fullWidth=fullWidth,
        fullHeight=fullHeight,
        ramp=ramp,
        verticalLegend=verticalLegend,
        roomLeftToLegend=roomLeftToLegend,
        roomBelowLegend=roomBelowLegend,
        legendShortSide=legendShortSide,
        legendLongSide=legendLongSide,
        interpolate=interpolate,
        selected=selected,
        tile=tile,
        )

    if selected is True:
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glPushMatrix()
        GL.glLoadIdentity()
        GL.glOrtho(0, float(fullWidth), 0, float(fullHeight), -1, 1)
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glPushMatrix()
        GL.glLoadIdentity()
        GL.glDisable( GL.GL_LIGHTING )
        GL.glPolygonMode(GL.GL_FRONT, GL.GL_FILL )
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glDepthMask(GL.GL_FALSE)
        GL.glDisable(GL.GL_LIGHTING)

        resizeSpot = []
        if verticalLegend:
            resizeSpot= [ roomLeftToLegend + legendShortSide,
                          roomBelowLegend + legendLongSide ]
        else:
            resizeSpot = [ roomLeftToLegend + legendLongSide,
                           roomBelowLegend + legendShortSide ]

        GL.glColor3fv(labelColor)
        GL.glBegin(GL.GL_QUADS)
        GL.glVertex2f(float(resizeSpot[0]+resizeSpotRadius),float(resizeSpot[1]-resizeSpotRadius))
        GL.glVertex2f(float(resizeSpot[0]+resizeSpotRadius),float(resizeSpot[1]+resizeSpotRadius))
        GL.glVertex2f(float(resizeSpot[0]-resizeSpotRadius),float(resizeSpot[1]+resizeSpotRadius))
        GL.glVertex2f(float(resizeSpot[0]-resizeSpotRadius),float(resizeSpot[1]-resizeSpotRadius))
        GL.glEnd()

        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glDepthMask(GL.GL_TRUE) 
#        GL.glEnable(GL.GL_LIGHTING)
        GL.glPopMatrix()
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glPopMatrix()
        GL.glMatrixMode(GL.GL_MODELVIEW)
    else:
        resizeSpot = None

    return [(lPt1[0], lPt1[1]), 
            (lPt2[0], lPt2[1]), 
            (lPt3[0], lPt3[1]), 
            (lPt4[0], lPt4[1]) ] , resizeSpot , verticalLegend


def drawLegendOnly(
                    fullWidth,
                    fullHeight,
                    ramp,
                    verticalLegend=True,
                    roomLeftToLegend=0,
                    roomBelowLegend=50,
                    legendShortSide=10,
                    legendLongSide=150,
                    interpolate=True,
                    selected=False,
                    tile=None,
                    ):

    if ramp is None or len(ramp) == 0:
        return

    if tile is not None and selected is True:
        selected = False

    # deducted values
    if verticalLegend is True:
        lRoomToLegendCloseLongSide = roomLeftToLegend
        lRoomToLegendCloseShortSide = roomBelowLegend
        if legendShortSide < 1:
            legendShortSide = 1
        if legendLongSide < 1:
            legendLongSide = 1
    else:
        lRoomToLegendCloseLongSide = roomBelowLegend 
        lRoomToLegendCloseShortSide = roomLeftToLegend 
        if legendShortSide < 1:
            legendShortSide = 1
        if legendLongSide < 1:
            legendLongSide = 1

    lRoomToLegendFarLongSide = lRoomToLegendCloseLongSide + legendShortSide 

    GL.glMatrixMode(GL.GL_PROJECTION)
    GL.glPushMatrix()
    GL.glLoadIdentity()
    if tile is None:
        GL.glOrtho(0, float(fullWidth), 0, float(fullHeight), -1, 1)
    else:
        GL.glOrtho(float(tile[0]), float(tile[1]), float(tile[2]), float(tile[3]), -1, 1)
    GL.glMatrixMode(GL.GL_MODELVIEW)
    GL.glPushMatrix()
    GL.glLoadIdentity()
    GL.glDisable( GL.GL_LIGHTING )

    GL.glPolygonMode(GL.GL_FRONT, GL.GL_FILL)
    GL.glDisable(GL.GL_DEPTH_TEST)
    GL.glDepthMask(GL.GL_FALSE)
    
    if len(ramp[0]) == 4: # there are alpha values, draw checkered bg
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        #draw a bunch of quads as background
        lCheckerSquareSize = legendShortSide * .5
        if lCheckerSquareSize != 0:
            nbquads = int(legendLongSide / lCheckerSquareSize)
        else:
            nbquads = 1
        c1 = ( 0.1, 0.1, 0.1 )
        c2 = ( 0.3, 0.3, 0.3 )
        c = c1
        x2 = None
        for i in range(nbquads+1):
            GL.glColor3fv(c)
            if i==nbquads and nbquads != 0:
                x1=x2
                x2=legendLongSide
            else:
                x1 = i*lCheckerSquareSize
                x2 = min (x1+lCheckerSquareSize, legendLongSide)

            GL.glBegin(GL.GL_QUADS)
            if verticalLegend is True:
                GL.glVertex2f(float(lRoomToLegendCloseLongSide),
                              float(lRoomToLegendCloseShortSide + x1))
                GL.glVertex2f(float(lRoomToLegendCloseLongSide+lCheckerSquareSize), 
                              float(lRoomToLegendCloseShortSide + x1))
                GL.glVertex2f(float(lRoomToLegendCloseLongSide+lCheckerSquareSize), 
                              float(lRoomToLegendCloseShortSide + x2))
                GL.glVertex2f(float(lRoomToLegendCloseLongSide), 
                              float(lRoomToLegendCloseShortSide + x2))
            else:
                GL.glVertex2f(float(lRoomToLegendCloseShortSide + x2),
                              float(lRoomToLegendCloseLongSide))
                GL.glVertex2f(float(lRoomToLegendCloseShortSide + x2),
                              float(lRoomToLegendCloseLongSide+lCheckerSquareSize))
                GL.glVertex2f(float(lRoomToLegendCloseShortSide + x1),
                              float(lRoomToLegendCloseLongSide+lCheckerSquareSize))
                GL.glVertex2f(float(lRoomToLegendCloseShortSide + x1),
                              float(lRoomToLegendCloseLongSide))
            GL.glEnd()
    
            if c==c1:
                c=c2
            else:
                c=c1
    
            GL.glColor3fv(c)
            GL.glBegin(GL.GL_QUADS)
            if verticalLegend is True:
                GL.glVertex2f(float(lRoomToLegendCloseLongSide+lCheckerSquareSize), 
                              float(lRoomToLegendCloseShortSide + x1))
                GL.glVertex2f(float(lRoomToLegendFarLongSide), 
                              float(lRoomToLegendCloseShortSide + x1))
                GL.glVertex2f(float(lRoomToLegendFarLongSide), 
                              float(lRoomToLegendCloseShortSide + x2))
                GL.glVertex2f(float(lRoomToLegendCloseLongSide+lCheckerSquareSize), 
                              float(lRoomToLegendCloseShortSide + x2))
            else:
                GL.glVertex2f(float(lRoomToLegendCloseShortSide + x2),
                              float(lRoomToLegendCloseLongSide+lCheckerSquareSize))
                GL.glVertex2f(float(lRoomToLegendCloseShortSide + x2),
                              float(lRoomToLegendFarLongSide))
                GL.glVertex2f(float(lRoomToLegendCloseShortSide + x1),
                              float(lRoomToLegendFarLongSide))
                GL.glVertex2f(float(lRoomToLegendCloseShortSide + x1),
                              float(lRoomToLegendCloseLongSide+lCheckerSquareSize))
            GL.glEnd()
        
    #interp = False
    if interpolate and (len(ramp) > 1): # we interpolate colors
        # draw a quad strip for the color map
        lDelta = legendLongSide/float(len(ramp)-1)
        GL.glBegin(GL.GL_QUAD_STRIP)
        for i in range(len(ramp)):
            if selected is True:
                c = deepcopy(ramp[i])
                for j in range(3):
                     c[j] = .35 + c[j]*.3
            else:
                c = ramp[i]
            if len(c)==3:
                GL.glColor3fv(c)
            elif len(c)==4:
                GL.glColor4fv(c)
            lStep = i*lDelta
            if verticalLegend is True:
                GL.glVertex2f(float(lRoomToLegendCloseLongSide), 
                              float(lRoomToLegendCloseShortSide + lStep))
                GL.glVertex2f(float(lRoomToLegendFarLongSide), 
                              float(lRoomToLegendCloseShortSide + lStep))
            else:
                GL.glVertex2f(float(lRoomToLegendCloseShortSide + lStep),
                              float(lRoomToLegendFarLongSide))
                GL.glVertex2f(float(lRoomToLegendCloseShortSide + lStep),
                              float(lRoomToLegendCloseLongSide))
        GL.glEnd()
    else: 
        # we draw a quad for each color
        lDelta = legendLongSide/float(len(ramp))
        for i in range(len(ramp)):
            if selected is True:
                c = deepcopy(ramp[i])
                for j in range(3):
                     c[j] = .35 + c[j]*.3
            else:
                c = ramp[i]
            if len(c)==3:
                GL.glColor3fv(c)
            elif len(c)==4:
                GL.glColor4fv(c)
            x1 = i*lDelta
            x2 = x1+lDelta
            GL.glBegin(GL.GL_QUADS)
            if verticalLegend is True:
                GL.glVertex2f(float(lRoomToLegendCloseLongSide),
                              float(lRoomToLegendCloseShortSide + x1))
                GL.glVertex2f(float(lRoomToLegendFarLongSide), 
                              float(lRoomToLegendCloseShortSide + x1))
                GL.glVertex2f(float(lRoomToLegendFarLongSide), 
                              float(lRoomToLegendCloseShortSide + x2))
                GL.glVertex2f(float(lRoomToLegendCloseLongSide), 
                              float(lRoomToLegendCloseShortSide + x2))
            else:
                GL.glVertex2f(float(lRoomToLegendCloseShortSide + x2),
                              float(lRoomToLegendCloseLongSide))
                GL.glVertex2f(float(lRoomToLegendCloseShortSide + x2),
                              float(lRoomToLegendFarLongSide))
                GL.glVertex2f(float(lRoomToLegendCloseShortSide + x1),
                              float(lRoomToLegendFarLongSide))
                GL.glVertex2f(float(lRoomToLegendCloseShortSide + x1),
                              float(lRoomToLegendCloseLongSide))
            GL.glEnd()
    
    if len(ramp[0])==4:
        GL.glDisable(GL.GL_BLEND)

    GL.glEnable(GL.GL_DEPTH_TEST)
    GL.glDepthMask(GL.GL_TRUE) 

#    GL.glEnable(GL.GL_LIGHTING)
    GL.glPopMatrix()
    GL.glMatrixMode(GL.GL_PROJECTION)
    GL.glPopMatrix()
    GL.glMatrixMode(GL.GL_MODELVIEW)
