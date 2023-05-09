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

"""
Python reimplementation of Brian Paul's TR (tile Rendering) library 
/* $Id: tileRenderer.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
/*
 * Tiled Rendering library
 * Version 1.1
 * Copyright (C) Brian Paul
 */
 *
 * This library allows one to render arbitrarily large images with OpenGL.
 * The basic idea is to break the image into tiles which are rendered one
 * at a time.  The tiles are assembled together to form the final, large
 * image.  Tiles and images can be of any size.
 *
 * Basic usage:
 *
 * 1. Allocate a tile rendering context:
 *       TRcontext t = trNew();
 *
 * 2. Specify the final image buffer and tile size:
 *       GLubyte image[W][H][4]
 *       trImageSize(t, W, H);
 *       trImageBuffer(t, GL_RGBA, GL_UNSIGNED_BYTE, (GLubyte *) image);
 *
 * 3. Setup your projection:
 *       trFrustum(t, left, right, bottom top, near, far);
 *    or
 *       trOrtho(t, left, right, bottom top, near, far);
 *    or
 *       trPerspective(t, fovy, aspect, near, far);
 *
 * 4. Render the tiles:
 *       do {
 *           trBeginTile(t);
 *           DrawMyScene();
 *       } while (trEndTile(t));
 *
 *    You provide the DrawMyScene() function which calls glClear() and
 *    draws all your stuff.
 *
 * 5. The image array is now complete.  Display it, write it to a file, etc.
 *
 * 6. Delete the tile rendering context when finished:
 *       trDelete(t);
"""

TR_VERSION = "1.1"
TR_MAJOR_VERSION = 1
TR_MINOR_VERSION = 1

TR_TILE_WIDTH = 100
TR_TILE_HEIGHT = 101
TR_TILE_BORDER = 102
TR_IMAGE_WIDTH = 103
TR_IMAGE_HEIGHT = 104
TR_ROWS = 105
TR_COLUMNS = 106
TR_CURRENT_ROW = 107
TR_CURRENT_COLUMN = 108
TR_CURRENT_TILE_WIDTH = 109
TR_CURRENT_TILE_HEIGHT = 110
TR_ROW_ORDER = 111
TR_TOP_TO_BOTTOM = 112
TR_BOTTOM_TO_TOP = 113

from opengltk.OpenGL import GL
from opengltk.OpenGL.GLU import gluProject
from opengltk.extent import _gllib as gllib

class TRctx:
   """Tile Rendring context object
"""

   DEFAULT_TILE_WIDTH = 256
   DEFAULT_TILE_HEIGHT = 256
   DEFAULT_TILE_BORDER = 0

   def __init__(self, zmin, zmax):

      self.zmin = zmin
      self.zmax = zmax

      # Final image parameters
      self.ImageWidth = 1000
      self.ImageHeight = 1000
      self.ImageFormat = None
      self.ImageType = None
      self.ImageBuffer = None

      # Tile parameters
      self.TileWidth = self.DEFAULT_TILE_WIDTH
      self.TileHeight = self.DEFAULT_TILE_HEIGHT
      self.TileWidthNB = self.DEFAULT_TILE_BORDER
      self.TileHeightNB = 0
      self.TileBorder = 0
      self.TileFormat = None
      self.TileType = None
      self.TileBuffer = None

      # tile projection params used to get global antialiasing to work
      self.tileleft = None
      self.tileright = None
      self.tilebottom = None
      self.tiletop = None

      # Projection parameters
      self.Perspective = True
      self.Left = -1.
      self.Right = 1.
      self.Bottom = -1.
      self.Top = 1.
      self.Near = 0.1
      self.Far = 50

      # Misc
      self.RowOrder = TR_BOTTOM_TO_TOP
      self.Rows = 2
      self.Columns = 2
      self.CurrentTile = -1
      self.CurrentTileWidth = 100
      self.CurrentTileHeight = 100
      self.CurrentRow = 0
      self.CurrentColumn = 0
      self.backBuffer = True   # render tiles using back buffer
      
      self.ViewportSave = (0,0,100,100)


   #Misc setup including computing number of tiles (rows and columns).

   def setup(self):
      self.Columns = (self.ImageWidth + self.TileWidthNB - 1) / \
                     self.TileWidthNB
      self.Rows = (self.ImageHeight + self.TileHeightNB - 1) / \
                  self.TileHeightNB;
      self.CurrentTile = 0

      assert(self.Columns >= 0)
      assert(self.Rows >= 0)


   def tileSize(self, width, height, border):

      assert(border >= 0)
      assert(width >= 1)
      assert(height >= 1)
      assert(width >= 2*border)
      assert(height >= 2*border)

      self.TileBorder = border
      self.TileWidth = width
      self.TileHeight = height
      self.TileWidthNB = width - 2 * border
      self.TileHeightNB = height - 2 * border
      self.setup()


   def tileBuffer(self, format, type, image):
      self.TileFormat = format
      self.TileType = type
      self.TileBuffer = image


   def imageSize(self, width, height):
      self.ImageWidth = width
      self.ImageHeight = height
      self.setup()


   def imageBuffer(self, format, type, image):
      self.ImageFormat = format
      self.ImageType = type
      self.ImageBuffer = image


   def get(self, param):
      if param==TR_TILE_WIDTH:
         return self.TileWidth
      elif param==TR_TILE_HEIGHT:
         return self.TileHeight
      elif param==TR_TILE_BORDER:
         return self.TileBorder
      elif param==TR_IMAGE_WIDTH:
         return self.ImageWidth
      elif param==TR_IMAGE_HEIGHT:
         return self.ImageHeight
      elif param==TR_ROWS:
         return self.Rows
      elif param==TR_COLUMNS:
         return self.Columns
      elif param==TR_CURRENT_ROW:
         if (self.CurrentTile<0):
            return -1
         else:
            return self.CurrentRow
      elif param==TR_CURRENT_COLUMN:
         if (self.CurrentTile<0):
            return -1
         else:
            return self.CurrentColumn
      elif param==TR_CURRENT_TILE_WIDTH:
         return self.CurrentTileWidth
      elif param==TR_CURRENT_TILE_HEIGHT:
         return self.CurrentTileHeight
      elif param==TR_ROW_ORDER:
         return int(self.RowOrder)
      else:
         return 0


   def rowOrder(self, order):
      if order==TR_TOP_TO_BOTTOM or order==TR_BOTTOM_TO_TOP:
         self.RowOrder = order


   def ortho(self, left, right, bottom, top, zNear, zFar):
      self.Perspective = False
      self.Left = left
      self.Right = right
      self.Bottom = bottom
      self.Top = top
      self.Near = zNear
      self.Far = zFar


   def frustum(self, left, right, bottom, top, zNear, zFar):
      self.Perspective = True
      self.Left = left
      self.Right = right
      self.Bottom = bottom
      self.Top = top
      self.Near = zNear
      self.Far = zFar


   def perspective(self, fovy, aspect, zNear, zFar):
      from math import tan
      ymax = zNear * tan(fovy * 3.14159265 / 360.0)
      ymin = -ymax
      xmin = ymin * aspect
      xmax = ymax * aspect
      self.frustum(xmin, xmax, ymin, ymax, zNear, zFar)


   def beginTile(self):

      if self.CurrentTile <= 0:
         self.setup()
      # Save user's viewport, will be restored after last tile rendered
      self.ViewportSave = GL.glGetIntegerv(GL.GL_VIEWPORT)

      # which tile (by row and column) we're about to render
      if self.RowOrder==TR_BOTTOM_TO_TOP:
         self.CurrentRow = self.CurrentTile / self.Columns
         self.CurrentColumn = self.CurrentTile % self.Columns
      elif self.RowOrder==TR_TOP_TO_BOTTOM:
         self.CurrentRow = self.Rows - (self.CurrentTile / self.Columns) - 1
         self.CurrentColumn = self.CurrentTile % self.Columns
      else:
         raise RuntimeError

      assert(self.CurrentRow < self.Rows)
      assert(self.CurrentColumn < self.Columns)

      border = self.TileBorder

      # Compute actual size of this tile with border
      if self.CurrentRow < self.Rows-1:
         tileHeight = self.TileHeight
      else:
         tileHeight = self.ImageHeight - (self.Rows-1) * \
                      (self.TileHeightNB) + 2 * border

      if self.CurrentColumn < self.Columns-1:
         tileWidth = self.TileWidth
      else:
         tileWidth = self.ImageWidth - (self.Columns-1) * \
                     (self.TileWidthNB) + 2 * border

      # Save tile size, with border
      self.CurrentTileWidth = tileWidth
      self.CurrentTileHeight = tileHeight

      GL.glViewport(0, 0, tileWidth, tileHeight) #tile size including border

      # save current matrix mode
      matrixMode = GL.glGetIntegerv(GL.GL_MATRIX_MODE)[0]
      GL.glMatrixMode(GL.GL_PROJECTION)
      GL.glLoadIdentity()

      # compute projection parameters
      self.tileleft = left = self.Left + (self.Right - self.Left) \
                      * (self.CurrentColumn * self.TileWidthNB - border) \
                      / self.ImageWidth
      self.tileright = right = left + (self.Right - self.Left) * \
                       self.TileWidth / self.ImageWidth
      self.tilebottom = bottom = self.Bottom + (self.Top - self.Bottom) \
                        * (self.CurrentRow * self.TileHeightNB - border) / \
                        self.ImageHeight
      self.tiletop = top = bottom + (self.Top - self.Bottom) * self.TileHeight / \
                     self.ImageHeight

      if self.Perspective:
         GL.glFrustum(float(left), float(right),
                      float(bottom), float(top),
                      float(self.Near), float(self.Far))
      else:
         GL.glOrtho(float(left), float(right),
                    float(bottom), float(top),
                    float(self.Near), float(self.Far))

      # restore user's matrix mode
      GL.glMatrixMode(int(matrixMode))


   def endTile(self):
      assert(self.CurrentTile>=0)

      # be sure OpenGL rendering is finished
      GL.glFinish() #was glFlush()

      # save current glPixelStore values
      prevRowLength = GL.glGetIntegerv(GL.GL_PACK_ROW_LENGTH)[0]
      prevSkipRows = GL.glGetIntegerv(GL.GL_PACK_SKIP_ROWS)[0]
      prevSkipPixels = GL.glGetIntegerv(GL.GL_PACK_SKIP_PIXELS)[0]
      #prevAlignment = GL.glGetIntegerv(GL_PACK_ALIGNMENT)[0]

      if self.TileBuffer is not None:
         srcX = self.TileBorder
         srcY = self.TileBorder
         srcWidth = self.TileWidthNB
         srcHeight = self.TileHeightNB
         GL.glReadPixels(srcX, srcY, srcWidth, srcHeight,
                         self.TileFormat, self.TileType, self.TileBuffer)

      if self.ImageBuffer is not None:
         srcX = self.TileBorder
         srcY = self.TileBorder
         srcWidth = self.CurrentTileWidth - 2 * self.TileBorder
         srcHeight = self.CurrentTileHeight - 2 * self.TileBorder
         destX = self.TileWidthNB * self.CurrentColumn
         destY = self.TileHeightNB * self.CurrentRow

##          #save single tile
##          # setup pixel store for glReadPixels
##          GL.glPixelStorei(GL.GL_PACK_ROW_LENGTH, self.CurrentTileWidth)
##          GL.glPixelStorei(GL.GL_PACK_SKIP_ROWS, 0)
##          GL.glPixelStorei(GL.GL_PACK_SKIP_PIXELS, 0)
##          #GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)

##          # read the tile into buffer
##          gllib.glReadPixels(srcX, srcY, srcWidth, srcHeight,
##                             self.ImageFormat, self.ImageType, self.oneTileRenderBuffer)
##          import Image, sys
##          im = Image.fromstring('RGB', (srcWidth, srcHeight),
##                                self.oneTileRenderBuffer) 
##          if sys.platform!='win32':
##             im = im.transpose(Image.FLIP_TOP_BOTTOM)
##          im.save('tile%d_%d.tif'%(self.CurrentRow, self.CurrentColumn))
         
         # setup pixel store for glReadPixels
         GL.glPixelStorei(GL.GL_PACK_ROW_LENGTH, self.ImageWidth)
         GL.glPixelStorei(GL.GL_PACK_SKIP_ROWS, destY)
         GL.glPixelStorei(GL.GL_PACK_SKIP_PIXELS, destX)
         #GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)

         # read the tile into the final image
         #print 'DEBUG Tile renderer'
         #print 'column, row:', self.CurrentColumn, self.CurrentRow
         #print 'srcWidth,srcHeight,destX, destY',srcWidth,srcHeight,destX, destY
         gllib.glReadPixels(srcX, srcY, srcWidth, srcHeight,
                            self.ImageFormat, self.ImageType, self.ImageBuffer)

      # restore previous glPixelStore values
      GL.glPixelStorei(GL.GL_PACK_ROW_LENGTH, int(prevRowLength))
      GL.glPixelStorei(GL.GL_PACK_SKIP_ROWS, int(prevSkipRows))
      GL.glPixelStorei(GL.GL_PACK_SKIP_PIXELS, int(prevSkipPixels))
      #GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, prevAlignment)

      # increment tile counter, return 1 if more tiles left to render
      self.CurrentTile+=1
      if self.CurrentTile >= self.Rows * self.Columns:
         # restore user's viewport
         GL.glViewport(int(self.ViewportSave[0]), int(self.ViewportSave[1]),
                       int(self.ViewportSave[2]), int(self.ViewportSave[3]))
         self.CurrentTile = -1  # all done
         return 0
      else:
         return 1


   def trRasterPos3f(self, x, y, z):
      """
Replacement for glRastePos3f() which avoids the problem with invalid
raster pos.
"""
      if self.CurrentTile<0:
         # not doing tile rendering right now.  Let OpenGL do this.
         GL.glRasterPos3f(float(x), float(y), float(z))
      else:
         # Get modelview, projection and viewport
         modelview = GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX)
         proj = GL.glGetDoublev(GL.GL_PROJECTION_MATRIX)
         viewport = [0, 0, self.CurrentTileWidth, self.CurrentTileHeight]

         ## FIXME, need to finish this
         # Project object coord to window coordinate
##          projpoint = gluProject(x, y, z, modelview, proj, viewport)
##          if gluProject(x, y, z, modelview, proj, viewport, &winX, &winY, &winZ)){

##          /* set raster pos to window coord (0,0) */
##          glMatrixMode(GL_MODELVIEW)
##          glPushMatrix()
##          glLoadIdentity()
##          glMatrixMode(GL_PROJECTION)
##          glPushMatrix()
##          glLoadIdentity()
##          glOrtho(0.0, self.CurrentTileWidth,
##                  0.0, self.CurrentTileHeight, 0.0, 1.0)
##          glRasterPos3f(0.0, 0.0, -winZ)

##          /* Now use empty bitmap to adjust raster position to (winX,winY) */
##          {
##             GLubyte bitmap[1] = {0}
##             glBitmap(1, 1, 0.0, 0.0, winX, winY, bitmap)
##          }

##          /* restore original matrices */
##          glPopMatrix() /*proj*/
##          glMatrixMode(GL_MODELVIEW)
##          glPopMatrix()


