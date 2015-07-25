from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtGui import QVector2D
from OpenGL import GL
from OpenGL import GLU
from water2d import Water2d
class GLWidget(QtOpenGL.QGLWidget):
    def __init__(self, parent=None):
        super(GLWidget, self).__init__(parent)
        self.lastPos = QtCore.QPoint()
        self.water2d = Water2d(100,50,100)
        self.startTimer(40)
    def timerEvent(self, event):
        self.water2d.Update(0.01)
        self.update()
    def initializeGL(self):
        self.qglClearColor(QtGui.QColor(255, 255, 255, 255))
        GL.glShadeModel(GL.GL_FLAT)
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_CULL_FACE)

    def paintGL(self):
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glLoadIdentity()
        GL.glColor3f(0,0,0)
        self.water2d.Render()
    def resizeGL(self, width, height):
        # //: means floor the divede result
        self.width=width
        self.height=height
        GL.glViewport(0, 0, width, height)
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        GLU.gluOrtho2D(0,100,0,100);
        GL.glMatrixMode(GL.GL_MODELVIEW)

    def mousePressEvent(self, event):
        self.lastPos = event.pos()
        x=event.x()*100.0/self.width
        y=(self.height-event.y())*100.0/self.height
        if event.buttons() & QtCore.Qt.RightButton:
            self.water2d.AddBox(QVector2D(x,y))
    def mouseMoveEvent(self, event):
        dx = event.x() - self.lastPos.x()
        dy = event.y() - self.lastPos.y()
        self.lastPos=event.pos()
        if event.buttons() & QtCore.Qt.LeftButton:
            x=event.x()
            y=event.y()
            self.water2d.HandleMouseMove(x*100.0/self.width,y*100.0/self.height,dx,dy)
