from PyQt4.QtCore import *
from PyQt4.QtGui import QPushButton, QVBoxLayout, QDialog, QApplication ,QMainWindow
from ui_mainwindow import Ui_MainWindow
from glwidget import GLWidget
class MainWindow(QMainWindow):
	def __init__(self, parent=None):
		super(MainWindow, self).__init__(parent)
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)
		self.gl_widget = GLWidget()
		self.ui.gl_layout.addWidget(self.gl_widget)
		#singal slot connect
		# self.connect(self.ui.apk1_open, SIGNAL('clicked()'),self.apk1_open_onclicked)
	# @pyqtSlot()
	# def apk1_open_onclicked(self):
	# 	self.emit(SIGNAL('open_apk1'))