# -*- coding: utf-8 -*-
import sys
from PyQt4.QtCore import *
from PyQt4.QtGui import QApplication,QFileDialog
from mainwindow import MainWindow
class Application(QObject):
	def Go(self,argv):
		app=QApplication(argv)
		self.mainwindow=MainWindow()
		self.__connect_mainwindow()
		self.mainwindow.show()
		sys.exit(app.exec_())
	def __connect_mainwindow(self):
		print "connect MainWindow"
