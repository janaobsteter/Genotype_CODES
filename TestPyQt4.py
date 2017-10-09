from PyQt4 import QtGui, QtCore, uic
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import sys

qtCreatorFile = '/home/jana/Genotipi/Genotipi_CODES/SelectionParameters.ui'
ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class OperatingWindow(QtGui.QMainWindow, ui_MainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        ui_MainWindow.__init__(self)

if __name__== "__main__":
    app = QtGui.QApplication(sys.argv)
    window = OperatingWindow()
    window.show()
