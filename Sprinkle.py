# -*- coding: utf-8 -*-

import sys
from PyQt4.QtCore import *
from PyQt4.QtGui import *

class MoviePlayer(QWidget):
    def __init__(self, gif1, gif2, parent=None):
        super(MoviePlayer, self).__init__(parent)

        self.setGeometry(200, 200, 400, 400)
        self.setWindowTitle("KDO SI????")

        self.movie_screen = QLabel()
        self.movie_screen.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.movie_screen.setAlignment(Qt.AlignCenter)
        self.gif1 = gif1
        self.gif2 = gif2

        btn_start = QLabel()
        btn_start.setText("Kdo si ti?")


        self.entered = QLineEdit()
        self.entered.returnPressed.connect(self.start)
        self.entered.textChanged.connect(self.stop)


        main_layout = QVBoxLayout()
        main_layout.addWidget(self.movie_screen)
        main_layout.addWidget(btn_start)
        main_layout.addWidget(self.entered)
        self.setLayout(main_layout)






    def start(self):
        """
        Start animation
        """
        if str(self.entered.text().toUtf8()).upper() == "MAVRICNA MRVICA":
            self.movie = QMovie(self.gif1, QByteArray(), self)
            self.movie.setCacheMode(QMovie.CacheAll)
            self.movie.setSpeed(100)
            self.movie_screen.setMovie(self.movie)
        if str(self.entered.text().toUtf8()).upper() != "MAVRICNA MRVICA":
            self.movie = QMovie(self.gif2, QByteArray(), self)
            self.movie.setCacheMode(QMovie.CacheAll)
            self.movie.setSpeed(100)
            self.movie_screen.setMovie(self.movie)
        self.movie.start()

    def stop(self):
        """
        Stop the animation
        """
        if str(self.entered.text().toUtf8()).upper() == "MAVRICNA MRVICA":
            self.movie = QMovie(self.gif1, QByteArray(), self)
            self.movie.setCacheMode(QMovie.CacheAll)
            self.movie.setSpeed(100)
            self.movie_screen.setMovie(self.movie)
        if str(self.entered.text().toUtf8()).upper() != "MAVRICNA MRVICA":
            self.movie = QMovie(self.gif2, QByteArray(), self)
            self.movie.setCacheMode(QMovie.CacheAll)
            self.movie.setSpeed(100)
            self.movie_screen.setMovie(self.movie)
        self.movie.stop()


app = QApplication(sys.argv)
player = MoviePlayer("/home/jana/Genotipi/Genotipi_CODES/sprinkle.gif", "/home/jana/Genotipi/Genotipi_CODES/No.gif", )
player.show()
sys.exit(app.exec_())