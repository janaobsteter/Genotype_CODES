# -*- coding: utf-8 -*-

import sys
from PyQt4.QtCore import *
from PyQt4.QtGui import *

class MoviePlayer(QWidget):
    def __init__(self, gif1, gif2, parent=None):
        super(MoviePlayer, self).__init__(parent)

        self.setGeometry(200, 200, 400, 400)
        self.setWindowTitle("Odgovori")

        self.movie_screen = QLabel()
        self.movie_screen.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.movie_screen.setAlignment(Qt.AlignCenter)
        self.gif1 = gif1
        self.gif2 = gif2

        btn_start = QLabel()
        btn_start.setText("Kaj gleda tvoje sefe?")


        self.entered = QLineEdit()
        self.entered.returnPressed.connect(self.start)
        self.entered.textEdited.connect(self.stop)


        main_layout = QVBoxLayout()
        main_layout.addWidget(self.movie_screen)
        main_layout.addWidget(btn_start)
        main_layout.addWidget(self.entered)
        self.setLayout(main_layout)






    def start(self):
        """
        Start animation
        """
        if str(self.entered.text().toUtf8()).upper() == "KURC":
            self.movie = QMovie(self.gif1, QByteArray(), self)
            self.movie.setCacheMode(QMovie.CacheAll)
            self.movie.setSpeed(100)
            self.movie_screen.setMovie(self.movie)
        if str(self.entered.text().toUtf8()).upper() != "KURC":
            self.movie = QMovie(self.gif2, QByteArray(), self)
            self.movie.setCacheMode(QMovie.CacheAll)
            self.movie.setSpeed(100)
            self.movie_screen.setMovie(self.movie)

        self.movie.start()

    def stop(self):
        """
        Stop the animation
        """
        if str(self.entered.text().toUtf8()).upper() != "kurc":
            self.movie = QMovie(self.gif1, QByteArray(), self)
            self.movie.setCacheMode(QMovie.CacheAll)
            self.movie.setSpeed(100)
            self.movie_screen.setMovie(self.movie)

        self.movie.stop()


app = QApplication(sys.argv)
player = MoviePlayer("/home/jana/Downloads/tenor.gif", "/home/jana/Downloads/tryagain.gif" )
player.show()
sys.exit(app.exec_())