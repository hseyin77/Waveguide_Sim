import sys
from math import sqrt, pi
from turtle import st
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QLineEdit, QPushButton, QRadioButton, QVBoxLayout
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QLabel, QLineEdit
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import special

# Define the cutoff frequency functions for TE and TM modes
def calc_te_cutoff(a, b, m, n):
    return c / sqrt(mu_0 * epsilon_0) * sqrt((m / a)**2 + (n / b)**2)

def calc_tm_cutoff(a, b, m, n):
    return c / sqrt(mu_0 / epsilon_0) * sqrt((m / a)**2 + (n / b)**2)
PI = math.pi
def TE_plot(m, n, a, b):
    x = np.linspace(0, a, 101)
    y = np.linspace(0, b, 101)
    X, Y = np.meshgrid(x, y)
    if m==0 and n ==0:
        st.error("m and n cannot be 0 at the same time")
        return
    u = np.cos(m * PI / a * X) * np.sin(n * PI / b * Y)
    v = -1 * np.sin(m * PI / a * X) * np.cos(n * PI / b * Y)
    fig, ax = plt.subplots()
    plt.streamplot(X, Y, u, v, color="xkcd:azure")
    plt.axis("scaled")
    plt.suptitle("E field")
    plt.xlim(0, a)
    plt.ylim(0, b)
    plt.show

    
    u = np.sin(m * PI / a * X) * np.cos(n * PI / b * Y)
    v = np.cos(m * PI / a * X) * np.sin(n * PI / b * Y)
    fig, ax = plt.subplots()
    plt.streamplot(x, y, u, v, color="red")
    plt.axis("scaled")
    plt.suptitle("H field")
    plt.xlim(0, a)
    plt.ylim(0, b)
    plt.show()
    

# Define the constants
c = 3 * 10**8  # speed of light in m/s
mu_0 = 4 * pi * 10**(-7)  # permeability of free space in Tm/A
epsilon_0 = 1 / (mu_0 * c**2)  # permittivity of free space in F/m

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Cutoff Frequency Calculator")
        self.setFixedSize(500, 500)

        pixmap = QPixmap("C:\\Users\Huseyin Sevuk\\Downloads\\space + time (1).png")
        self.background_label = QLabel(self)
        self.background_label.setPixmap(pixmap)
        self.background_label.setGeometry(0, 0, 500, 500)

        # Create the input boxes for a and b
        self.a_label = QLabel("a (mm):", self)
        self.a_label.move(40, 55)
        self.a_entry = QLineEdit(self)
        self.a_entry.move(110, 55)
        self.a_entry.setFixedWidth(100)
        self.round_corners(self.a_entry)

        self.b_label = QLabel("b (mm):", self)
        self.b_label.move(40, 105)
        self.b_entry = QLineEdit(self)
        self.b_entry.move(110, 105)
        self.b_entry.setFixedWidth(100)
        self.round_corners(self.b_entry)

        # Create the selection button for TE and TM modes
        self.mode_label = QLabel("Choose mode:", self)
        self.mode_label.move(350, 50)
        self.te_button = QRadioButton("TE", self,)
        self.te_button.move(345, 80)
        self.tm_button = QRadioButton("TM", self)
        self.tm_button.move(395, 80)
        self.te_button.setChecked(True)

        # Create the input boxes for m and n
        self.m_label = QLabel("m:", self)
        self.m_label.move(335, 110)
        self.m_entry = QLineEdit(self)
        self.m_entry.move(360, 110)
        self.m_entry.setFixedWidth(25)
        self.round_corners(self.m_entry)

        self.n_label = QLabel("n:", self)
        self.n_label.move(400, 110)
        self.n_entry = QLineEdit(self)
        self.n_entry.move(420, 110)
        self.n_entry.setFixedWidth(25)
        self.round_corners(self.n_entry)

        # Create the button to trigger the calculation
        self.calc_button = QPushButton("Calculate", self)
        self.calc_button.move(200, 345)
        self.calc_button.setFixedWidth(100)
        self.calc_button.setFixedHeight(30)
        self.calc_button.setStyleSheet("border-radius: 10px; background-color: #6777d1; color: white; font-family = Berlin Sans FB; font-size: 15px; font-weight: bold;")
        self.calc_button.clicked.connect(self.calc_cutoff)
        

        # Create the result label
        self.result_label = QLabel("", self)
        self.result_label.move(120, 430)
        self.result_label.setFixedWidth(300)

        self.set_font_properties()

    def set_font_properties(self):
        # Set font properties for labels
        font = self.font()
        font.setFamily("Berlin Sans FB")      
        font.setPointSize(11)

        self.a_label.setFont(font)
        self.b_label.setFont(font)
        self.mode_label.setFont(font)
        self.m_label.setFont(font)
        self.n_label.setFont(font)
        self.result_label.setFont(font)

        self.a_label.setStyleSheet("color: #002060;")
        self.b_label.setStyleSheet("color: #002060;")
        self.mode_label.setStyleSheet("color: #002060;")
        self.m_label.setStyleSheet("color: #002060;")
        self.n_label.setStyleSheet("color: #002060;")
        self.result_label.setStyleSheet("color: red;")

    def round_corners(self, widget):
        widget.setStyleSheet("border-radius: 10px;")

    def calc_cutoff(self):
        a = float(self.a_entry.text())
        b = float(self.b_entry.text())
        m = int(self.m_entry.text())
        n = int(self.n_entry.text())
        mode = "TE" if self.te_button.isChecked() else "TM"
        
        if mode == "TE":
            freq = calc_te_cutoff(a, b, m, n)
            TE_plot(a, b, m, n)  
        elif mode == "TM":
            freq = calc_tm_cutoff(a, b, m, n)
        self.result_label.setText(f"Cutoff frequency: {freq:.2f} GHz")
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
