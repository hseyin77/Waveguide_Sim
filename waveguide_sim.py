import sys
from math import sqrt, pi
from turtle import st
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QLineEdit, QPushButton, QRadioButton, QVBoxLayout, QButtonGroup, QDesktopWidget, QMessageBox
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QLabel, QLineEdit
from PyQt5.QtGui import QPixmap, QMovie
from PyQt5.QtCore import Qt, QTimer
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import jv, jn_zeros
PI = math.pi


c = 3 * 10**8  # m/s cinsinden ışık hızı
def show_error_and_reset_entries(msg_text, *entries):
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setText(msg_text)
    msg.setWindowTitle("Error")
    msg.exec_()
    for entry in entries:
        entry.clear()


# dairesel ve dikdörtgen dalga kılavuzuna ait kesim frekansı formüllerini içeren fonksiyonlar
def calc_circ_cutoff(a):
    return ((1.8412*c) / (2*PI*a)) / 10000000


def calc_te_cutoff(a, b, m, n, mu_0, epsilon_0):
   
    if mu_0 == 0 and epsilon_0 == 0:
        return((c / 2) * sqrt((m / a)**2 + (n / b)**2)) / 10000000
    else:
        return(c / (2*sqrt(mu_0 * epsilon_0)) * sqrt((m / a)**2 + (n / b)**2)) / 10000000

def calc_tm_cutoff(a, b, m, n,  mu_0, epsilon_0):
    if mu_0 == 0 and epsilon_0 == 0:
        return((c / 2) * sqrt((m / a)**2 + (n / b)**2)) / 10000000
    else:
        return(c / (2*sqrt(mu_0 * epsilon_0)) * sqrt((m / a)**2 + (n / b)**2)) / 10000000
def circular_plot_TE(a, m, n):
    r = np.linspace(0, a, 101)
    t = np.linspace(0, 2*PI, 101)
    T, RAD = np.meshgrid(t, r)

    X = special.jnp_zeros(m, n) 
    U = special.jv(m, X[-1].round(3) / a * RAD) * np.sin(m * T)
    V = special.jvp(m, X[-1].round(3) / a * RAD) * np.cos(m * T)
    plt.axis("scaled")
    fig, ax = plt.subplots()
    #plt.polar(2 * PI * a)
    plt.streamplot(T, RAD, V, U, color="xkcd:azure")
    plt.axis("scaled")
    plt.suptitle("E field")
    plt.show()

    
    U = -1 * special.jv(m, X[-1].round(3) / a * RAD) * np.cos(m * T)
    V = special.jv(m, X[-1].round(3) / a * RAD) * np.sin(m * T)
    fig, ax = plt.subplots()
    #plt.polar(2*PI*a)
    plt.streamplot(T, RAD, V, U, color="red")
    plt.axis("scaled")
    plt.suptitle("H field")
    plt.show()    
#TE ve TM modlarına ait grafiklerin bastırılması için gerekli formülleri içeren fonsiyonlar
def circular_plot_TM(a, m, n):
     # Constants
    E0 = 1  # Electric field magnitude

    # Calculate the first zero of the Bessel function for the TM mode
    p = n  # Root number
    bessel_zeros = jn_zeros(m, p)
    kc = bessel_zeros[n - 1] / a  # Cutoff wavenumber

    # Create the grid
    resolution = 50  # Number of points in the grid
    r = np.linspace(0, a, resolution)
    theta = np.linspace(0, 2 * np.pi, resolution)
    R, Theta = np.meshgrid(r, theta)

    # Calculate the electric field component for the TM mode
    Er = E0 * jv(m, kc * R)

    # Normalize the vectors for uniform arrow length
    Er_magnitude = np.sqrt(Er**2)
    Er_normalized = Er / Er_magnitude

    # Convert to Cartesian coordinates for plotting
    X = R * np.cos(Theta)
    Y = R * np.sin(Theta)

    # Plot the electric field distribution with uniformly sized arrows
    plt.figure(figsize=(8, 8))
    plt.quiver(X[::5, ::5], Y[::5, ::5], Er_normalized[::5, ::5] * np.cos(Theta[::5, ::5]), Er_normalized[::5, ::5] * np.sin(Theta[::5, ::5]), color='r', scale=20)
    plt.title(f'Electric Field Distribution for TM{m}{n} Mode (XY Plane)')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.show()

    # Calculate the magnetic field component for the TM mode
    H_phi = (1j * kc * E0 / (kc**2)) * jv(m, kc * R)
    H_phi_magnitude = np.abs(H_phi)
    H_phi_normalized = H_phi / H_phi_magnitude

    # Plot the magnetic field distribution with uniformly sized arrows
    plt.figure(figsize=(8, 8))
    plt.quiver(X[::5, ::5], Y[::5, ::5], -H_phi_normalized[::5, ::5] * np.cos(Theta[::5, ::5]), -H_phi_normalized[::5, ::5] * np.sin(Theta[::5, ::5]), color='b', scale=20)
    plt.title(f'Magnetic Field Distribution for TM{m}{n} Mode (XY Plane)')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.show()
    

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
    
def TM_plot(m, n, a, b):
    x = np.linspace(0, a, 101)
    y = np.linspace(0, b, 101)
    X, Y = np.meshgrid(x, y)
    if m==0 and n ==0:
        st.error("m and n cannot be 0 at the same time!")
        return
    u = np.cos(m * PI / a * X) * np.sin(n * PI / b * Y)
    v = np.sin(m * PI / a * X) * np.cos(n * PI / b * Y)
    
    fig, ax = plt.subplots()
    plt.streamplot(X, Y, u, v, color="xkcd:azure")
    plt.axis("scaled")
    plt.suptitle("E field")
    plt.xlim(0, a)
    plt.ylim(0, b)
    plt.show

    
    u = np.sin(m * PI / a * X) * np.cos(n * PI / b * Y)
    v = -1 * np.cos(m * PI / a * X) * np.sin(n * PI / b * Y)
    fig, ax = plt.subplots()
    plt.streamplot(x, y, u, v, color="red")
    plt.axis("scaled")
    plt.suptitle("H field")
    plt.xlim(0, a)
    plt.ylim(0, b)
    plt.show()
#program başlatılırken oluşan yükleme ekranına ait sınıf
class LoadingScreen(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Loading...")
        self.setFixedSize(500, 500)
        self.setWindowFlags(Qt.FramelessWindowHint)
        
        
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())


        # başlangıç ekranının arka plan resmini ekledim
        self.background_label = QLabel(self)
        pixmap = QPixmap("intro.png") 
        self.background_label.setPixmap(pixmap)
        self.background_label.setGeometry(0, 0, 500, 500)
        self.background_label.setScaledContents(True) 

        self.movie = QMovie("Cube@1x-1.4s-200px-200px.gif")
        self.loading_label = QLabel(self)
        self.loading_label.setMovie(self.movie)
        self.loading_label.move(150,150)
        self.movie.start()
#programın ana penceresi
class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Waveguide Designer")
        self.setFixedSize(500, 500)

        pixmap = QPixmap("arkaplan.png")
        self.background_label = QLabel(self)
        self.background_label.setPixmap(pixmap)
        self.background_label.setGeometry(0, 0, 500, 500)

        

        #burada kullanıcıdan giriş alacağımız input kutularının oluşturulması ve konumlandırılmasına ait parametreler mevcut
        self.a_label = QLabel("a (cm):", self)
        self.a_label.move(40, 55)
        self.a_entry = QLineEdit(self)
        self.a_entry.move(100, 55)
        self.a_entry.setFixedWidth(100)
        self.round_corners(self.a_entry)

        self.b_label = QLabel("b (cm):", self)
        self.b_label.move(40, 105)
        self.b_entry = QLineEdit(self)
        self.b_entry.move(100, 105)
        self.b_entry.setFixedWidth(100)
        self.round_corners(self.b_entry)

        #mod seçimi
        self.mode_label = QLabel("Operation mode", self)
        self.mode_label.move(335, 50)
        self.mode_group = QButtonGroup(self)  # Button group for operation mode
        self.te_button = QRadioButton("TE", self)
        self.te_button.move(345, 80)
        self.tm_button = QRadioButton("TM", self)
        self.tm_button.move(395, 80)
        self.mode_group.addButton(self.te_button)
        self.mode_group.addButton(self.tm_button)
        self.te_button.setChecked(True)

        #mod seçimi
        self.waveguide_label = QLabel("Waveguide Type", self)
        self.waveguide_label.move(335, 150)
        self.waveguide_group = QButtonGroup(self) 
        self.rectangular_button = QRadioButton("Rectangular", self)
        self.rectangular_button.move(350, 180)
        self.circular_button = QRadioButton("Circular", self)
        self.circular_button.move(350, 200)
        self.waveguide_group.addButton(self.rectangular_button)
        self.waveguide_group.addButton(self.circular_button)

        self.r_label = QLabel("**a label keeps the R value for circular waveguide", self)
        self.r_label.move(40, 240)
        self.warning_label = QLabel("**In Circular mode can only visualize for TM01 mode",self)
        self.warning_label.move(40,260)
        
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

        
        self.mu_label = QLabel("μ Permeability:", self)
        self.mu_label.move(40, 155)
        self.mu_entry = QLineEdit(self)
        self.mu_entry.move(160, 155)
        self.mu_entry.setFixedWidth(50)
        self.round_corners(self.mu_entry)

        self.epsilon_label = QLabel("ε Permittivity:", self)
        self.epsilon_label.move(40, 205)
        self.epsilon_entry = QLineEdit(self)
        self.epsilon_entry.move(160, 205)
        self.epsilon_entry.setFixedWidth(50)
        self.round_corners(self.epsilon_entry)
        

        #hesapla butonunu oluşturdum
        self.calc_button = QPushButton("Calculate", self)
        self.calc_button.move(200, 345)
        self.calc_button.setFixedWidth(100)
        self.calc_button.setFixedHeight(30)
        self.calc_button.setStyleSheet("border-radius: 10px; background-color: #6777d1; color: white; font-family = Berlin Sans FB; font-size: 15px; font-weight: bold;")
        self.calc_button.clicked.connect(self.calc_cutoff)
        
        #reset butonu
        self.reset_button = QPushButton("Reset", self)
        self.reset_button.move(350, 15)
        self.reset_button.setFixedWidth(100)
        self.reset_button.setFixedHeight(30)
        self.reset_button.setStyleSheet(
            "border-radius: 10px; background-color: #d16a6a; color: white; font-family = Berlin Sans FB; font-size: 15px; font-weight: bold;")
        self.reset_button.clicked.connect(self.reset_fields)

        #sonucun gösterileceği pencere
        self.result_label = QLabel("", self)
        self.result_label.move(150, 420)
        self.result_label.setFixedWidth(300)
        

        self.set_font_properties()

    def set_font_properties(self):
        
        font = self.font()
        font.setFamily("Berlin Sans FB")      
        font.setPointSize(11)
        
        self.r_label.setFont(font)
        self.a_label.setFont(font)
        self.b_label.setFont(font)
        self.mode_label.setFont(font)
        self.waveguide_label.setFont(font)
        self.m_label.setFont(font)
        self.n_label.setFont(font)
        self.result_label.setFont(font)
        self.epsilon_label.setFont(font)
        self.mu_label.setFont(font)
        self.warning_label.setFont(font)

        self.warning_label.setStyleSheet("color: #002060;")
        self.r_label.setStyleSheet("color: #002060;")
        self.a_label.setStyleSheet("color: #002060;")
        self.b_label.setStyleSheet("color: #002060;")
        self.mode_label.setStyleSheet("color: #002060;")
        self.waveguide_label.setStyleSheet("color: #002060;")
        self.m_label.setStyleSheet("color: #002060;")
        self.n_label.setStyleSheet("color: #002060;")
        self.result_label.setStyleSheet("color: white;")
        self.epsilon_label.setStyleSheet("color: #002060;")
        self.mu_label.setStyleSheet("color: #002060;")


    def round_corners(self, widget):
        widget.setStyleSheet("border-radius: 10px;")

    def reset_fields(self):
        self.a_entry.clear()
        self.b_entry.clear()
        self.m_entry.clear()
        self.n_entry.clear()
        self.mu_entry.clear()
        self.epsilon_entry.clear()
        self.result_label.clear()
        plt.close('all')

    def calc_cutoff(self):
        
        waveguide = "Rectangular" if self.rectangular_button.isChecked() else "Circular"
        if waveguide == "Circular":
            a = float(self.a_entry.text())
            n = int(self.n_entry.text())
            m = int(self.m_entry.text())

            if m ==0 and n == 0:
                show_error_and_reset_entries("m and n cannot be 0 at the same time!", self.m_entry, self.n_entry)
                return
            elif n == 0:
                show_error_and_reset_entries("n cannot be 0!", self.n_entry)
                return
            freq = calc_circ_cutoff(a)
            mode = "TE" if self.te_button.isChecked() else "TM"
            if mode == "TE":
                circular_plot_TE(a, m, n)
            else:
                circular_plot_TM(a, m, n)
            
        elif waveguide == "Rectangular":
            a = float(self.a_entry.text())
            b = float(self.b_entry.text())
            m = int(self.m_entry.text())
            n = int(self.n_entry.text())
            mu_0 = float(self.mu_entry.text())
            epsilon_0 = float(self.epsilon_entry.text())
            if a == 0 or b == 0:
                show_error_and_reset_entries("a and b cannot be 0!", self.a_entry, self.b_entry)
                return
            if m == 0 and n == 0:
                show_error_and_reset_entries("m and n cannot be 0 at the same time!", self.m_entry, self.n_entry)
                return
            mode = "TE" if self.te_button.isChecked() else "TM"
            if mode == "TE":
                freq = calc_te_cutoff(a, b, m, n, mu_0, epsilon_0)
                TE_plot(m, n, a, b)  
            else:
                freq = calc_tm_cutoff(a, b, m, n, mu_0, epsilon_0)
                TM_plot(m, n, a, b)
        self.result_label.setText(f"Cutoff frequency: {freq:.2f} GHz")
if __name__ == "__main__":
    
    app = QApplication(sys.argv)
    
    #yükleme ekranını göster
    loading_screen = LoadingScreen()
    loading_screen.show()

    # ana pencereyi aç
    main_window = MainWindow()

    def show_main_window():
        loading_screen.close()
        main_window.show()

    # yükleme ekraının ne kadar açık kalacağını göster
    QTimer.singleShot(10000, show_main_window)

    sys.exit(app.exec_())
