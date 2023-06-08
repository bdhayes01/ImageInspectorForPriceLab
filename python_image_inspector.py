# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 15:58:09 2021

@author: wtd14
"""

import os
import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import (QApplication, QFileDialog, QListWidgetItem, QMainWindow)
import PyQt5.uic as uic
from PyQt5 import QtCore
import numpy as np
from scipy import signal
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from copy import copy
import math
import csv
from matplotlib.figure import Figure

# imported file
import ROI_class as roi  # ROI_class.py

MW_width = 1491
MW_height = 951
isIM = None


def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)


mainWindow_ui_path = resource_path("python_image_inspector_layout.ui")
mmcWindow_ui_path = resource_path("python_image_inspector_multicomp.ui")
DESI_ID_path = resource_path("DESI_ID_File.csv")

loaded_ui_multicomp = uic.loadUiType(mmcWindow_ui_path)[0]


class MultiMapCompareobject(QMainWindow, loaded_ui_multicomp):
    def __init__(self):  # Initialization of the code
        super(MultiMapCompareobject, self).__init__()
        self.setWindowTitle("Multi Map Compare")
        self.setFixedWidth(int(MW_width * 0.9))
        self.setFixedHeight(871)
        self.setupUi(self)

        self.Map_listbox_mmcWindow.itemClicked.connect(self.MultiMapCompare_Map_listbox_pickitem)

        self.pickitem = ''

        # canvas
        self.map_canvas1 = None
        self.map_canvas2 = None
        self.map_canvas3 = None
        self.map_canvas4 = None
        self.map_canvas5 = None
        self.map_canvas6 = None

        # plots
        self.map_ax1 = None
        self.map_ax2 = None
        self.map_ax3 = None
        self.map_ax4 = None
        self.map_ax5 = None
        self.map_ax6 = None

        # images
        self.map_img1 = None
        self.map_img2 = None
        self.map_img3 = None
        self.map_img4 = None
        self.map_img5 = None
        self.map_img6 = None

        # toolbars
        self.map_toolbar1 = None
        self.map_toolbar2 = None
        self.map_toolbar3 = None
        self.map_toolbar4 = None
        self.map_toolbar5 = None
        self.map_toolbar6 = None

        # colorbars
        self.map_cbar1 = None
        self.map_cbar2 = None
        self.map_cbar3 = None
        self.map_cbar4 = None
        self.map_cbar5 = None
        self.map_cbar6 = None

        self.map_packet = {
            1: [self.map_canvas1, self.map_ax1, self.map_img1, self.map_toolbar1, self.map_cbar1, self.slot1_plot,
                self.slot1_name],
            2: [self.map_canvas2, self.map_ax2, self.map_img2, self.map_toolbar2, self.map_cbar2, self.slot2_plot,
                self.slot2_name],
            3: [self.map_canvas3, self.map_ax3, self.map_img3, self.map_toolbar3, self.map_cbar3, self.slot3_plot,
                self.slot3_name],
            4: [self.map_canvas4, self.map_ax4, self.map_img4, self.map_toolbar4, self.map_cbar4, self.slot4_plot,
                self.slot4_name],
            5: [self.map_canvas5, self.map_ax5, self.map_img5, self.map_toolbar5, self.map_cbar5, self.slot5_plot,
                self.slot5_name],
            6: [self.map_canvas6, self.map_ax6, self.map_img6, self.map_toolbar6, self.map_cbar6, self.slot6_plot,
                self.slot6_name]}

    def MultiMapCompare_Map_listbox_pickitem(self, item):
        self.pickitem = item

    def display_mmcWindow(self):
        self.show()


loaded_ui_main = uic.loadUiType(mainWindow_ui_path)[0]


class MainGUIobject(QtWidgets.QMainWindow, loaded_ui_main):
    def __init__(self, parent=None):  # Initialization of the code
        QtWidgets.QMainWindow.__init__(self, parent)
        super(MainGUIobject, self).__init__()
        self.annotation = None
        self.mzVals = None
        self.intensity = None
        self.drifts = None
        self.chosenDataIso = None
        self.view = None
        self.viewPlusOne = None
        self.viewPlusTwo = None
        self.con_cbar = None
        self.setupUi(self)
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowSystemMenuHint)

        # load in default ID file
        self.ids_pd = pd.read_csv(DESI_ID_path)

        # multi map compare instance
        self.mmcWindow = MultiMapCompareobject()

        # Function connection
        self.find_file.clicked.connect(self.find_file_Callback)
        self.find_file_mzOI.clicked.connect(self.find_file_mzOI_Callback)
        self.start_cube.clicked.connect(self.start_cube_Callback)
        self.pick_point.clicked.connect(self.pick_point_Callback)
        self.mass_up.clicked.connect(self.mass_up_Callback)
        self.mass_down.clicked.connect(self.mass_down_Callback)
        self.start.returnPressed.connect(self.start_Callback)
        self.msindex.returnPressed.connect(self.msindex_Callback)
        self.zmax.sliderMoved.connect(self.zmax_Callback)
        self.zmax.valueChanged.connect(self.zmax_Callback)
        self.temp_max.returnPressed.connect(self.temp_max_Callback)
        self.zmax_isotope.sliderMoved.connect(self.zmax_isotope_Callback)
        self.zmax_isotope.valueChanged.connect(self.zmax_isotope_Callback)
        self.zmin_isotope.sliderMoved.connect(self.zmin_isotope_Callback)
        self.zmin_isotope.valueChanged.connect(self.zmin_isotope_Callback)
        self.ROI_select.clicked.connect(self.ROI_select_Callback_mask)
        self.ROI_process.clicked.connect(self.ROI_select_Callback_process)
        self.exportROI.clicked.connect(self.exportROI_Callback)
        self.ROI_listbox.itemDoubleClicked.connect(self.ROI_listbox_Callback)
        self.importROI.clicked.connect(self.importROI_Callback)
        self.deleteROIbutton.clicked.connect(self.deleteROIbutton_Callback)
        self.exportROI_val.clicked.connect(self.exportROI_spectra_val_Callback)
        self.find_IDlist.clicked.connect(self.find_IDlist_Callback)
        self.Map_listbox.itemDoubleClicked.connect(self.Map_listbox_Callback)
        self.exportConcMap.clicked.connect(self.export_ConcMap_Callback)
        self.exportIsotopeMap.clicked.connect(self.export_IsotopeMap_Callback)
        self.deleteMapbutton.clicked.connect(self.deleteMapbutton_Callback)
        self.clearMapListboxbutton.clicked.connect(self.clearMapbutton_Callback)
        self.extract_Map_mzOI.clicked.connect(self.mzOI_extractMap_Callback)

        # imButton = self.find_el
        self.IMDataButton.clicked.connect(self.setIM)
        self.MSDataButton.clicked.connect(self.setMS)
        # self.MSDataButton.clicked(isIM = False)
        # self.IMDataButton.clicked(isIM = True)

        self.multiMW.clicked.connect(self.MultiMapCompare_Display_Callback)
        # self.mmcWindow.slot1_load.clicked.connect(self.MultiMapCompare_LoadMap_Callback)
        self.mmcWindow.slot1_load.clicked.connect(
            lambda: self.MultiMapCompare_LoadMap_Callback(self.mmcWindow.slot1_load))
        self.mmcWindow.slot2_load.clicked.connect(
            lambda: self.MultiMapCompare_LoadMap_Callback(self.mmcWindow.slot2_load))
        self.mmcWindow.slot3_load.clicked.connect(
            lambda: self.MultiMapCompare_LoadMap_Callback(self.mmcWindow.slot3_load))
        self.mmcWindow.slot4_load.clicked.connect(
            lambda: self.MultiMapCompare_LoadMap_Callback(self.mmcWindow.slot4_load))
        self.mmcWindow.slot5_load.clicked.connect(
            lambda: self.MultiMapCompare_LoadMap_Callback(self.mmcWindow.slot5_load))
        self.mmcWindow.slot6_load.clicked.connect(
            lambda: self.MultiMapCompare_LoadMap_Callback(self.mmcWindow.slot6_load))
        self.mmcWindow.exportMapData.clicked.connect(self.MultiMapCompare_exportMapData_Callback)

        # function image_inspector_OpeningFcn
        self.has_data = 0
        self.ptflag = False
        self.rgflag = False
        self.rgflagROI = False
        self.rgflagROIexportMS = False
        self.exSpecflag = False
        self.exConcflag = False
        self.exIsotopeflag = False
        self.fName_mzOI_flag = False
        self.micrometer.setChecked(True)
        self.massplusone.setChecked(True)
        self.ideal_ratio.setText('1')
        QListWidgetItem('ROI list appears here', self.ROI_listbox)
        QListWidgetItem('Map list appears here', self.Map_listbox)

        self.fName = ""
        self.fName_mzOI = ""
        self.cubefilename = ""
        self.ROI = {}
        self.ROIplots = {}
        self.ROI_img_mean = {}
        self.ROIcount = 0
        self.Maps = {}
        self.Mapcount = 0
        self.Map_listselect_text = ""
        self.x = 0
        self.x_end = 0
        self.y = 0
        self.y_end = 0
        self.z = 0
        self.intens = 0
        self.a = 0
        self.b = 0
        self.c = 0
        self.img_mean = 0
        self.reshaped = 0
        self.PS_Peak_Intensity = 0
        self.PS_Sum_Intensity = 0
        self.PS_Index = 0
        self.PS_Mean_Intensity = 0
        self.Background = 0
        self.Noise = 0
        self.img_std = 0
        self.z_max = 0
        self.z_min = 0
        self.t_max = 0
        self.areas = 0
        self.includemassplustwo = False
        self.mintratio = 0
        self.good = 0
        self.checkpoint_maxima = 10
        self.err_multp = 1e6
        self.mzOI_index = 0

        self.resolution = 1000

        self.mintratio = 0
        self.mnormintratio = 0
        self.mintratiowithmassplustwo = 0
        self.mnormintratiowithmassplustwo = 0
        self.mplusoneintratio = 0
        self.mplusonenormintratio = 0
        self.mplusoneintratiowithmassplustwo = 0
        self.mplusonenormintratiowithmassplustwo = 0
        self.mplustwointratiowithmassplustwo = 0
        self.mplustwonormintratiowithmassplustwo = 0
        self.i_max = 0
        self.i_min = 0
        self.iso_max = 0
        self.iso_min = 0

        self.x_picked = 0
        self.y_picked = 0

        # canvas
        self.spectra_canvas = None
        self.con_canvas = None
        self.kin_canvas = None

        # plots
        self._spectra_ax = None
        self._con_ax = None
        self._kin_ax = None

        # ROI
        self.h = None
        self.mask = None
        self.ROI_listselect_text = ""
        self.ROI_listselect_array = []

        # multipmap compare variables
        self.ConcMapData = 0  # map data, x_end, y_end, ()
        self.IsotopeMapData = 0  # # map data, x_end, y_end, (iso_min, iso_max)

    # --- Executes on button press in find_IDlist.
    # can load a new ID list the consists of two columns
    # must be .csv file
    # the 1st column must be named "m/z"
    # the 2nd column must be named "Lipid ID"
    def find_IDlist_Callback(self):
        # But this code does not handle .h5 or .mat files
        self.fName_IDlist = QFileDialog.getOpenFileName(self, 'Pick ID List', filter='*.csv')
        self.IDlist_name.setText(self.fName_IDlist[0])
        self.ids_pd = pd.read_csv(self.fName_IDlist[0])

    # --- Executes on button press in mass_up.
    # When the button is clicked, the mass index is incremented by 1
    # and the image for the new mass is displayed
    def mass_up_Callback(self):
        if isIM:
            val = float(self.start.text()) + 1.0
            self.start.setText(str(val))
            self.im_point()
            return 0

        if self.has_data:
            self.index = self.index + 1
            self.refresh_image()
            # Update the mass and index boxes
            self.start.setText(str(self.z[self.index]))
            self.msindex.setText(str(self.index))
            self.Noise_Output_Box.setText(str(self.img_std[self.index]))
            # Update spectra annotation
            x_mass_up_val = self.z[self.index]
            y_mass_up_val = self.img_mean[self.index]
            x_mass_up_type = '%s' % float('%.6g' % x_mass_up_val)
            y_mass_up_type = '%s' % float('%.6g' % y_mass_up_val)

            if (pd.notnull(self.spectra_df['max'][self.index])) and (self.spectra_df['max'][self.index] != 0):
                threshold = float(self.pick_IDthreshold.text())
                for i in range(len(self.ids_pd['m/z'])):
                    err = abs((x_mass_up_val - self.ids_pd['m/z'][i]) / self.ids_pd['m/z'][i]) * self.err_multp
                    # print(err)
                    if (err < threshold):
                        ID_in = self.ids_pd['Lipid ID'][i]
                        break
                    else:
                        ID_in = 'not defined'
                self.ID_Output_Box.setText(str(ID_in))
                # Annotate
                self.annotate_spectra_ID(x_mass_up_type, y_mass_up_type, ID_in)
            else:
                self.annotate_spectra(x_mass_up_type, y_mass_up_type)

    # --- Executes on button press in mass_down.
    # When the button is clicked, the mass index is decremented by 1
    # and the image for the new mass is displayed
    def mass_down_Callback(self):
        if isIM:
            val = float(self.start.text()) - 1.0
            self.start.setText(str(val))
            self.im_point()
            return 0
        if self.has_data:
            self.index = self.index - 1
            self.refresh_image()
            # Update the mass and index boxes
            self.start.setText(str(self.z[self.index]))
            self.msindex.setText(str(self.index))
            self.Noise_Output_Box.setText(str(self.img_std[self.index]))
            # Update spectra annotation
            x_mass_down_val = self.z[self.index]
            y_mass_down_val = self.img_mean[self.index]
            x_mass_down_type = '%s' % float('%.6g' % x_mass_down_val)
            y_mass_down_type = '%s' % float('%.6g' % y_mass_down_val)

            if (pd.notnull(self.spectra_df['max'][self.index])) and (self.spectra_df['max'][self.index] != 0):
                threshold = float(self.pick_IDthreshold.text())
                for i in range(len(self.ids_pd['m/z'])):
                    err = abs((x_mass_down_val - self.ids_pd['m/z'][i]) / self.ids_pd['m/z'][i]) * self.err_multp
                    # print(err)
                    if (err < threshold):
                        ID_in = self.ids_pd['Lipid ID'][i]
                        break
                    else:
                        ID_in = 'not defined'
                self.ID_Output_Box.setText(str(ID_in))
                # Annotate
                self.annotate_spectra_ID(x_mass_down_type, y_mass_down_type, ID_in)
            else:
                self.annotate_spectra(x_mass_down_type, y_mass_down_type)

    # This function reads a mass input by the user, and then finds
    # the mass nearest to the input mass in the z data array.
    # It displays that mass and its index, then displays the image
    # at that mass.
    def start_Callback(self):
        if self.has_data:
            set_mass = self.start.text()  # returns contents of start as a double
            isNum = True
            try:
                float(set_mass)
            except ValueError:
                isNum = False
            if isNum == False:
                self.index = 1
                self.msindex.setText(str(self.index))
                print("Wrong Input Mass")
            else:
                set_mass = float(set_mass)
                # finds the index of the closest value to the set_mass value
                self.index = np.where(abs(self.z - set_mass) == (abs(self.z - set_mass)).min())[0][0]
            self.start.setText(str(self.z[self.index]))
            self.msindex.setText(str(self.index))
            self.refresh_image()
            # Update the Noise boxes
            self.Noise_Output_Box.setText(str(self.img_std[self.index]))
            # Update spectra annotation
            x_start_val = self.z[self.index]
            y_start_val = self.img_mean[self.index]
            x_start_type = '%s' % float('%.6g' % x_start_val)
            y_start_type = '%s' % float('%.6g' % y_start_val)

            if (pd.notnull(self.spectra_df['max'][self.index])) and (self.spectra_df['max'][self.index] != 0):
                threshold = float(self.pick_IDthreshold.text())
                for i in range(len(self.ids_pd['m/z'])):
                    err = abs((x_start_val - self.ids_pd['m/z'][i]) / self.ids_pd['m/z'][i]) * self.err_multp
                    # print(err)
                    if (err < threshold):
                        ID_in = self.ids_pd['Lipid ID'][i]
                        break
                    else:
                        ID_in = 'not defined'
                self.ID_Output_Box.setText(str(ID_in))
                # Annotate
                self.annotate_spectra_ID(x_start_type, y_start_type, ID_in)
            else:
                self.annotate_spectra(x_start_type, y_start_type)

    # This function gets an index input by the user and updates to image to
    # correspond to that index.
    def msindex_Callback(self):
        if self.has_data:
            temp = self.msindex.text()
            isNum = True
            try:
                float(temp)
            except ValueError:
                isNum = False
            if isNum == False:
                self.index = 1
                self.msindex.setText(str(self.index))
            else:
                self.index = np.int32(temp)
            if self.index > len(self.z):
                self.index = len(self.z - 1)
                self.msindex.setText(str(self.index))
            self.start.setText(str(self.z[self.index]))
            self.refresh_image()
            # Update the Noise boxes
            self.Noise_Output_Box.setText(str(self.img_std[self.index]))
            # Update spectra annotation
            x_msindex_val = self.z[self.index]
            y_msindex_val = self.img_mean[self.index]
            x_msindex_type = '%s' % float('%.6g' % x_msindex_val)
            y_msindex_type = '%s' % float('%.6g' % y_msindex_val)

            if (pd.notnull(self.spectra_df['max'][self.index])) and (self.spectra_df['max'][self.index] != 0):
                threshold = float(self.pick_IDthreshold.text())
                for i in range(len(self.ids_pd['m/z'])):
                    err = abs((x_msindex_val - self.ids_pd['m/z'][i]) / self.ids_pd['m/z'][i]) * self.err_multp
                    # print(err)
                    if (err < threshold):
                        ID_in = self.ids_pd['Lipid ID'][i]
                        break
                    else:
                        ID_in = 'not defined'
                self.ID_Output_Box.setText(str(ID_in))
                # Annotate
                self.annotate_spectra_ID(x_msindex_type, y_msindex_type, ID_in)
            else:
                self.annotate_spectra(x_msindex_type, y_msindex_type)

    def temp_max_Callback(self):
        if self.has_data:
            temp = self.temp_max.text()
            isNum = True
            try:
                float(temp)
            except ValueError:
                isNum = False
            if isNum == False:
                self.temp_max.setText(str(self.zmax))
                print("Wrong Input Max")
            else:
                set_temp = int(temp)
                self.t_max = set_temp
            self.temp_max.setText(str(self.t_max))
            self.zmax.setValue(math.ceil(self.t_max))
            self.scale_image()

    # This function updates the image to the current index
    def refresh_image(self):
        # finds the matrix of areas for the specified index
        bkg = np.zeros((len(self.y), len(self.x)))
        summarea = np.zeros((len(self.y), len(self.x)))
        self.areas = np.zeros((len(self.y), len(self.x)))
        numberpoints = 12
        for i in range(len(self.y)):
            for k in range(len(self.x)):
                sumeach = np.sum(self.intens[i, k, np.arange(self.index - numberpoints, (self.index + 70 + 1))])
                summarea[i, k] = sumeach
                sumback = np.sum(
                    self.intens[i, k, np.arange(self.index + numberpoints, (self.index + 2 * numberpoints + 1))])
                bkg[i, k] = sumback
                self.areas[i, k] = (sumeach - sumback)
        # pull the maximum value from the image and put it in the boxes
        self.z_max, self.z_min = self.arrlims(self.areas)
        self.max_int.setText(str(self.z_max))
        self.temp_max.setText(str(self.z_max))
        # set up the slider
        self.zmax.setMinimum(0)
        self.zmax.setMaximum(math.ceil(self.z_max))
        self.zmax.setValue(math.ceil(self.z_max))
        self.zmax.setSingleStep(round(0.01 * self.z_max))
        self.ConcMapData = [np.flip(self.areas, axis=0), self.x_end, self.y_end]
        # plot the image
        if self.con_canvas:
            self.plot_con.removeWidget(self.con_canvas)

            # self.con_cbar.remove()
            # Remove the colorbar
        if self._con_ax:
            self._con_ax.cla()

            # self.con_cbar.remove()  # Could the error be that there isn't a position set beforehand?
            # I don't know why this error is ocurring
            # Could it be that the colorbar should be replaced, not removed?
            self.con_img = self._con_ax.imshow(np.flip(self.areas, axis=0), cmap='jet', aspect='auto',
                                               extent=[0, self.x_end, 0, self.y_end])
            self._con_ax.set_xlabel('x, mm')
            self._con_ax.set_ylabel('y, mm')
            self.con_cbar = plt.colorbar(self.con_img)
            self._con_ax.invert_yaxis()
            self._con_ax.set_aspect('equal')
            self.con_canvas.draw()
            self.exConcflag = True

    # this function updates the isotope image and all the isotope ratio boxes
    # when the pick_point button is selected
    def refresh_isotoperatio(self):
        def good_column(inputarr):  # not in MATLAB
            temp = np.zeros(np.shape(self.good)[0] * np.shape(self.good)[1])
            A = np.reshape(inputarr, -1, order='F')
            B = np.reshape(self.good, -1, order='F')
            for i in range(len(B)):
                if B[i] == 1:
                    temp[i] = A[i]
                else:
                    temp[i] = 0
            temp = temp[temp != 0]
            return temp

        # Controls the resolution of sliders
        # The slider will have the equal amount of pieces between the max and min
        # ex. resolution =1000 will split the slider in to 1000 pieces
        if self.includemassplustwo:
            # setting the isotope ratio boxes to include the mass plus two values
            # setting the m ratio boxes
            good_mintratiowithmassplustwo = good_column(
                self.mintratiowithmassplustwo)  # creates a column of the good ratios from areas
            mean_mintratiowithmassplustwo = np.mean(
                good_mintratiowithmassplustwo)  # mean of the column of good ratios from areas
            std_dev_mintratiowithmassplustwo = np.std(good_mintratiowithmassplustwo,
                                                      ddof=1)  # st dev of the column of good ratios
            numberpts = len(good_mintratiowithmassplustwo)
            self.Msumratio.setText(
                str('%s' % float('%.5g' % mean_mintratiowithmassplustwo)))  # set the isotoperatio box
            self.Msumstandard_error.setText(
                str('%s' % float('%.5g' % std_dev_mintratiowithmassplustwo)))  # set the standard error box
            self.numberpoints.setText(str(numberpts))  # sets the number of points averaged
            # setting the mplusone ratio boxes
            good_mplusoneintratiowithmassplustwo = good_column(
                self.mplusoneintratiowithmassplustwo)  # creates a column of the good ratios from areas
            mean_mplusoneintratiowithmassplustwo = np.mean(
                good_mplusoneintratiowithmassplustwo)  # mean of the column of good ratios from areas
            std_dev_mplusoneintratiowithmassplustwo = np.std(good_mplusoneintratiowithmassplustwo,
                                                             ddof=1)  # st dev of the column of good ratios
            self.Mplusonesumratio.setText(
                str('%s' % float('%.5g' % mean_mplusoneintratiowithmassplustwo)))  # set the isotoperatio box
            self.Mplusonesumstandard_error.setText(
                str('%s' % float('%.5g' % std_dev_mplusoneintratiowithmassplustwo)))  # set the standard error box
            # setting the mplustwo ratio boxes
            good_mplustwointratiowithmassplustwo = good_column(
                self.mplustwointratiowithmassplustwo)  # creates a column of the good ratios from areas
            mean_mplustwointratiowithmassplustwo = np.mean(
                good_mplustwointratiowithmassplustwo)  # mean of the column of good ratios from areas
            std_dev_mplustwointratiowithmassplustwo = np.std(
                good_mplustwointratiowithmassplustwo)  # st dev of the column of good ratios
            self.Mplustwosumratio.setText(
                str('%s' % float('%.5g' % mean_mplustwointratiowithmassplustwo)))  # set the isotoperatio box
            self.Mplustwosumstandard_error.setText(
                str('%s' % float('%.5g' % std_dev_mplustwointratiowithmassplustwo)))  # set the standard error box
            # pull the maximum value from the image and put it in the boxes
            self.i_max, self.i_min = self.arrlims(self.mplusonenormintratiowithmassplustwo)
            self.max_int_iso.setText(str('%s' % float('%.5g' % self.i_max)))
            self.max_iso.setText(str('%s' % float('%.5g' % self.i_max)))
            self.min_int_iso.setText(str('%s' % float('%.5g' % self.i_min)))
            self.min_iso.setText(str('%s' % float('%.5g' % self.i_min)))
            # set up the two sliders
            self.zmax_truearr = np.linspace(self.i_min, self.i_max, self.resolution)
            self.zmax_isotope.setMinimum(0)
            self.zmax_isotope.setMaximum(self.resolution - 1)
            self.zmax_isotope.setValue(self.resolution - 1)
            self.zmax_isotope.setSingleStep(1)
            self.zmin_isotope.setMinimum(0)
            self.zmin_isotope.setMaximum(self.resolution - 1)
            self.zmin_isotope.setValue(0)
            self.zmin_isotope.setSingleStep(1)
            # isotope map storing variable
            self.IsotopeMapData = [np.flip(self.mplusonenormintratiowithmassplustwo, axis=0), self.x_end, self.y_end]
            # now tell the isotopeaxis to plot the image
            if self.kin_canvas:
                self._kin_ax.cla()
                self.kin_cbar.remove()
                self.kin_img = self._kin_ax.imshow(np.flip(self.mplusonenormintratiowithmassplustwo, axis=0),
                                                   cmap='jet', aspect='auto', extent=[0, self.x_end, 0, self.y_end])
                self._kin_ax.set_xlabel('x, mm')
                self._kin_ax.set_ylabel('y, mm')
                self.kin_cbar = self.kin_canvas.figure.colorbar(self.kin_img)
                self._kin_ax.invert_yaxis()
                self._kin_ax.set_aspect('equal')
                self.kin_canvas.draw()
                self.exIsotopeflag = True
            else:
                self.kin_canvas = FigureCanvas(plt.figure(tight_layout=True))
                self.kin_toolbar = NavigationToolbar(self.kin_canvas, self)
                self.plot_kin.addWidget(self.kin_toolbar)
                self.plot_kin.addWidget(self.kin_canvas)
                self._kin_ax = self.kin_canvas.figure.subplots()
                self.kin_img = self._kin_ax.imshow(np.flip(self.mplusonenormintratiowithmassplustwo, axis=0),
                                                   cmap='jet', aspect='auto', extent=[0, self.x_end, 0, self.y_end])
                self._kin_ax.set_xlabel('x, mm')
                self._kin_ax.set_ylabel('y, mm')
                self.kin_cbar = self.kin_canvas.figure.colorbar(self.kin_img)
                self._kin_ax.invert_yaxis()
                self._kin_ax.set_aspect('equal')
                self.exIsotopeflag = True
        else:
            # setting the isotope ratio boxes to include the mass plus two values
            # setting the m ratio boxes
            # good_mintratio = self.mintratio[self.good] # creates a column of the good ratios from areas
            good_mintratio = good_column(self.mintratio)  # creates a column of the good ratios from areas
            mean_mintratio = np.mean(good_mintratio)  # mean of the column of good ratios from areas
            std_dev_mintratio = np.std(good_mintratio, ddof=1)  # st dev of the column of good ratios
            numberpts = len(good_mintratio)
            self.Msumratio.setText(str('%s' % float('%.5g' % mean_mintratio)))  # set the isotoperatio box
            self.Msumstandard_error.setText(str('%s' % float('%.5g' % std_dev_mintratio)))  # set the standard error box
            self.numberpoints.setText(str(numberpts))  # sets the number of points averaged
            # setting the mplusone ratio boxes
            good_mplusoneintratio = good_column(self.mplusoneintratio)  # creates a column of the good ratios from areas
            mean_mplusoneintratio = np.mean(good_mplusoneintratio)  # mean of the column of good ratios from areas
            std_dev_mplusoneintratio = np.std(good_mplusoneintratio, ddof=1)  # st dev of the column of good ratios
            self.Mplusonesumratio.setText(str('%s' % float('%.5g' % mean_mplusoneintratio)))  # set the isotoperatio box
            self.Mplusonesumstandard_error.setText(
                str('%s' % float('%.5g' % std_dev_mplusoneintratio)))  # set the standard error box
            # setting the mplustwo ratio boxes
            self.Mplustwosumratio.setText("N/A")  # set the isotoperatio box
            self.Mplustwosumstandard_error.setText("N/A")  # set the standard error box
            # pull the maximum value from the image and put it in the boxes
            self.i_max, self.i_min = self.arrlims(self.mplusonenormintratio)
            self.max_int_iso.setText(str('%s' % float('%.5g' % self.i_max)))
            self.max_iso.setText(str('%s' % float('%.5g' % self.i_max)))
            self.min_int_iso.setText(str('%s' % float('%.5g' % self.i_min)))
            self.min_iso.setText(str('%s' % float('%.5g' % self.i_min)))
            # set up the two sliders
            self.zmax_truearr = np.linspace(self.i_min, self.i_max, self.resolution)
            self.zmax_isotope.setMinimum(0)
            self.zmax_isotope.setMaximum(self.resolution - 1)
            self.zmax_isotope.setValue(self.resolution - 1)
            self.zmax_isotope.setSingleStep(1)
            self.zmin_isotope.setMinimum(0)
            self.zmin_isotope.setMaximum(self.resolution - 1)
            self.zmin_isotope.setValue(0)
            self.zmin_isotope.setSingleStep(1)
            self.IsotopeMapData = [np.flip(self.mplusonenormintratio, axis=0), self.x_end, self.y_end]
            if self.kin_canvas:
                self._kin_ax.cla()
                self.kin_cbar.remove()
                self.kin_img = self._kin_ax.imshow(np.flip(self.mplusonenormintratio, axis=0), cmap='jet',
                                                   aspect='auto',
                                                   extent=[0, self.x_end, 0, self.y_end])
                self._kin_ax.set_xlabel('x, mm')
                self._kin_ax.set_ylabel('y, mm')
                self.kin_cbar = self.kin_canvas.figure.colorbar(self.kin_img)
                self._kin_ax.invert_yaxis()
                self._kin_ax.set_aspect('equal')
                self.kin_canvas.draw()
                self.exIsotopeflag = True
            else:
                self.kin_canvas = FigureCanvas(plt.figure(tight_layout=True))
                self.kin_toolbar = NavigationToolbar(self.kin_canvas, self)
                self.plot_kin.addWidget(self.kin_toolbar)
                self.plot_kin.addWidget(self.kin_canvas)
                self._kin_ax = self.kin_canvas.figure.subplots()
                self.kin_img = self._kin_ax.imshow(np.flip(self.mplusonenormintratio, axis=0), cmap='jet',
                                                   aspect='auto', extent=[0, self.x_end, 0, self.y_end])
                self._kin_ax.set_xlabel('x, mm')
                self._kin_ax.set_ylabel('y, mm')
                self.kin_cbar = self.kin_canvas.figure.colorbar(self.kin_img)
                self._kin_ax.invert_yaxis()
                self._kin_ax.set_aspect('equal')
                self.exIsotopeflag = True

    def arrlims(self, input_array):
        # this function locates the maximum and minimum values in a 2-d array
        minval = np.min(input_array)
        arrmin = np.min(minval)
        maxval = np.max(input_array)
        arrmax = np.max(maxval)
        return arrmax, arrmin

    # This function allows the user to select a region of interest,
    # then plots the mass spectrum averaged over that ROI
    # --- Executes on button press in ROI_select.
    def ROI_select_Callback_mask(self):
        if self.massplusone.isChecked():
            self.includemassplustwo = False
        elif self.massplustwo.isChecked():
            self.includemassplustwo = True
        if self.has_data:
            if self.h:
                self.h.disconnect()

            self.h = roi.new_ROI(self.con_img)

    # --- Executes on button press in ROI_process.
    def ROI_select_Callback_process(self):
        if self.has_data:
            if isIM:
                self.ROI_select_IM_Callback()
                return 0
            self.binI = self.h.get_mask().astype(int)
            # get average spectrum
            # this finds all the pixels (i.e. rows) that are part of the ROI

            # Below is a random ROI for testing purposes, keep it commented
            # self.binI = np.zeros((48, 60))
            # self.binI[9:13, 4] = 1
            # self.binI[8:14, 5] = 1
            # self.binI[7:15, 6] = 1
            # self.binI[7:16, 7] = 1
            # self.binI[7:15, 8] = 1
            # self.binI[7:15, 8] = 1
            # self.binI[7:16, 9] = 1
            # self.binI[6:17, 10:17] = 1
            # self.binI[5:17, 17:25] = 1
            # self.binI[5:17, 25:30] = 1
            # self.binI[5:16, 30:34] = 1

            self.f = np.argwhere(np.ravel(self.binI, order='F'))[:, 0]

            if self.massbox.text() != 0 and self.massbox.text() != '':
                massvalue = float(self.massbox.text())
                massindex = np.where(abs(self.z - (massvalue)) == (abs(self.z - (massvalue))).min())[0][0]
                massplusoneindex = np.where(abs(self.z - (massvalue + 1)) == (abs(self.z - (massvalue + 1))).min())[0][
                    0]
                massplustwoindex = np.where(abs(self.z - (massvalue + 2)) == (abs(self.z - (massvalue + 2))).min())[0][
                    0]
                # finds the index of the closest value to the massplusone, massplustwo and mass value
                mplusone = self.intens[:, :, massplusoneindex]  # this gets the massplusone plane
                mplustwo = self.intens[:, :, massplustwoindex]  # this gets the massplustwo plane
                m = self.intens[:, :, massindex]  # this gets mass plane
                zmax = np.max(m[:])  # finds the max overall in mass selected plane
                threshold = zmax / 7  # sets the threshold for the peaks to be real above a fifth of maximum
                good = (m > threshold).astype(
                    int)  # creates a matrix with ones and zeros for if the statement is true or not for each value

                self.mintratio = np.zeros((len(self.y), len(self.x)))  # creates matrix for integrated ratios
                self.mintratiowithmassplustwo = np.zeros(
                    (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                self.mplusoneintratio = np.zeros((len(self.y), len(self.x)))  # creates matrix for integrated ratios
                self.mplusoneintratiowithmassplustwo = np.zeros(
                    (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                self.mplustwointratiowithmassplustwo = np.zeros(
                    (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                bkg = np.zeros((len(self.y), len(self.x)))  # creates matrix for handles.Background
                summarr = np.zeros((len(self.y), len(self.x)))  # creates matrix for the area sums
                summplusonearr = np.zeros((len(self.y), len(self.x)))  # creates matrix for the sumplusone area sums
                summplustwoarr = np.zeros((len(self.y), len(self.x)))  # creates matrix for the sumplustwo area sums
                numberpoints = 12

                for i in range(len(self.y)):
                    for k in range(len(self.x)):
                        if self.good[i, k]:
                            summ = np.sum(self.intens[i, k, np.arange(self.index - numberpoints,
                                                                      (self.index + numberpoints + 1))])
                            summarr[i, k] = summ
                            summplusone = np.sum(self.intens[i, k, np.arange(massplusoneindex - numberpoints,
                                                                             (massplusoneindex + numberpoints + 1))])
                            summplusonearr[i, k] = summplusone
                            summplustwo = np.sum(self.intens[i, k, np.arange(massplustwoindex - numberpoints,
                                                                             (massplustwoindex + numberpoints + 1))])
                            summplustwoarr[i, k] = summplustwo
                            sumback = np.sum(self.intens[i, k, np.arange(self.index + numberpoints,
                                                                         (self.index + 2 * numberpoints + 1))])
                            bkg[i, k] = sumback
                            # getting all the integrated ratios for the main mass peak
                            self.mintratio[i, k] = (summ - sumback) / ((summ - sumback) + (summplusone - sumback))
                            self.mintratiowithmassplustwo[i, k] = (summ - sumback) / (
                                    (summ - sumback) + (summplusone - sumback) + (summplustwo - sumback))
                            # getting all the integrated ratios for the mass peak plus one
                            self.mplusoneintratio[i, k] = (summplusone - sumback) / (
                                    (summ - sumback) + (summplusone - sumback))
                            self.mplusoneintratiowithmassplustwo[i, k] = (summplusone - sumback) / (
                                    (summ - sumback) + (summplusone - sumback) + (summplustwo - sumback))
                            # getting all the integrated ratios for the mass peak plus two
                            self.mplustwointratiowithmassplustwo[i, k] = (summplustwo - sumback) / (
                                    (summ - sumback) + (summplusone - sumback) + (summplustwo - sumback))
                        else:
                            self.mintratio[i, k] = 0
                            self.mintratiowithmassplustwo[i, k] = 0
                            self.mplusoneintratio[i, k] = 0
                            self.mplusoneintratiowithmassplustwo[i, k] = 0
                            self.mplustwointratiowithmassplustwo[i, k] = 0
                # checked up to this point that everything is correct
                if self.includemassplustwo:
                    ROI_mintratiowithmassplustwo = np.zeros(
                        (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    ROI_mplusoneintratiowithmassplustwo = np.zeros(
                        (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    ROI_mplustwointratiowithmassplustwo = np.zeros(
                        (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    for i in range(len(self.y)):
                        for k in range(len(self.x)):
                            if self.binI[i, k]:
                                ROI_mintratiowithmassplustwo[i, k] = self.mintratiowithmassplustwo[i, k]
                                ROI_mplusoneintratiowithmassplustwo[i, k] = self.mplusoneintratiowithmassplustwo[i, k]
                                ROI_mplustwointratiowithmassplustwo[i, k] = self.mplustwointratiowithmassplustwo[i, k]
                            else:
                                ROI_mintratiowithmassplustwo[i, k] = 0
                                ROI_mplusoneintratiowithmassplustwo[i, k] = 0
                                ROI_mplustwointratiowithmassplustwo[i, k] = 0
                    ROImnonzeroswithmassplustwo = ROI_mintratiowithmassplustwo[np.nonzero(ROI_mintratiowithmassplustwo)]
                    ROImplusonenonzeroswithmassplustwo = ROI_mplusoneintratiowithmassplustwo[
                        np.nonzero(ROI_mplusoneintratiowithmassplustwo)]
                    ROImplustwononzeroswithmassplustwo = ROI_mplustwointratiowithmassplustwo[
                        np.nonzero(ROI_mplustwointratiowithmassplustwo)]
                    ROImaveragewithmassplustwo = np.mean(ROImnonzeroswithmassplustwo)
                    ROImplusoneaveragewithmassplustwo = np.mean(ROImplusonenonzeroswithmassplustwo)
                    ROImplustwoaveragewithmassplustwo = np.mean(ROImplustwononzeroswithmassplustwo)
                    ROImstdevwithmassplustwo = np.std(ROImnonzeroswithmassplustwo, ddof=1)
                    ROImplusonestdevwithmassplustwo = np.std(ROImplusonenonzeroswithmassplustwo, ddof=1)
                    ROImplustwostdevwithmassplustwo = np.std(ROImplustwononzeroswithmassplustwo, ddof=1)
                    NumberinROI = len(ROImnonzeroswithmassplustwo)
                    self.Msumratio.setText(str('%s' % float('%.5g' % ROImaveragewithmassplustwo)))  # set isotope ratio
                    self.Msumstandard_error.setText(
                        str('%s' % float('%.5g' % ROImstdevwithmassplustwo)))  # set error box
                    self.numberpoints.setText(str(NumberinROI))  # set number of points
                    self.Mplusonesumratio.setText(
                        str('%s' % float('%.5g' % ROImplusoneaveragewithmassplustwo)))  # set isotope ratio
                    self.Mplusonesumstandard_error.setText(
                        str('%s' % float('%.5g' % ROImplusonestdevwithmassplustwo)))  # set error box
                    self.Mplustwosumratio.setText(
                        str('%s' % float('%.5g' % ROImplustwoaveragewithmassplustwo)))  # set isotope ratio
                    self.Mplustwosumstandard_error.setText(
                        str('%s' % float('%.5g' % ROImplustwostdevwithmassplustwo)))  # set error box

                else:
                    ROI_mintratio = np.zeros((len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    ROI_mplusoneintratio = np.zeros((len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    for i in range(len(self.y)):
                        for k in range(len(self.x)):
                            if self.binI[i, k]:
                                ROI_mintratio[i, k] = self.mintratio[i, k]
                                ROI_mplusoneintratio[i, k] = self.mplusoneintratio[i, k]
                            else:
                                ROI_mintratio[i, k] = 0
                                ROI_mplusoneintratio[i, k] = 0
                    ROImnonzeros = ROI_mintratio[np.nonzero(ROI_mintratio)]
                    ROImplusonenonzeros = ROI_mplusoneintratio[np.nonzero(ROI_mplusoneintratio)]
                    ROImaverage = np.mean(ROImnonzeros)
                    ROImplusoneaverage = np.mean(ROImplusonenonzeros)
                    ROImstdev = np.std(ROImnonzeros, ddof=1)
                    ROImplusonestdev = np.std(ROImplusonenonzeros, ddof=1)
                    NumberinROI = len(ROImnonzeros)
                    self.Msumratio.setText(str('%s' % float('%.5g' % ROImaverage)))  # set the isotoperatio box
                    self.Msumstandard_error.setText(str('%s' % float('%.5g' % ROImstdev)))  # set error box
                    self.numberpoints.setText(str(NumberinROI))  # set number of points
                    self.Mplusonesumratio.setText(str('%s' % float('%.5g' % ROImplusoneaverage)))  # set isotope ratio
                    self.Mplusonesumstandard_error.setText(
                        str('%s' % float('%.5g' % ROImplusonestdev)))  # set error box
                    self.Mplustwosumratio.setText("N/A")  # set isotope ratio
                    self.Mplustwosumstandard_error.setText("N/A")  # set error box
            else:
                print("No data loaded")
            self.img_mean = np.mean(self.reshaped[self.f, :], axis=0)
            self.img_std = self.Noise
            self.spectra_df = pd.DataFrame({'m/z': self.z, 'intensity': self.img_mean})
            self.spectra_df['max'] = self.spectra_df.iloc[
                signal.argrelextrema(self.spectra_df.intensity.values, np.greater_equal, order=self.checkpoint_maxima)[
                    0]]['intensity']
            # Plot the spectrum averaged over the ROI
            if self.spectra_canvas:
                self._spectra_ax.cla()
                self._spectra_ax.plot(self.z, self.img_mean, 'k', linewidth=0.3, picker=True)
                self._spectra_ax.set_title('Average Spectrum Across Selected Region')
                self._spectra_ax.set_xlabel('m/z')
                self._spectra_ax.set_ylabel('intensity')
                self.rgflag = True  # this is so that a MS spectrum can be exported
                self.rgflagROI = True  # this is so that the ROI can be saved to listbox
            else:
                print("No data loaded")
        else:
            print("No data loaded")

    def ROI_select_IM_Callback(self):
        self.binI = self.h.get_mask().astype(int)
        self.binI = np.flipud(self.binI)
        f = np.argwhere(np.ravel(self.binI, order='F'))[:, 0]

        x = self.chosenData

        # raveled = np.ravel(self.chosenData, order='F')

        theList = []

        i = 0

        for line in x:
            for frame in line:
                if i in f and frame != 0:
                    theList.append(frame)
                i += 1
            i += 1

        filtered = []
        for val in theList:
            for val2 in val:
                filtered.append(val2)

        self.ROIData = filtered

        # with open("ROI.csv", 'w', newline='') as file:
        #     writer = csv.writer(file)
        #     writer.writerow(["M/Z Value", "Intensity", "Drift Time", "Line", "Frame Num"])
        #
        #     for line in filtered:
        #         writer.writerow(line)
        # print("Done writing to .csv")
        # Need to find the m/z and not just the intensity at every single point. how to do this?
        # Have the whole m/z spectra and then be able to select just 1 m/z

    # This function plots a mass spectrum corresponing to a selected point
    # on the displayed image, or an image corresponding to a selected point
    # on the mass spectrum
    # --- Executes on button press in pick_point.
    def pick_point_Callback(self):  # this function only considers "line" referring to MATLAB. Not sure what it means
        if isIM:
            self.im_point()
            # self.refresh_isotoperatio()
            return 0

        # choose what masspluswhat peak
        if self.massplusone.isChecked():
            self.includemassplustwo = False
        elif self.massplustwo.isChecked():
            self.includemassplustwo = True
        self.start.setText(str(self.z[self.index]))
        self.msindex.setText(str(self.index))
        self.Noise_Output_Box.setText(str(self.img_std[self.index]))
        self.refresh_image()
        # calculating the area and plotting the isotope image
        m = self.intens[:, :, self.index]  # this index gets the selected plane of values for this index
        # now find all the M+1 and M+2 planes
        massplusonemass = self.z[self.index] + 1
        massplustwomass = self.z[self.index] + 2
        massplusoneindex = np.where(abs(self.z - massplusonemass) == (abs(self.z - massplusonemass)).min())[0][0]
        massplustwoindex = np.where(abs(self.z - massplustwomass) == (abs(self.z - massplustwomass)).min())[0][0]
        # finds the index of the closest value to the massplusonemass and massplustwomass value
        mplusone = self.intens[:, :, massplusoneindex]  # this gets the massplusone plane
        mplustwo = self.intens[:, :, massplustwoindex]  # this gets the massplustwo plane
        self.massbox.setText(str(self.z[self.index]))  # sets the mass box
        zmax = np.max(m)  # finds the max overall in mass selected plane
        threshold = zmax / 7  # sets the threshold for the peaks to be real above a fifth of maximum
        self.good = (m > threshold).astype(int)  # creates a matrix with ones and zeros for if the
        # statement is true or not for each value
        normfactor = int(self.ideal_ratio.text())  # gets the normalization factor from ideal_ratio
        self.mintratio = np.zeros((len(self.y), len(self.x)))  # creates matrix for integrated ratios
        self.mnormintratio = np.ones((len(self.y), len(self.x)))  # creates matrix for normalized ratios
        self.mintratiowithmassplustwo = np.zeros((len(self.y), len(self.x)))  # creates matrix for integrated ratios
        self.mnormintratiowithmassplustwo = np.ones((len(self.y), len(self.x)))  # creates matrix for normalized ratios
        self.mplusoneintratio = np.zeros((len(self.y), len(self.x)))  # creates matrix for integrated ratios
        self.mplusonenormintratio = np.ones((len(self.y), len(self.x)))  # creates matrix for normalized ratios
        self.mplusoneintratiowithmassplustwo = np.zeros(
            (len(self.y), len(self.x)))  # creates matrix for integrated ratios
        self.mplusonenormintratiowithmassplustwo = np.ones(
            (len(self.y), len(self.x)))  # creates matrix for normalized ratios
        self.mplustwointratiowithmassplustwo = np.zeros(
            (len(self.y), len(self.x)))  # creates matrix for integrated ratios
        self.mplustwonormintratiowithmassplustwo = np.ones(
            (len(self.y), len(self.x)))  # creates matrix for normalized ratios
        bkg = np.zeros((len(self.y), len(self.x)))  # creates matrix for handles.Background
        summarr = np.zeros((len(self.y), len(self.x)))  # creates matrix for the area sums
        summplusonearr = np.zeros((len(self.y), len(self.x)))  # creates matrix for the sumplusone area sums
        summplustwoarr = np.zeros((len(self.y), len(self.x)))  # creates matrix for the sumplustwo area sums
        numberpoints = 12
        for i in range(len(self.y)):
            for k in range(len(self.x)):
                if self.good[i, k]:
                    summ = np.sum(
                        self.intens[i, k, np.arange(self.index - numberpoints, (self.index + numberpoints + 1))])
                    summarr[i, k] = summ
                    summplusone = np.sum(self.intens[i, k, np.arange(massplusoneindex - numberpoints,
                                                                     (massplusoneindex + numberpoints + 1))])
                    summplusonearr[i, k] = summplusone
                    summplustwo = np.sum(self.intens[i, k, np.arange(massplustwoindex - numberpoints,
                                                                     (massplustwoindex + numberpoints + 1))])
                    summplustwoarr[i, k] = summplustwo
                    sumback = np.sum(
                        self.intens[i, k, np.arange(self.index + numberpoints, (self.index + 2 * numberpoints + 1))])
                    bkg[i, k] = sumback
                    # getting all the integrated ratios for the main mass peak
                    self.mintratio[i, k] = (summ - sumback) / ((summ - sumback) + (summplusone - sumback))
                    self.mintratiowithmassplustwo[i, k] = (summ - sumback) / (
                            (summ - sumback) + (summplusone - sumback) + (summplustwo - sumback))
                    self.mnormintratio[i, k] = self.mintratio[i, k] / normfactor
                    self.mnormintratiowithmassplustwo[i, k] = self.mintratiowithmassplustwo[i, k] / normfactor
                    # getting all the integrated ratios for the mass peak plus one
                    self.mplusoneintratio[i, k] = (summplusone - sumback) / ((summ - sumback) + (summplusone - sumback))
                    self.mplusoneintratiowithmassplustwo[i, k] = (summplusone - sumback) / (
                            (summ - sumback) + (summplusone - sumback) + (summplustwo - sumback))
                    self.mplusonenormintratio[i, k] = self.mplusoneintratio[i, k] / normfactor
                    self.mplusonenormintratiowithmassplustwo[i, k] = self.mplusoneintratiowithmassplustwo[
                                                                         i, k] / normfactor
                    # getting all the integrated ratios for the mass peak plus two
                    self.mplustwointratiowithmassplustwo[i, k] = (summplustwo - sumback) / (
                            (summ - sumback) + (summplusone - sumback) + (summplustwo - sumback))
                    self.mplustwonormintratiowithmassplustwo[i, k] = self.mplustwointratiowithmassplustwo[
                                                                         i, k] / normfactor
                else:
                    self.mintratio[i, k] = 0
                    self.mintratiowithmassplustwo[i, k] = 0
                    self.mnormintratio[i, k] = 1
                    self.mnormintratiowithmassplustwo[i, k] = 1
                    self.mplusoneintratio[i, k] = 0
                    self.mplusoneintratiowithmassplustwo[i, k] = 0
                    self.mplusonenormintratio[i, k] = 1
                    self.mplusonenormintratiowithmassplustwo[i, k] = 1
                    self.mplustwointratiowithmassplustwo[i, k] = 0
                    self.mplustwonormintratiowithmassplustwo[i, k] = 1
        self.refresh_isotoperatio()

    def im_point(self):
        # NOTE: to Brian . Okay, here's the deal, when building this the code was so confusing,
        # so we are going to do everything here. In reality, we shouldn't do this, and we won't, in the long run.
        # When cleaning up, this must be factored into several different functions.

        # A couple more notes:
        # 1. How far apart should we accept? Surely not the only val that was selected? +- .5 m/z??
        # 2. You probably shouldn't do all of these calculations over again,
        # probably should pull them out of chosen data.

        fileID = open(self.cubefilename)
        data = np.fromfile(fileID, dtype=np.float32)
        frameNum = data[0]
        fileNum = data[1]

        numFrames = 0
        numFiles = 0
        i = 2

        chosenData = []
        theChosenData = []
        chosenDataPlusOne = []
        theChosenDataPlusOne = []
        chosenDataPlusTwo = []
        theChosenDataPlusTwo = []
        frameDone = False

        maxIntensity = 0
        maxIntensityPlusOne = 0
        maxIntensityPlusTwo = 0

        while numFiles < fileNum:
            if frameDone:
                chosenData.append(lineData)
                theChosenData.append(otherLine)
                chosenDataPlusOne.append(lineDataPlusOne)
                theChosenDataPlusOne.append(otherLinePlusOne)
                chosenDataPlusTwo.append(lineDataPlusTwo)
                theChosenDataPlusTwo.append(otherLinePlusTwo)
                numFiles += 1
                numFrames = 0
                frameDone = False
                continue
            frameDone = False
            lineData = []
            lineDataPlusOne = []
            lineDataPlusTwo = []
            otherLine = []
            valAdded = False
            valAddedPlusOne = False
            valAddedPlusTwo = False
            theVal = 0
            theValPlusOne = 0
            theValPlusTwo = 0
            otherLine = []
            otherLinePlusOne = []
            otherLinePlusTwo = []
            while numFrames < frameNum:
                totalDriftBins = data[i]
                frameDone = True
                currdriftBin = 0
                i += 1
                while currdriftBin < totalDriftBins:
                    numValues = data[i]
                    driftTime = data[i + 1]
                    i += 2
                    for a in range(int(numValues)):
                        if (data[i] >= (float(self.start.text()) - .5)) and (
                                data[i] < (float(self.start.text()) + .5)):
                            theVal += data[i + 1]
                            otherVal.append([data[i], data[i + 1], driftTime, numFiles, numFrames])
                            valAdded = True
                            if data[i + 1] > maxIntensity:
                                maxIntensity = data[i + 1]
                        if (data[i] >= (float(self.start.text()) + .5)) and (
                                data[i] < (float(self.start.text()) + 1.5)):
                            theValPlusOne += data[i + 1]
                            otherValPlusOne.append([data[i], data[i + 1], driftTime, numFiles, numFrames])
                            valAddedPlusOne = True
                            if data[i + 1] > maxIntensityPlusOne:
                                maxIntensityPlusOne = data[i + 1]
                        if (data[i] >= (float(self.start.text()) + 1.5)) and (
                                data[i] < (float(self.start.text()) + 2.5)):
                            theValPlusTwo += data[i + 1]
                            otherValPlusTwo.append([data[i], data[i + 1], driftTime, numFiles, numFrames])
                            valAddedPlusTwo = True
                            if data[i + 1] > maxIntensityPlusTwo:
                                maxIntensityPlusTwo = data[i + 1]
                        i += 2
                    currdriftBin += 1

                if not valAdded:
                    lineData.append(0)
                    otherLine.append(0)
                elif valAdded:
                    lineData.append(theVal)
                    otherLine.append(otherVal)
                if not valAddedPlusOne:
                    lineDataPlusOne.append(0)
                    otherLinePlusOne.append(0)
                elif valAddedPlusOne:
                    lineDataPlusOne.append(theValPlusOne)
                    otherLinePlusOne.append(otherValPlusOne)
                if not valAddedPlusTwo:
                    lineDataPlusTwo.append(0)
                    otherLinePlusTwo.append(0)
                elif valAddedPlusTwo:
                    lineDataPlusTwo.append(theValPlusTwo)
                    otherLinePlusTwo.append(otherValPlusTwo)
                valAdded = False
                valAddedPlusOne = False
                valAddedPlusTwo = False
                otherVal = []
                otherValPlusOne = []
                otherValPlusTwo = []
                theVal = 0
                theValPlusOne = 0
                theValPlusTwo = 0
                numFrames += 1

        numY = len(chosenData)
        numX = len(chosenData[0])
        xend = numX * .075
        yend = numY * .15

        self.max_int.setText(str(maxIntensity))
        self.temp_max.setText(str(0))
        self.zmax.setMinimum(0)
        self.zmax.setMaximum(int(maxIntensity))

        if self.view:
            self.plot_con.removeWidget(self.view)

        self.view = FigureCanvas(Figure(figsize=(5, 3)))
        self.axes = self.view.figure.subplots()
        self.toolbar = NavigationToolbar(self.view, self)
        self.plot_con.addWidget(self.view)
        self.con_img = self.axes.imshow(chosenData, cmap='jet',
                                        aspect=(yend / xend), extent=[0, xend, 0, yend])
        plt.colorbar(self.con_img)
        self.view.draw()

        if self.massplusone.isChecked():
            if self.viewPlusOne:
                self.plot_kin.removeWidget(self.viewPlusOne)
            if self.viewPlusTwo:
                self.plot_kin.removeWidget(self.viewPlusTwo)

            self.viewPlusOne = FigureCanvas(Figure(figsize=(5, 3)))
            self.axes = self.viewPlusOne.figure.subplots()
            self.toolbar = NavigationToolbar(self.viewPlusOne, self)
            self.plot_kin.addWidget(self.viewPlusOne)
            self.con_img2 = self.axes.imshow(chosenDataPlusOne, cmap='inferno',
                                             aspect=(yend / xend), extent=[0, xend, 0, yend])
            plt.colorbar(self.con_img2)
            self.viewPlusOne.draw()

            self.max_iso.setText(str(maxIntensityPlusOne))
            self.max_int_iso.setText(str(maxIntensityPlusOne))
            self.zmax_isotope.setMinimum(0)
            self.zmax_isotope.setMaximum(int(maxIntensityPlusOne))

            self.min_iso.setText(str(0))
            self.min_int_iso.setText(str(0))
            self.zmin_isotope.setMinimum(0)
            self.zmin_isotope.setMaximum(int(maxIntensityPlusOne))
            self.zmax_isotope.setValue(int(maxIntensityPlusOne))

        elif self.massplustwo.isChecked():
            if self.viewPlusTwo:
                self.plot_kin.removeWidget(self.viewPlusTwo)
            if self.viewPlusOne:
                self.plot_kin.removeWidget(self.viewPlusOne)

            self.viewPlusTwo = FigureCanvas(Figure(figsize=(5, 3)))
            self.axes = self.viewPlusTwo.figure.subplots()
            self.toolbar = NavigationToolbar(self.viewPlusTwo, self)
            self.plot_kin.addWidget(self.viewPlusTwo)
            self.con_img2 = self.axes.imshow(chosenDataPlusTwo, cmap='inferno', interpolation='gaussian',
                                             aspect=(yend / xend), extent=[0, xend, 0, yend])
            plt.colorbar(self.con_img2)
            self.viewPlusTwo.draw()

            self.max_iso.setText(str(maxIntensityPlusTwo))
            self.max_int_iso.setText(str(maxIntensityPlusTwo))
            self.zmax_isotope.setMinimum(0)
            self.zmax_isotope.setMaximum(int(maxIntensityPlusTwo))

            self.min_iso.setText(str(0))
            self.min_int_iso.setText(str(0))
            self.zmin_isotope.setMinimum(0)
            self.zmin_isotope.setMaximum(int(maxIntensityPlusTwo))
            self.zmax_isotope.setValue(int(maxIntensityPlusTwo))

        self.chosenData = theChosenData
        if self.massplusone:
            self.chosenDataIso = theChosenDataPlusOne
        else:
            self.chosenDataIso = theChosenDataPlusTwo

        self.zmax.setValue(0)
        self.zmin_isotope.setValue(0)

    # --- Executes on button press in find_file.
    def find_file_Callback(self):
        # But this code does not handle .h5 or .mat files
        self.fName = QFileDialog.getOpenFileName(self, 'Pick Data Cube', filter='*.mat, *.h5 *.bin')
        self.wspc_name.setText(self.fName[0])

        # --- Executes on button press in start_cube.

    def start_cube_Callback(self):
        # hObject    handle to start_cube (see GCBO)
        # eventdata  reserved - to be defined in a future version of MATLAB
        # handles    structure with handles and user data (see GUIDATA)
        # The following statements load the data cube into the data structure
        # shared by the functions that serve the GUI
        # These lines set determine which units were used for the x and y
        # coordinates. The default is mm. Other possible units include
        # micrometers, entered as either 'um' or 'microns', and cm.
        # self.IMDataButton

        self.ppm_min.setMinimum(1)
        self.ppm_min.setMaximum(10000)
        self.ppm_min.setValue(100)

        self.cubefilename = self.fName[0]
        filename = self.cubefilename
        print("Working to read datacube")
        if filename.endswith('.mat'):
            print("This code can't process .mat files")
            return
        elif filename.endswith('.h5'):
            print("This code can't process .h5 files")
            return
        elif filename.endswith('.bin'):
            print("File extension: .bin")

            if isIM == False:
                self.cubeAsMSData(filename)
            elif isIM:
                self.cubeAsIMData(filename)
            else:
                print("Please select whether the file is IM or MS Data")
                return

        else:
            print('Unexpected file extension')
            return

    def cubeAsIMData(self, filename):
        fileID = open(filename)
        data = np.fromfile(fileID, dtype=np.float32)
        frameNum = data[0]
        fileNum = data[1]

        numFrames = 0
        numFiles = 0
        i = 2

        mzVals = []
        intensity = []
        drifts = []
        chosenData = []
        frameDone = False

        while numFiles < fileNum:
            if frameDone:
                chosenData.append(lineData)
                numFiles += 1
                numFrames = 0
                frameDone = False
                continue
            frameDone = False
            lineData = []
            valAdded = False
            theVal = 0
            while numFrames < frameNum:
                totalDriftBins = data[i]
                frameDone = True
                currdriftBin = 0
                i += 1
                while currdriftBin < totalDriftBins:
                    numValues = data[i]
                    i += 1
                    driftTime = data[i]
                    i += 1
                    for a in range(int(numValues)):
                        theVal += data[i + 1]
                        valAdded = True
                        intensity.append(data[i + 1])
                        drifts.append(driftTime)
                        mzVals.append(data[i])
                        i += 2
                    currdriftBin += 1

                if not valAdded:
                    lineData.append(0)
                elif valAdded:
                    lineData.append(theVal)
                valAdded = False
                theVal = 0
                numFrames += 1

        if self._spectra_ax:
            self.plot_spectra.removeWidget(self.spectra_toolbar)
            self.plot_spectra.removeWidget(self.spectra_canvas)
        if self.viewPlusOne:
            self.plot_kin.removeWidget(self.viewPlusOne)
            self.plot_kin.addWidget(FigureCanvas(plt.figure(tight_layout=True)))
        elif self.viewPlusTwo:
            self.plot_kin.removeWidget(self.viewPlusTwo)
            self.plot_kin.addWidget(FigureCanvas(plt.figure(tight_layout=True)))

        self.spectra_canvas = FigureCanvas(plt.figure(tight_layout=True))
        self.spectra_canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.spectra_canvas.setFocus()
        self.spectra_toolbar = NavigationToolbar(self.spectra_canvas, self)
        self.plot_spectra.addWidget(self.spectra_toolbar)
        self.plot_spectra.addWidget(self.spectra_canvas)
        self._spectra_ax = self.spectra_canvas.figure.subplots()
        x = self._spectra_ax.scatter(mzVals, intensity, s=.01, c=drifts, cmap="Greens", alpha=0.75, picker=True)
        plt.colorbar(x).set_label('Drift times')
        self._spectra_ax.set_title('Points In Selected Region')
        self._spectra_ax.set_xlabel('m/z')
        self._spectra_ax.set_ylabel('intensity')
        self.spectra_canvas.mpl_connect('pick_event', self.data_cursor_click)
        self.spectra_canvas.mpl_connect('key_press_event', self.data_cursor_key)
        self.exSpecflag = True
        # plt.yscale('log')
        # it is x, y
        plt.ylabel('intensity')
        plt.xlabel('m/z')

        self.mzVals = mzVals
        self.intensity = intensity
        self.drifts = drifts

        numY = len(chosenData)
        numX = len(chosenData[0])
        xend = numX * .075
        yend = numY * .15

        if self.view:
            self.plot_con.removeWidget(self.view)
        self.view = FigureCanvas(Figure(figsize=(5, 3)))
        self.axes = self.view.figure.subplots()
        self.toolbar = NavigationToolbar(self.view, self)
        self.plot_con.addWidget(self.view)
        self.con_img = self.axes.imshow(chosenData, cmap='jet', interpolation='gaussian',
                                        aspect=(yend / xend), extent=[0, xend, 0, yend])
        plt.colorbar(self.con_img)
        # self.axes.set_title('Points In Selected Region')
        # self.axes.set_xlabel('m/z')
        # self.axes.set_ylabel('intensity')
        self.view.draw()

        self.has_data = 1
        # plt.show()
        # print("Process the IM Data")

    def cubeAsMSData(self, filename):
        fileID = open(filename)
        data = np.fromfile(fileID, dtype=np.float32)
        x = int(data[0])
        y = int(data[1])
        z = int(data[2])
        end = len(data)
        imgZ = data[end - z:end]
        img = data[3:x * y * z + 3]
        del data
        img = img.reshape(z, x, y, order='F')
        img = img.transpose(2, 1, 0)
        img = np.flip(img, 0)
        # These two lines of code set the pixel size to 75*150
        imgX = np.arange(37.5, 75 * x + 37.5, 75)
        imgY = np.arange(75, 150 * y + 75, 150)
        fileID.close()

        print("Working to display data\n")

        # if statement that checks if ROI with a certain value exists. Not sure how it accepts ROI values yet
        # below belongs to the "else" part
        ROIcount = 0
        if (self.ROI and self.ROIcount):
            self.ROIcount = ROIcount
            self.ROIcountbox.setText(str(self.ROIcount))
        else:
            ROIcount = 0
            self.ROIcount = ROIcount
            self.ROIcountbox.setText(str(self.ROIcount))

        self.has_data = 1
        if self.micrometer.isChecked():
            scalefact = 1e3
        elif self.millimeter.isChecked():
            scalefact = 1
        elif self.centimeter.isChecked():
            scalefact = 0.1
        # Note that x and y units are converted to mm, and the
        # origin for the plot is set to the corner of the image.
        # Note also that the x and y arrays for images have to be
        # double precision for the data cursor to work.
        self.x = abs((imgX - imgX[0]) / scalefact)
        xm, xn = 1, len(imgX)
        if xm > xn:
            self.x = self.x.conj().T
        self.x_end = self.x[len(self.x) - 1]
        self.y = abs((imgY - imgY[0]) / scalefact)
        ym, yn = 1, len(imgY)
        if ym < yn:
            self.y = self.y.conj().T
        self.y_end = self.y[len(self.y) - 1]
        self.z = imgZ
        zm, zn = len(self.z), 1
        if zm < zn:
            self.z = self.z.conj().T
        self.ptflag = False
        self.rgflag = False
        self.rgflagROI = False
        self.intens = img
        # this transforms the datacube into a 2D matrix,
        # where the rows are the pixels and columns are the m/z values
        self.a, self.b, self.c = np.shape(self.intens)[0:3]
        self.reshaped = np.reshape(self.intens, (self.a * self.b, self.c), order="F")
        # Plot the mass spectrum averaged over the entire image.
        self.img_mean = np.mean(self.reshaped, axis=0)

        # Now, we get the handles.Background handles.Noise by: First, finding the PS peak in the
        # data. Then, we check if a scan is handles.Background by comparing it to the
        # average PS peak (greater than 1/8 average PS intensity is brain, less than
        # 1/8 average PS intensity is handles.Background) We can then find the handles.Noise at each
        # m/z by taking the standard deviation of all the handles.Background intensities at
        # each m/z.
        self.PS_Peak_Intensity = 0
        self.PS_Sum_Intensity = 0
        imgYSize = int(len(imgY.astype(np.int32)) / 2)
        imgXSize = int(len(imgX.astype(np.int32)) / 2)
        for i in range(len(imgZ)):
            if ((imgZ[i] > 834.4) and (imgZ[i] < 834.7)):  # PS is usually 834.55 and the
                # Agilent should be accurate within a few hundreths of a m/z.
                if img[imgYSize - 1, imgXSize - 1, i] > self.PS_Peak_Intensity:
                    # I selected the middle point to test where the PS peak was in
                    # hopes that the middle point is on the brain. If it is not,
                    # this may cause us to filter by a contaminant peak not PS.
                    self.PS_Peak_Intensity = img[imgYSize - 1, imgXSize - 1, i]
                    self.PS_Index = i  # PS_Index should be smaller by 1 than that of in MATLAB
        # this adds up all PS intensities in the whole image then divides by the
        # number of pixels to get the average intensity.
        for i in range(len(imgY)):
            for j in range(len(imgX)):
                self.PS_Sum_Intensity = self.PS_Sum_Intensity + img[i, j, self.PS_Index]
        self.PS_Mean_Intensity = self.PS_Sum_Intensity / (len(imgY) * len(imgX))
        # Now, we will create an array of the handles.Background by comparing each point to
        # 1/8 of the handles.PS_Mean_Intensity. If it is a lower value, we will add it to
        # the handles.Background array.
        print("Calculating self.Noise\n")
        self.Background = np.zeros(((len(imgY) * len(imgX)), len(imgZ)))
        WhereIsBackground = np.zeros((len(imgY), len(imgX)))
        k = 0
        for i in range(len(imgY)):
            for j in range(len(imgX)):
                if img[i, j, self.PS_Index] < (self.PS_Mean_Intensity / 8):
                    self.Background[k, :] = img[i, j, :]
                    WhereIsBackground[i, j] = 1
                    k = k + 1
        self.Background = np.delete(self.Background, np.s_[k:(len(imgY) * len(imgX)) + 1], 0)
        # Now we calculate the handles.Noise as the standard deviation of all intensities at
        # an m/z in the handles.Noise. An arry of stdev's is created at each m/z

        self.Noise = copy(imgZ)
        for i in range(len(imgZ)):
            self.Noise[i] = np.std(self.Background[:, i], ddof=1)
            if (math.isnan(self.Noise[i])):
                self.Noise[i] = 0
        # Now we calculate the baseline signal at each m/z by averaging the
        # handles.Background intensities.
        self.Background_Signal = np.mean(self.Background, 0)
        self.spectra_df = pd.DataFrame({'m/z': self.z, 'intensity': self.img_mean})
        self.spectra_df['max'] = self.spectra_df.iloc[
            signal.argrelextrema(self.spectra_df.intensity.values, np.greater_equal, order=self.checkpoint_maxima)[0]][
            'intensity']
        self.img_std = self.Noise
        if self.spectra_canvas:
            self._spectra_ax.cla()
            self._spectra_ax.plot(self.z, self.img_mean, 'k', linewidth=0.3, picker=True)
            self._spectra_ax.set_title('Average Spectrum Across Selected Region')
            self._spectra_ax.set_xlabel('m/z')
            self._spectra_ax.set_ylabel('intensity')
        else:
            self.spectra_canvas = FigureCanvas(plt.figure(tight_layout=True))
            self.spectra_canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
            self.spectra_canvas.setFocus()
            self.spectra_toolbar = NavigationToolbar(self.spectra_canvas, self)
            self.plot_spectra.addWidget(self.spectra_toolbar)
            self.plot_spectra.addWidget(self.spectra_canvas)
            self._spectra_ax = self.spectra_canvas.figure.subplots()
            self._spectra_ax.plot(self.z, self.img_mean, 'k', linewidth=0.3, picker=True)
            self._spectra_ax.set_title('Average Spectrum Across Selected Region')
            self._spectra_ax.set_xlabel('m/z')
            self._spectra_ax.set_ylabel('intensity')
            self.spectra_canvas.mpl_connect('pick_event', self.data_cursor_click)
            self.spectra_canvas.mpl_connect('key_press_event', self.data_cursor_key)
            self.exSpecflag = True

        # Pick a point in the middle of the data set to start
        # This means that the initial image will typically be junk
        self.index = np.uint32(len(imgZ) / 2) - 1
        # Put the corresponding values in the mass, index, and handles.Noise windows
        self.start.setText(str(imgZ[self.index]))
        self.msindex.setText(str(self.index))
        self.Noise_Output_Box.setText(str(self.img_std[self.index]))
        # pull the maximum value from the image and put it in the boxes
        self.z_max, self.z_min = self.arrlims(self.intens[:, :, (round(len(imgZ) / 2 - 1))])
        self.max_int.setText(str(self.z_max))
        self.temp_max.setText(str(self.z_max))
        # set up the slider (scrollbar)
        self.zmax.setMinimum(0)
        self.zmax.setMaximum(math.ceil(self.z_max))
        self.zmax.setValue(math.ceil(self.z_max))
        self.zmax.setSingleStep(round(0.01 * self.z_max))
        # concentration map storing variable
        self.ConcMapData = [np.flip(img[:, :, round(len(imgZ) / 2) - 1], axis=0), self.x_end, self.y_end]
        # set up the initial image
        if self.con_canvas:
            self.plot_con.removeWidget(self.con_canvas)
            # self.con_cbar.remove()
        if self.con_canvas:
            self._con_ax.cla()
            # self.con_cbar.remove()
            self.con_img = self._con_ax.imshow(np.flip(img[:, :, round(len(imgZ) / 2) - 1], axis=0), cmap='jet',
                                               aspect='auto', extent=[0, self.x_end, 0, self.y_end])
            self._con_ax.set_xlabel('x, mm')
            self._con_ax.set_ylabel('y, mm')
            self.con_cbar = plt.colorbar(self.con_img)
            self._con_ax.invert_yaxis()
            self._con_ax.set_aspect('equal')
            self.con_canvas.draw()
            self.exConcflag = True
        else:
            self.con_canvas = FigureCanvas(plt.figure(tight_layout=True))
            self.con_toolbar = NavigationToolbar(self.con_canvas, self)
            self.plot_con.addWidget(self.con_toolbar)
            self.plot_con.addWidget(self.con_canvas)
            self._con_ax = self.con_canvas.figure.subplots()
            self.con_img = self._con_ax.imshow(np.flip(img[:, :, round(len(imgZ) / 2) - 1], axis=0), cmap='jet',
                                               aspect='auto', extent=[0, self.x_end, 0, self.y_end])
            self._con_ax.set_xlabel('x, mm')
            self._con_ax.set_ylabel('y, mm')
            self.con_cbar = plt.colorbar(self.con_img)
            self._con_ax.invert_yaxis()
            self._con_ax.set_aspect('equal')
            self.exConcflag = True

    # --- Executes on slider movement.
    def zmax_Callback(self):
        if isIM:
            self.temp_max.setText(str(self.zmax.sliderPosition()))
            self.scale_image()
            return 0
        self.t_max = self.zmax.sliderPosition()
        self.temp_max.setText(str(self.t_max))
        self.scale_image()

    # This function updates the image to the current index
    def scale_image(self):
        if isIM:
            x = self.chosenData

            data = []
            position = self.zmax.sliderPosition()

            for line in x:
                newLine = []
                for frame in line:
                    newFrame = 0
                    valAdded = False
                    if frame != 0:
                        for val in frame:
                            if val[1] > position:
                                valAdded = True
                                newFrame += val[1]
                        if not valAdded:
                            newFrame = 0
                    newLine.append(newFrame)
                data.append(newLine)

            numY = len(data)
            numX = len(data[0])
            xend = numX * .075
            yend = numY * .15

            if self.view:
                self.plot_con.removeWidget(self.view)
            self.view = FigureCanvas(Figure(figsize=(5, 3)))
            self.axes = self.view.figure.subplots()
            self.toolbar = NavigationToolbar(self.view, self)
            self.plot_con.addWidget(self.view)
            self.con_img = self.axes.imshow(data, cmap='jet', interpolation='gaussian',
                                            aspect=(yend / xend), extent=[0, xend, 0, yend])
            plt.colorbar(self.con_img)
            self.view.draw()
            return 0
        # plot the image
        bkg = np.zeros((len(self.y), len(self.x)))  # creates matrix for self.Background
        summarea = np.zeros((len(self.y), len(self.x)))  # creates matrix for the area sums
        self.areas = np.zeros((len(self.y), len(self.x)))
        numberpoints = 12
        for i in range(len(self.y)):
            for k in range(len(self.x)):
                sumeach = np.sum(self.intens[i, k, np.arange(self.index - numberpoints, (self.index + 70 + 1))])
                summarea[i, k] = sumeach
                sumback = np.sum(
                    self.intens[i, k, np.arange(self.index + numberpoints, (self.index + 2 * numberpoints + 1))])
                bkg[i, k] = sumback
                self.areas[i, k] = (sumeach - sumback)
        self.ConcMapData = [np.flip(self.areas, axis=0), self.x_end, self.y_end, self.z_min, self.t_max]
        if self._con_ax:
            self._con_ax.cla()
            # if self.con_cbar:
            #     print("here")
            # self.con_cbar.remove()
            if self.con_canvas:
                self.plot_con.removeWidget(self.con_canvas)
                # self.con_cbar.remove()
            clims = np.array([self.z_min, self.t_max])
            self.con_img = self._con_ax.imshow(np.flip(self.areas, axis=0), cmap='jet', aspect='auto', vmin=clims[0],
                                               vmax=clims[1], extent=[0, self.x_end, 0, self.y_end])
            self._con_ax.set_xlabel('x, mm')
            self._con_ax.set_ylabel('y, mm')
            self.con_cbar = plt.colorbar(self.con_img)
            # self.con_cbar = self.con_canvas.figure.colorbar(self.con_img)
            self._con_ax.invert_yaxis()
            self._con_ax.set_aspect('equal')
            self.con_canvas.draw()
            self.exConcflag = True

    # --- Executes on slider movement.
    def zmax_isotope_Callback(self):
        if isIM:
            self.max_iso.setText(str(self.zmax_isotope.sliderPosition()))
            self.scale_iso_image()
            return 0
        self.iso_max = self.zmax_truearr[self.zmax_isotope.sliderPosition()]
        self.iso_min = self.zmax_truearr[self.zmin_isotope.sliderPosition()]  # TODO: Is this line necessary?
        self.max_iso.setText(str(('%s' % float('%.5g' % self.iso_max))))
        self.scale_iso_image()

    # --- Executes on slider movement.
    def zmin_isotope_Callback(self):
        if isIM:
            self.min_iso.setText(str(self.zmin_isotope.sliderPosition()))
            self.scale_iso_image()
            return 0
        self.iso_max = self.zmax_truearr[self.zmax_isotope.sliderPosition()]
        self.iso_min = self.zmax_truearr[self.zmin_isotope.sliderPosition()]
        self.min_iso.setText(str(('%s' % float('%.5g' % self.iso_min))))
        self.scale_iso_image()

    # This function updates the image to the current index
    def scale_iso_image(self):
        if isIM:
            if self.chosenDataIso is None:  # This will be triggered only at the beginning
                return 0
            x = self.chosenDataIso
            data = []
            highest = self.zmax_isotope.sliderPosition()
            lowest = self.zmin_isotope.sliderPosition()

            for line in x:
                newLine = []
                for frame in line:
                    newFrame = 0
                    valAdded = False
                    if frame != 0:
                        for val in frame:
                            if highest >= val[1] >= lowest:
                                valAdded = True
                                newFrame += val[1]
                        if not valAdded:
                            newFrame = 0
                    newLine.append(newFrame)
                data.append(newLine)

            numY = len(data)
            numX = len(data[0])
            xend = numX * .075
            yend = numY * .15

            if self.viewPlusOne:
                self.plot_kin.removeWidget(self.viewPlusOne)
            if self.viewPlusTwo:
                self.plot_kin.removeWidget(self.viewPlusTwo)

            self.viewPlusOne = FigureCanvas(Figure(figsize=(5, 3)))
            self.axes = self.viewPlusOne.figure.subplots()
            self.toolbar = NavigationToolbar(self.view, self)
            self.plot_kin.addWidget(self.viewPlusOne)
            self.con_img2 = self.axes.imshow(data, cmap='inferno',
                                             aspect=(yend / xend), extent=[0, xend, 0, yend])
            plt.colorbar(self.con_img2)
            self.viewPlusOne.draw()
            return 0
        if self.massplustwo.isChecked():
            if self._kin_ax:
                self.IsotopeMapData = [np.flip(self.mplusonenormintratiowithmassplustwo, axis=0), self.x_end,
                                       self.y_end, self.iso_min, self.iso_max]
                self._kin_ax.cla()
                self.kin_cbar.remove()
                clims = np.array([self.iso_min, self.iso_max])
                self.kin_img = self._kin_ax.imshow(np.flip(self.mplusonenormintratiowithmassplustwo, axis=0),
                                                   cmap='jet', aspect='auto', vmin=clims[0], vmax=clims[1],
                                                   extent=[0, self.x_end, 0, self.y_end])
                self._kin_ax.set_xlabel('x, mm')
                self._kin_ax.set_ylabel('y, mm')
                self.kin_cbar = self.kin_canvas.figure.colorbar(self.kin_img)
                self._kin_ax.invert_yaxis()
                self._kin_ax.set_aspect('equal')
                self.kin_canvas.draw()
                self.exIsotopeflag = True
            else:
                pass
        else:
            if self._kin_ax:
                self.IsotopeMapData = [np.flip(self.mplusonenormintratio, axis=0), self.x_end, self.y_end, self.iso_min,
                                       self.iso_max]
                self._kin_ax.cla()
                self.kin_cbar.remove()
                clims = np.array([self.iso_min, self.iso_max])
                self.kin_img = self._kin_ax.imshow(np.flip(self.mplusonenormintratio, axis=0), cmap='jet',
                                                   aspect='auto', vmin=clims[0], vmax=clims[1],
                                                   extent=[0, self.x_end, 0, self.y_end])
                self._kin_ax.set_xlabel('x, mm')
                self._kin_ax.set_ylabel('y, mm')
                self.kin_cbar = self.kin_canvas.figure.colorbar(self.kin_img)
                self._kin_ax.invert_yaxis()
                self._kin_ax.set_aspect('equal')
                self.kin_canvas.draw()
                self.exIsotopeflag = True
            else:
                pass

    def export_ConcMap_Callback(self):
        # when clicking on the exportConcMap button, it will save the filename and concentration the map
        if self.exConcflag:
            self.Maps[self.exportConcMapname.text()] = self.ConcMapData
            self.refreshMaplistbox()
            self.Mapcount += 1

    def export_IsotopeMap_Callback(self):
        # when clicking on the exportConcMap button, it will save the filename and concentration the map
        if self.exIsotopeflag:
            self.Maps[self.exportIsotopeMapname.text()] = self.IsotopeMapData
            self.refreshMaplistbox()
            self.Mapcount += 1

    def Map_listbox_Callback(self, item):
        self.Map_listselect_text = item.text()

    def deleteMapbutton_Callback(self):
        if self.Mapcount == 0:
            print('There are no maps in the listbox')
        else:
            if self.Map_listselect_text:
                del self.Maps[self.Map_listselect_text]
                self.refreshMaplistbox()
                self.Mapcount -= 1
                print('Map removed')
            else:
                print('Item does not exist. Please double click on another item')

    def clearMapbutton_Callback(self):
        if self.Mapcount == 0:
            print('There are no maps in the listbox')
        else:
            self.Maps.clear()
            self.Mapcount = 0
            self.refreshMaplistbox()

    def refreshMaplistbox(self):
        if len(self.Maps) == 0:
            # set box with default text
            self.Map_listbox.clear()
            if "Map_listselect_text" in locals():
                del self.Map_listselect_text
            QListWidgetItem('Map list appears here', self.Map_listbox)
        else:
            self.Map_listbox.clear()
            if "Map_listselect_text" in locals():
                del self.Map_listselect_text
            listboxitems = list(self.Maps.keys())
            for i in range(len(listboxitems)):
                self.Map_listbox.addItem(listboxitems[i])

    def exportROI_Callback(self):
        if isIM:
            self.ROIcount = self.ROIcount + 1
            self.ROIcountbox.setText(str(self.ROIcount))
            self.ROI[self.exportROIfilename.text()] = self.binI
            ROI = self.ROI
            ROIcount = self.ROIcount
            self.refreshROIlistbox()
            return 0

        # when clicking on the exportROI button, it will increase the ROIcount by 1
        # and it will save the filename and binI to the ROI cell variable
        # then it will save them to the actual cubefilename to be retrievable next
        # time the cube is opened
        if self.rgflagROI:  # only export a ROI if one has been selected
            self.ROIcount = self.ROIcount + 1
            self.ROIcountbox.setText(str(self.ROIcount))
            self.ROI[self.exportROIfilename.text()] = self.binI
            # Would I really need to save it?
            ROI = self.ROI
            ROIcount = self.ROIcount
            self.refreshROIlistbox()

    def refreshROIlistbox(self):
        if len(self.ROI) == 0:
            # set box with default text
            self.ROI_listbox.clear()
            del self.ROI_listselect_text
            QListWidgetItem('ROI list appears here', self.ROI_listbox)
        else:
            self.ROI_listbox.clear()
            listboxitems = list(self.ROI.keys())
            for i in range(len(listboxitems)):
                self.ROI_listbox.addItem(listboxitems[i])

    # --- Executes on selection change in ROI_listbox
    def ROI_listbox_Callback(self, item):
        # this should tell me which number in list is selected
        # now I need the outline of the ROI to appear
        if isIM:
            self.ROI_listselect_text = item.text()
            self.ROI_listselect_array = self.ROI[item.text()]
            self.ROI_img_mean[item.text()] = self.img_mean
            self._con_ax = 1
            self.ROIplots[self.ROI_listselect_text] = item  # What am I doing here??
            # self.ROI[self.ROI_listselect_text] =
            return 0

        self.ROI_listselect_text = item.text()
        self.ROI_listselect_array = self.ROI[item.text()]
        self.ROI_img_mean[self.ROI_listselect_text] = self.img_mean
        structboundaries = self.boundary_tracer(
            self.ROI_listselect_array)  # collects the indices of the pixels on the boundaries and holes
        xy = structboundaries
        pixelx = self.x_end / (len(self.x) - 1)  # mm pixel dimension
        pixely = self.y_end / (len(self.y) - 1)  # mm pixel dimension
        X = (xy[:, 1] * pixelx) - pixelx  # makes image appear in right spot
        Y = (xy[:, 0] * pixely) - pixely  # makes image appear in right spot

        if self.massplusone.isChecked():
            self.includemassplustwo = False
        elif self.massplustwo.isChecked():
            self.includemassplustwo = True
        if self.has_data:
            # get average spectrum
            # this finds all the pixels (i.e. rows) that are part of the ROI
            self.f = np.argwhere(np.ravel(self.ROI_listselect_array, order='F'))[:, 0]  ### check
            # this extracts only those rows that are part of the ROI from the reshaped data, and takes their mean
            if (self.massbox.text()) != 0:
                massvalue = float(self.massbox.text())
                massindex = np.where(abs(self.z - (massvalue)) == (abs(self.z - (massvalue))).min())[0][0]
                massplusoneindex = np.where(abs(self.z - (massvalue + 1)) == (abs(self.z - (massvalue + 1))).min())[0][
                    0]
                massplustwoindex = np.where(abs(self.z - (massvalue + 2)) == (abs(self.z - (massvalue + 2))).min())[0][
                    0]
                # finds the index of the closest value to the massplusone, massplustwo and mass value
                mplusone = self.intens[:, :, massplusoneindex]  # this gets the massplusone plane
                mplustwo = self.intens[:, :, massplustwoindex]  # this gets the massplustwo plane
                m = self.intens[:, :, massindex]  # this gets mass plane
                zmax = np.max(m[:])  # finds the max overall in mass selected plane
                threshold = zmax / 7  # sets the threshold for the peaks to be real above a fifth of maximum
                good = (m > threshold).astype(
                    int)  # creates a matrix with ones and zeros for if the statement is true or not for each value

                self.mintratio = np.zeros((len(self.y), len(self.x)))  # creates matrix for integrated ratios
                self.mintratiowithmassplustwo = np.zeros(
                    (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                self.mplusoneintratio = np.zeros((len(self.y), len(self.x)))  # creates matrix for integrated ratios
                self.mplusoneintratiowithmassplustwo = np.zeros(
                    (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                self.mplustwointratiowithmassplustwo = np.zeros(
                    (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                bkg = np.zeros((len(self.y), len(self.x)))  # creates matrix for handles.Background
                summarr = np.zeros((len(self.y), len(self.x)))  # creates matrix for the area sums
                summplusonearr = np.zeros((len(self.y), len(self.x)))  # creates matrix for the sumplusone area sums
                summplustwoarr = np.zeros((len(self.y), len(self.x)))  # creates matrix for the sumplustwo area sums
                numberpoints = 12

                for i in range(len(self.y)):
                    for k in range(len(self.x)):
                        if self.good[i, k]:
                            summ = np.sum(
                                self.intens[
                                    i, k, np.arange(self.index - numberpoints, (self.index + numberpoints + 1))])
                            summarr[i, k] = summ
                            summplusone = np.sum(self.intens[i, k, np.arange(massplusoneindex - numberpoints,
                                                                             (massplusoneindex + numberpoints + 1))])
                            summplusonearr[i, k] = summplusone
                            summplustwo = np.sum(self.intens[i, k, np.arange(massplustwoindex - numberpoints,
                                                                             (massplustwoindex + numberpoints + 1))])
                            summplustwoarr[i, k] = summplustwo
                            sumback = np.sum(
                                self.intens[
                                    i, k, np.arange(self.index + numberpoints, (self.index + 2 * numberpoints + 1))])
                            bkg[i, k] = sumback
                            # getting all the integrated ratios for the main mass peak
                            self.mintratio[i, k] = (summ - sumback) / ((summ - sumback) + (summplusone - sumback))
                            self.mintratiowithmassplustwo[i, k] = (summ - sumback) / (
                                    (summ - sumback) + (summplusone - sumback) + (summplustwo - sumback))
                            # getting all the integrated ratios for the mass peak plus one
                            self.mplusoneintratio[i, k] = (summplusone - sumback) / (
                                    (summ - sumback) + (summplusone - sumback))
                            self.mplusoneintratiowithmassplustwo[i, k] = (summplusone - sumback) / (
                                    (summ - sumback) + (summplusone - sumback) + (summplustwo - sumback))
                            # getting all the integrated ratios for the mass peak plus two
                            self.mplustwointratiowithmassplustwo[i, k] = (summplustwo - sumback) / (
                                    (summ - sumback) + (summplusone - sumback) + (summplustwo - sumback))
                        else:
                            self.mintratio[i, k] = 0
                            self.mintratiowithmassplustwo[i, k] = 0
                            self.mplusoneintratio[i, k] = 0
                            self.mplusoneintratiowithmassplustwo[i, k] = 0
                            self.mplustwointratiowithmassplustwo[i, k] = 0
                # checked up to this point that everything is correct
                if self.includemassplustwo:
                    ROI_mintratiowithmassplustwo = np.zeros(
                        (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    ROI_mplusoneintratiowithmassplustwo = np.zeros(
                        (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    ROI_mplustwointratiowithmassplustwo = np.zeros(
                        (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    for i in range(len(self.y)):
                        for k in range(len(self.x)):
                            if self.ROI_listselect_array[i, k]:
                                ROI_mintratiowithmassplustwo[i, k] = self.mintratiowithmassplustwo[i, k]
                                ROI_mplusoneintratiowithmassplustwo[i, k] = self.mplusoneintratiowithmassplustwo[i, k]
                                ROI_mplustwointratiowithmassplustwo[i, k] = self.mplustwointratiowithmassplustwo[i, k]
                            else:
                                ROI_mintratiowithmassplustwo[i, k] = 0
                                ROI_mplusoneintratiowithmassplustwo[i, k] = 0
                                ROI_mplustwointratiowithmassplustwo[i, k] = 0
                    ROImnonzeroswithmassplustwo = ROI_mintratiowithmassplustwo[np.nonzero(ROI_mintratiowithmassplustwo)]
                    ROImplusonenonzeroswithmassplustwo = ROI_mplusoneintratiowithmassplustwo[
                        np.nonzero(ROI_mplusoneintratiowithmassplustwo)]
                    ROImplustwononzeroswithmassplustwo = ROI_mplustwointratiowithmassplustwo[
                        np.nonzero(ROI_mplustwointratiowithmassplustwo)]
                    ROImaveragewithmassplustwo = np.mean(ROImnonzeroswithmassplustwo)
                    ROImplusoneaveragewithmassplustwo = np.mean(ROImplusonenonzeroswithmassplustwo)
                    ROImplustwoaveragewithmassplustwo = np.mean(ROImplustwononzeroswithmassplustwo)
                    ROImstdevwithmassplustwo = np.std(ROImnonzeroswithmassplustwo, ddof=1)
                    ROImplusonestdevwithmassplustwo = np.std(ROImplusonenonzeroswithmassplustwo, ddof=1)
                    ROImplustwostdevwithmassplustwo = np.std(ROImplustwononzeroswithmassplustwo, ddof=1)
                    NumberinROI = len(ROImnonzeroswithmassplustwo)
                    self.Msumratio.setText(str('%s' % float('%.5g' % ROImaveragewithmassplustwo)))  # set isotope ratio
                    self.Msumstandard_error.setText(
                        str('%s' % float('%.5g' % ROImstdevwithmassplustwo)))  # set error box
                    self.numberpoints.setText(str(NumberinROI))  # set number of points
                    self.Mplusonesumratio.setText(
                        str('%s' % float('%.5g' % ROImplusoneaveragewithmassplustwo)))  # set isotope ratio
                    self.Mplusonesumstandard_error.setText(
                        str('%s' % float('%.5g' % ROImplusonestdevwithmassplustwo)))  # set error box
                    self.Mplustwosumratio.setText(
                        str('%s' % float('%.5g' % ROImplustwoaveragewithmassplustwo)))  # set isotope ratio
                    self.Mplustwosumstandard_error.setText(
                        str('%s' % float('%.5g' % ROImplustwostdevwithmassplustwo)))  # set error box

                else:
                    ROI_mintratio = np.zeros((len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    ROI_mplusoneintratio = np.zeros((len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    for i in range(len(self.y)):
                        for k in range(len(self.x)):
                            if self.ROI_listselect_array[i, k]:
                                ROI_mintratio[i, k] = self.mintratio[i, k]
                                ROI_mplusoneintratio[i, k] = self.mplusoneintratio[i, k]
                            else:
                                ROI_mintratio[i, k] = 0
                                ROI_mplusoneintratio[i, k] = 0
                    ROImnonzeros = ROI_mintratio[np.nonzero(ROI_mintratio)]
                    ROImplusonenonzeros = ROI_mplusoneintratio[np.nonzero(ROI_mplusoneintratio)]
                    ROImaverage = np.mean(ROImnonzeros)
                    ROImplusoneaverage = np.mean(ROImplusonenonzeros)
                    ROImstdev = np.std(ROImnonzeros, ddof=1)
                    ROImplusonestdev = np.std(ROImplusonenonzeros, ddof=1)
                    NumberinROI = len(ROImnonzeros)
                    self.Msumratio.setText(str('%s' % float('%.5g' % ROImaverage)))  # set the isotoperatio box
                    self.Msumstandard_error.setText(str('%s' % float('%.5g' % ROImstdev)))  # set error box
                    self.numberpoints.setText(str(NumberinROI))  # set number of points
                    self.Mplusonesumratio.setText(str('%s' % float('%.5g' % ROImplusoneaverage)))  # set isotope ratio
                    self.Mplusonesumstandard_error.setText(
                        str('%s' % float('%.5g' % ROImplusonestdev)))  # set error box
                    self.Mplustwosumratio.setText("N/A")  # set isotope ratio
                    self.Mplustwosumstandard_error.setText("N/A")  # set error box
            else:
                print("No data loaded")

        if self._con_ax:
            if self.ROI_listselect_text in self.ROIplots:
                pass
            else:
                self.ROIplots[self.ROI_listselect_text] = self._con_ax.plot(X, Y, 'wo', ms=2)

    # loads the spectra from the selected ROI in the ROI_listbox
    # --- Executes on button press in 'Load selected ROI spectra' button
    def importROI_Callback(self):
        if (self.ROI_listselect_text == ""):
            print("No item selected")
        else:
            self.rgflag = False
            self.ptflag = False
            self.rgflagROI = False
            self.rgflagROIexportMS = True

            structboundaries = self.boundary_tracer(
                self.ROI_listselect_array)  # collects the indices of the pixels on the boundaries and holes
            xy = structboundaries
            pixelx = self.x_end / (len(self.x) - 1)  # mm pixel dimension
            pixely = self.y_end / (len(self.y) - 1)  # mm pixel dimension
            X = (xy[:, 1] * pixelx) - pixelx  # makes image appear in right spot
            Y = (xy[:, 0] * pixely) - pixely  # makes image appear in right spot

            if self._con_ax:
                if self.ROI_listselect_text in self.ROIplots:
                    print('')
                else:
                    self.ROIplots[self.ROI_listselect_text] = self._con_ax.plot(X, Y, 'wo', ms=2)
                    self.ROIplots_plotselect = self.ROIplots[self.ROI_listselect_text]

                # now I need the mass spectrum average and the ratios to be given over to the right
            if self.massplusone.isChecked():
                self.includemassplustwo = False
            elif self.massplustwo.isChecked():
                self.includemassplustwo = True
            if self.has_data:
                # get average spectrum
                # this finds all the pixels (i.e. rows) that are part of the ROI
                self.f = np.argwhere(np.ravel(self.ROI_listselect_array, order='F'))[:, 0]  ### check
                # this extracts only those rows that are part of the ROI from the reshaped data, and takes their mean
                if (self.massbox.text()) != 0:
                    massvalue = float(self.massbox.text())
                    massindex = np.where(abs(self.z - (massvalue)) == (abs(self.z - (massvalue))).min())[0][0]
                    massplusoneindex = \
                        np.where(abs(self.z - (massvalue + 1)) == (abs(self.z - (massvalue + 1))).min())[0][0]
                    massplustwoindex = \
                        np.where(abs(self.z - (massvalue + 2)) == (abs(self.z - (massvalue + 2))).min())[0][0]
                    # finds the index of the closest value to the massplusone, massplustwo and mass value
                    mplusone = self.intens[:, :, massplusoneindex]  # this gets the massplusone plane
                    mplustwo = self.intens[:, :, massplustwoindex]  # this gets the massplustwo plane
                    m = self.intens[:, :, massindex]  # this gets mass plane
                    zmax = np.max(m[:])  # finds the max overall in mass selected plane
                    threshold = zmax / 7  # sets the threshold for the peaks to be real above a fifth of maximum
                    good = (m > threshold).astype(
                        int)  # creates a matrix with ones and zeros for if the statement is true or not for each value

                    self.mintratio = np.zeros((len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    self.mintratiowithmassplustwo = np.zeros(
                        (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    self.mplusoneintratio = np.zeros((len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    self.mplusoneintratiowithmassplustwo = np.zeros(
                        (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    self.mplustwointratiowithmassplustwo = np.zeros(
                        (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                    bkg = np.zeros((len(self.y), len(self.x)))  # creates matrix for handles.Background
                    summarr = np.zeros((len(self.y), len(self.x)))  # creates matrix for the area sums
                    summplusonearr = np.zeros((len(self.y), len(self.x)))  # creates matrix for the sumplusone area sums
                    summplustwoarr = np.zeros((len(self.y), len(self.x)))  # creates matrix for the sumplustwo area sums
                    numberpoints = 12

                    for i in range(len(self.y)):
                        for k in range(len(self.x)):
                            if self.good[i, k]:
                                summ = np.sum(
                                    self.intens[
                                        i, k, np.arange(self.index - numberpoints, (self.index + numberpoints + 1))])
                                summarr[i, k] = summ
                                summplusone = np.sum(self.intens[i, k, np.arange(massplusoneindex - numberpoints,
                                                                                 (
                                                                                         massplusoneindex + numberpoints + 1))])
                                summplusonearr[i, k] = summplusone
                                summplustwo = np.sum(self.intens[i, k, np.arange(massplustwoindex - numberpoints,
                                                                                 (
                                                                                         massplustwoindex + numberpoints + 1))])
                                summplustwoarr[i, k] = summplustwo
                                sumback = np.sum(
                                    self.intens[i, k, np.arange(self.index + numberpoints,
                                                                (self.index + 2 * numberpoints + 1))])
                                bkg[i, k] = sumback
                                # getting all the integrated ratios for the main mass peak
                                self.mintratio[i, k] = (summ - sumback) / ((summ - sumback) + (summplusone - sumback))
                                self.mintratiowithmassplustwo[i, k] = (summ - sumback) / (
                                        (summ - sumback) + (summplusone - sumback) + (summplustwo - sumback))
                                # getting all the integrated ratios for the mass peak plus one
                                self.mplusoneintratio[i, k] = (summplusone - sumback) / (
                                        (summ - sumback) + (summplusone - sumback))
                                self.mplusoneintratiowithmassplustwo[i, k] = (summplusone - sumback) / (
                                        (summ - sumback) + (summplusone - sumback) + (summplustwo - sumback))
                                # getting all the integrated ratios for the mass peak plus two
                                self.mplustwointratiowithmassplustwo[i, k] = (summplustwo - sumback) / (
                                        (summ - sumback) + (summplusone - sumback) + (summplustwo - sumback))
                            else:
                                self.mintratio[i, k] = 0
                                self.mintratiowithmassplustwo[i, k] = 0
                                self.mplusoneintratio[i, k] = 0
                                self.mplusoneintratiowithmassplustwo[i, k] = 0
                                self.mplustwointratiowithmassplustwo[i, k] = 0
                    # checked up to this point that everything is correct
                    if self.includemassplustwo:
                        ROI_mintratiowithmassplustwo = np.zeros(
                            (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                        ROI_mplusoneintratiowithmassplustwo = np.zeros(
                            (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                        ROI_mplustwointratiowithmassplustwo = np.zeros(
                            (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                        for i in range(len(self.y)):
                            for k in range(len(self.x)):
                                if self.ROI_listselect_array[i, k]:
                                    ROI_mintratiowithmassplustwo[i, k] = self.mintratiowithmassplustwo[i, k]
                                    ROI_mplusoneintratiowithmassplustwo[i, k] = self.mplusoneintratiowithmassplustwo[
                                        i, k]
                                    ROI_mplustwointratiowithmassplustwo[i, k] = self.mplustwointratiowithmassplustwo[
                                        i, k]
                                else:
                                    ROI_mintratiowithmassplustwo[i, k] = 0
                                    ROI_mplusoneintratiowithmassplustwo[i, k] = 0
                                    ROI_mplustwointratiowithmassplustwo[i, k] = 0
                        ROImnonzeroswithmassplustwo = ROI_mintratiowithmassplustwo[
                            np.nonzero(ROI_mintratiowithmassplustwo)]
                        ROImplusonenonzeroswithmassplustwo = ROI_mplusoneintratiowithmassplustwo[
                            np.nonzero(ROI_mplusoneintratiowithmassplustwo)]
                        ROImplustwononzeroswithmassplustwo = ROI_mplustwointratiowithmassplustwo[
                            np.nonzero(ROI_mplustwointratiowithmassplustwo)]
                        ROImaveragewithmassplustwo = np.mean(ROImnonzeroswithmassplustwo)
                        ROImplusoneaveragewithmassplustwo = np.mean(ROImplusonenonzeroswithmassplustwo)
                        ROImplustwoaveragewithmassplustwo = np.mean(ROImplustwononzeroswithmassplustwo)
                        ROImstdevwithmassplustwo = np.std(ROImnonzeroswithmassplustwo, ddof=1)
                        ROImplusonestdevwithmassplustwo = np.std(ROImplusonenonzeroswithmassplustwo, ddof=1)
                        ROImplustwostdevwithmassplustwo = np.std(ROImplustwononzeroswithmassplustwo, ddof=1)
                        NumberinROI = len(ROImnonzeroswithmassplustwo)
                        self.Msumratio.setText(
                            str('%s' % float('%.5g' % ROImaveragewithmassplustwo)))  # set isotope ratio
                        self.Msumstandard_error.setText(
                            str('%s' % float('%.5g' % ROImstdevwithmassplustwo)))  # set error box
                        self.numberpoints.setText(str(NumberinROI))  # set number of points
                        self.Mplusonesumratio.setText(
                            str('%s' % float('%.5g' % ROImplusoneaveragewithmassplustwo)))  # set isotope ratio
                        self.Mplusonesumstandard_error.setText(
                            str('%s' % float('%.5g' % ROImplusonestdevwithmassplustwo)))  # set error box
                        self.Mplustwosumratio.setText(
                            str('%s' % float('%.5g' % ROImplustwoaveragewithmassplustwo)))  # set isotope ratio
                        self.Mplustwosumstandard_error.setText(
                            str('%s' % float('%.5g' % ROImplustwostdevwithmassplustwo)))  # set error box

                    else:
                        ROI_mintratio = np.zeros((len(self.y), len(self.x)))  # creates matrix for integrated ratios
                        ROI_mplusoneintratio = np.zeros(
                            (len(self.y), len(self.x)))  # creates matrix for integrated ratios
                        for i in range(len(self.y)):
                            for k in range(len(self.x)):
                                if self.ROI_listselect_array[i, k]:
                                    ROI_mintratio[i, k] = self.mintratio[i, k]
                                    ROI_mplusoneintratio[i, k] = self.mplusoneintratio[i, k]
                                else:
                                    ROI_mintratio[i, k] = 0
                                    ROI_mplusoneintratio[i, k] = 0
                        ROImnonzeros = ROI_mintratio[np.nonzero(ROI_mintratio)]
                        ROImplusonenonzeros = ROI_mplusoneintratio[np.nonzero(ROI_mplusoneintratio)]
                        ROImaverage = np.mean(ROImnonzeros)
                        ROImplusoneaverage = np.mean(ROImplusonenonzeros)
                        ROImstdev = np.std(ROImnonzeros, ddof=1)
                        ROImplusonestdev = np.std(ROImplusonenonzeros, ddof=1)
                        NumberinROI = len(ROImnonzeros)
                        self.Msumratio.setText(str('%s' % float('%.5g' % ROImaverage)))  # set the isotoperatio box
                        self.Msumstandard_error.setText(str('%s' % float('%.5g' % ROImstdev)))  # set error box
                        self.numberpoints.setText(str(NumberinROI))  # set number of points
                        self.Mplusonesumratio.setText(
                            str('%s' % float('%.5g' % ROImplusoneaverage)))  # set isotope ratio
                        self.Mplusonesumstandard_error.setText(
                            str('%s' % float('%.5g' % ROImplusonestdev)))  # set error box
                        self.Mplustwosumratio.setText("N/A")  # set isotope ratio
                        self.Mplustwosumstandard_error.setText("N/A")  # set error box
                else:
                    print("No data loaded")
                self.img_mean = np.mean(self.reshaped[self.f, :], axis=0)
                self.img_std = np.std(self.reshaped[self.f, :], ddof=1)
                self.ROI_img_mean[self.ROI_listselect_text] = self.img_mean
                # Plot the spectrum averaged over the ROI
                self.spectra_df = pd.DataFrame({'m/z': self.z, 'intensity': self.img_mean})
                self.spectra_df['max'] = self.spectra_df.iloc[
                    signal.argrelextrema(self.spectra_df.intensity.values, np.greater_equal,
                                         order=self.checkpoint_maxima)[0]]['intensity']
                if self.spectra_canvas:
                    self._spectra_ax.cla()
                    self._spectra_ax.plot(self.z, self.img_mean, 'k', linewidth=0.3, picker=True)
                    self._spectra_ax.set_title(
                        'Average Spectrum Across Selected Region from ROI listbox: ' + self.ROI_listselect_text)
                    self._spectra_ax.set_xlabel('m/z')
                    self._spectra_ax.set_ylabel('intensity')
                else:
                    print("No data loaded")
            else:
                print("No data loaded")

    def exportROI_spectra_val_Callback(self):
        if (self.ROI_listselect_text == ""):
            print("No item selected. Double click an item in the ROI listbox")
        else:
            img_mean = self.ROI_img_mean[self.ROI_listselect_text]
            path = QFileDialog.getSaveFileName(self, 'Save CSV', os.getenv('HOME'), 'CSV(*.csv)')
            if path[0] != '':
                if isIM:
                    with open(path[0], 'w', newline='') as csv_file:
                        writer = csv.writer(csv_file, dialect='excel', delimiter=',', lineterminator='\n')
                        writer.writerow(["M/Z Value", "Intensity", "Drift Time", "Line", "Frame Num"])
                        for line in self.ROIData:
                            writer.writerow(line)
                else:
                    with open(path[0], 'w', newline='') as csv_file:
                        writer = csv.writer(csv_file, dialect='excel', delimiter=',', lineterminator='\n')
                        writer.writerow(('m/z', 'intensity'))
                        for row in range(len(self.z)):
                            writer.writerow((self.z[row], img_mean[row]))
            else:
                print("No location selected. Please select a location")

    def deleteROIbutton_Callback(self):
        if self.ROIcount == 0:
            print('There are no ROIs in the listbox')
        else:
            if self._con_ax:
                if self.ROI_listselect_text in self.ROIplots:
                    # How to get the plot removed from the graph?
                    if not isIM:
                        self.ROIplots[self.ROI_listselect_text].pop(0).remove()
                    del self.ROIplots[self.ROI_listselect_text]
                    del self.ROI[self.ROI_listselect_text]
                    del self.ROI_img_mean[self.ROI_listselect_text]
                    self.ROIcount = self.ROIcount - 1
                    self.ROIcountbox.setText(str(self.ROIcount))
                    self.refreshROIlistbox()
                    print('plot removed')
                    if isIM:
                        self.im_point()
                        # This is in order to remove the blue lines from the figure.
                        # If the blue lines aren't a problem, this can be taken out.
                        # for i in range(2):
                        #     self.h.lines.pop(-1).remove()
                        #     self.h.xcoords.pop()
                        #     self.h.ycoords.pop()
                        #     self.h.previous_point = self.h.xcoords[-1], self.h.ycoords[-1]
                        #     self.h.line = self.h.lines[-1]
                        #     if len(self.h.lines) == 0:
                        #         self.h.line = None
                        #         return
                        # self.h.lines.pop(-1).remove()
                        # self.h.xcoords.pop()
                        # self.h.ycoords.pop()
                        # self.h.disconnect()
                else:
                    print('Item does not exist. Please double click on another item')
            # if isIM:
            #     if self.ROI_listselect_text in self.ROIplots:
            #         self.ROIplots[self.ROI_listselect_text].pop(0).remove()
            #         del self.ROIplots[self.ROI_listselect_text]
            #         del self.ROI[self.ROI_listselect_text]
            #         del self.ROI_img_mean[self.ROI_listselect_text]
            #         self.ROIcount = self.ROIcount - 1
            #         self.ROIcountbox.setText(str(self.ROIcount))
            #         self.refreshROIlistbox()

    # --- Executes on button press in find_file_mzOI.
    def find_file_mzOI_Callback(self):
        self.fName_mzOI = QFileDialog.getOpenFileName(self, 'Pick list: m/z of interest', filter='*.csv')
        self.mzOI_listname.setText(self.fName_mzOI[0])
        self.fName_mzOI_flag = True

    def mzOI_extractMap_Callback(self):
        if self.exSpecflag:
            if self.fName_mzOI_flag:
                filename_mzOI = self.fName_mzOI[0]
                print("Working to read list: m/z of interest")
                fileID_mzOI = open(filename_mzOI)
                data = pd.read_csv(fileID_mzOI).to_numpy()[:, 0]
                fileID_mzOI.close()
                err = 0
                self.mzOI_index = np.zeros(len(data))
                temp_mzOI_index = []
                threshold = float(self.pick_mzthreshold.text())
                for i in range(len(data)):
                    for j in range(len(self.spectra_df['max'])):
                        if (pd.notnull(self.spectra_df['max'][j])) and (self.spectra_df['max'][j] != 0):
                            err = abs(
                                (data[i] - self.spectra_df['m/z'][j]) / self.spectra_df['m/z'][j]) * self.err_multp
                            if (err < threshold):
                                temp_mzOI_index.append(j)
                    if (len(temp_mzOI_index) != 0):
                        self.mzOI_index[i] = np.min(temp_mzOI_index)
                    else:
                        self.mzOI_index[i] = np.nan
                    temp_mzOI_index.clear()
                print("index of mzOI is prepared")
                mzOI_index_display = pd.DataFrame({'m/z of interest': data, 'index': self.mzOI_index})
                print(mzOI_index_display)

                if (len(self.mzOI_index) != 0):
                    for i in range(len(self.mzOI_index)):
                        if (self.mzOI_index[i] != np.nan):
                            self.index = int(self.mzOI_index[i])
                            bkg = np.zeros((len(self.y), len(self.x)))
                            summarea = np.zeros((len(self.y), len(self.x)))
                            self.areas = np.zeros((len(self.y), len(self.x)))
                            numberpoints = 12
                            for j in range(len(self.y)):
                                for k in range(len(self.x)):
                                    sumeach = np.sum(
                                        self.intens[j, k, np.arange(self.index - numberpoints, (self.index + 70 + 1))])
                                    summarea[j, k] = sumeach
                                    sumback = np.sum(
                                        self.intens[
                                            j, k, np.arange(self.index + numberpoints,
                                                            (self.index + 2 * numberpoints + 1))])
                                    bkg[j, k] = sumback
                                    self.areas[j, k] = (sumeach - sumback)
                            # pull the maximum value from the image and put it in the boxes
                            self.z_max, self.z_min = self.arrlims(self.areas)
                            self.ConcMapData = [np.flip(self.areas, axis=0), self.x_end, self.y_end]
                            self.Maps[str(data[i]) + ' Conc'] = self.ConcMapData
                            self.Mapcount += 1
                    self.refreshMaplistbox()
                else:
                    print("There is no matching m/z")
            else:
                print("Load a .csv file of m/z of interest")
        else:
            print("spectra doesn't exist")

    # not in original MATLAB code
    def data_cursor_click(self, event):
        if isinstance(event.artist, Line2D):
            thisline = event.artist
            self.x_picked = thisline.get_xdata()
            self.y_picked = thisline.get_ydata()
            ind = event.ind
            self.index = np.where(self.z == self.x_picked[ind][0])[0][0]
            self.start.setText(str(self.z[self.index]))
            self.msindex.setText(str(self.index))
            self.Noise_Output_Box.setText(str(self.img_std[self.index]))

            if (pd.notnull(self.spectra_df['max'][self.index])) and (self.spectra_df['max'][self.index] != 0):
                threshold = float(self.pick_IDthreshold.text())
                for i in range(len(self.ids_pd['m/z'])):
                    err = abs((self.x_picked[ind][0] - self.ids_pd['m/z'][i]) / self.ids_pd['m/z'][i]) * self.err_multp
                    # print(err)
                    if (err < threshold):
                        ID_in = self.ids_pd['Lipid ID'][i]
                        break
                    else:
                        ID_in = 'not defined'
                self.ID_Output_Box.setText(str(ID_in))
                # Annotate
                x_sig = '%s' % float('%.6g' % self.x_picked[ind][0])
                y_sig = '%s' % float('%.6g' % self.y_picked[ind][0])
                self.annotate_spectra_ID(x_sig, y_sig, ID_in)
            else:
                x_sig = '%s' % float('%.6g' % self.x_picked[ind][0])
                y_sig = '%s' % float('%.6g' % self.y_picked[ind][0])
                self.annotate_spectra(x_sig, y_sig)
        else:
            indexes = event.ind
            i = indexes[0]
            theXmOverZ = event.mouseevent.lastevent.xdata
            theY = event.mouseevent.lastevent.ydata

            self.start.setText("%.5f" % theXmOverZ)
            if isIM:
                self.IM_spectra_annotation(theXmOverZ, theY)
                return 0

    # This is a function that calculates how far away a value must be to be defined as a different lipid.
    def ppm_calc(self, mzVal):
        if self.ppm_min.value() == 0:
            return 0  # is this correct??
        frac = (mzVal * self.ppm_min.value()) / 1000000  # This is the ppm function engineered to find a value.

        # Before I had this, but it wasn't what I needed
        # y = mzVal - frac
        # return y
        return frac  # Is this right??

    def IM_spectra_annotation(self, mz, intensity):
        lipid_map = {}
        lipid_id = "not defined"
        if self.ids_pd is not None:
            mz_vals = self.ids_pd['m/z']
            lipid_ids = self.ids_pd['Lipid ID']
            for i in range(mz_vals.size):
                x = mz_vals[i]
                if (mz + .5) > x > (mz - .5):
                    diff = abs(mz - x)
                    lipid_map[diff] = lipid_ids[i]
        if len(lipid_map.keys()) > 0:
            key = min(lipid_map.keys())
            lipid_id = lipid_map[key]

        mz_vals = np.asarray(self.mzVals)
        drifts = np.asarray(self.drifts)

        diff = self.ppm_calc(mz)
        vals = np.where(mz + diff > mz_vals, mz_vals, 0)
        vals = np.where(mz - diff < vals, drifts, 0)
        # TODO: These +- need to be the line between m/z.
        #  How far apart should they be to be considered different lipids?

        abc = vals.nonzero()

        finals = drifts[abc]
        if len(finals) != 0:
            range_low = min(finals)
            range_high = max(finals)
        else:
            range_low = 0
            range_high = "Error. Try a lower minimum ppm."

        self.ID_Output_Box.setText(lipid_id)

        if self.annotation is not None:
            self.annotation.remove()

        self.annotation = self._spectra_ax.annotate(
            "X = {0:.4f}\nY = {1:.4f}\nID = {2}\nDrift Range = {3}-{4}".format
            (mz, intensity, lipid_id, range_low, range_high), xy=(mz, intensity), xycoords='data',
            va='bottom', ha='left',
            bbox=dict(boxstyle='square, pad=0.3', facecolor='white'))
        self.spectra_canvas.draw()

    def data_cursor_key(self, event):
        if (event.key == 'left'):
            self.index = self.index - 1
            self.start.setText(str(self.z[self.index]))
            self.msindex.setText(str(self.index))
            self.Noise_Output_Box.setText(str(self.img_std[self.index]))
            x_val = self.z[self.index]
            y_val = self.img_mean[self.index]
            x_key = '%s' % float('%.6g' % x_val)
            y_key = '%s' % float('%.6g' % y_val)

            if (pd.notnull(self.spectra_df['max'][self.index])) and (self.spectra_df['max'][self.index] != 0):
                threshold = float(self.pick_IDthreshold.text())
                for i in range(len(self.ids_pd['m/z'])):
                    err = abs((x_val - self.ids_pd['m/z'][i]) / self.ids_pd['m/z'][i]) * self.err_multp
                    # print(err)
                    if (err < threshold):
                        ID_in = self.ids_pd['Lipid ID'][i]
                        break
                    else:
                        ID_in = 'not defined'
                self.ID_Output_Box.setText(str(ID_in))
                # Annotate
                self.annotate_spectra_ID(x_key, y_key, ID_in)
            else:
                self.annotate_spectra(x_key, y_key)

        elif (event.key == 'right'):
            self.index = self.index + 1
            self.start.setText(str(self.z[self.index]))
            self.msindex.setText(str(self.index))
            self.Noise_Output_Box.setText(str(self.img_std[self.index]))
            x_val = self.z[self.index]
            y_val = self.img_mean[self.index]
            x_key = '%s' % float('%.6g' % x_val)
            y_key = '%s' % float('%.6g' % y_val)

            if (pd.notnull(self.spectra_df['max'][self.index])) and (self.spectra_df['max'][self.index] != 0):
                threshold = float(self.pick_IDthreshold.text())
                for i in range(len(self.ids_pd['m/z'])):
                    err = abs((x_val - self.ids_pd['m/z'][i]) / self.ids_pd['m/z'][i]) * self.err_multp
                    # print(err)
                    if (err < threshold):
                        ID_in = self.ids_pd['Lipid ID'][i]
                        break
                    else:
                        ID_in = 'not defined'
                self.ID_Output_Box.setText(str(ID_in))
                # Annotate
                self.annotate_spectra_ID(x_key, y_key, ID_in)
            else:
                self.annotate_spectra(x_key, y_key)

    def annotate_spectra(self, x_in, y_in):
        if self._spectra_ax:
            self._spectra_ax.cla()
            self._spectra_ax.plot(self.z, self.img_mean, 'k', linewidth=0.3, picker=True)
            self._spectra_ax.set_title('Average Spectrum Across Selected Region')
            self._spectra_ax.set_xlabel('m/z')
            self._spectra_ax.set_ylabel('intensity')
            self.spectra_ann = self._spectra_ax.annotate("X = {0}\nY = {1}".format(x_in, y_in),
                                                         xy=(self.z[self.index], self.img_mean[self.index]),
                                                         xycoords='data',
                                                         va="bottom", ha="left",
                                                         bbox=dict(boxstyle="square, pad=0.3", facecolor="white"),
                                                         )
            self.spectra_canvas.draw()

    def annotate_spectra_ID(self, x_in, y_in, ID_in):
        if self._spectra_ax:
            self._spectra_ax.cla()
            self._spectra_ax.plot(self.z, self.img_mean, 'k', linewidth=0.3, picker=True)
            self._spectra_ax.set_title('Average Spectrum Across Selected Region')
            self._spectra_ax.set_xlabel('m/z')
            self._spectra_ax.set_ylabel('intensity')
            self.spectra_ann = self._spectra_ax.annotate("X = {0}\nY = {1}\nID: {2}".format(x_in, y_in, ID_in),
                                                         xy=(self.z[self.index], self.img_mean[self.index]),
                                                         xycoords='data',
                                                         va="bottom", ha="left",
                                                         bbox=dict(boxstyle="square, pad=0.3", facecolor="white"),
                                                         )
            self.spectra_canvas.draw()

    def boundary_tracer(self, arr):
        indices_list = []
        for i in range(np.shape(arr)[0]):
            for j in range(np.shape(arr)[1]):
                if arr[i, j] == 1:
                    if i == 0 and j == 0:
                        if (arr[i + 1, j] == 0) or (arr[i, j + 1] == 0):
                            indices_list.append([i, j])
                    elif (i == 0) and (j == np.shape(arr)[1] - 1):
                        if (arr[i, j - 1] == 0) or (arr[i + 1, j] == 0):
                            indices_list.append([i, j])
                    elif (i == (np.shape(arr)[0] - 1)) and j == 0:
                        if (arr[i - 1, j] == 0) or (arr[i, j + 1] == 0):
                            indices_list.append([i, j])
                    elif (i == np.shape(arr)[0] - 1) and (j == np.shape(arr)[1] - 1):
                        if (arr[i - 1, j] == 0) or (arr[i, j - 1] == 0):
                            indices_list.append([i, j])
                    elif (i in range(1, np.shape(arr)[0] - 1)) and (j == 0):
                        if (arr[i - 1, j] == 0) or (arr[i, j + 1] == 0) or (arr[i + 1, j] == 0):
                            indices_list.append([i, j])
                    elif (i in range(1, np.shape(arr)[0] - 1)) and (j == np.shape(arr)[1] - 1):
                        if (arr[i - 1, j] == 0) or (arr[i, j - 1] == 0) or (arr[i + 1, j] == 0):
                            indices_list.append([i, j])
                    elif (i == 0) and (j in range(1, np.shape(arr)[1] - 1)):
                        if (arr[i, j - 1] == 0) or (arr[i, j + 1] == 0) or (arr[i + 1, j] == 0):
                            indices_list.append([i, j])
                    elif (i == np.shape(arr)[0] - 1) and (j in range(1, np.shape(arr)[1] - 1)):
                        if (arr[i - 1, j] == 0) or (arr[i, j - 1] == 0) or (arr[i, j + 1] == 0):
                            indices_list.append([i, j])
                    else:
                        if (arr[i - 1, j] == 0) or (arr[i + 1, j] == 0) or (arr[i, j - 1] == 0) or (arr[i, j + 1] == 0):
                            indices_list.append([i, j])
        indicies_array = np.array(indices_list)
        x_all = indicies_array[:, 1]
        x_init_bw = np.min(np.where(x_all == np.min(x_all)))
        init_coor = np.reshape(indicies_array[x_init_bw, :], (1, 2))
        indicies_array = np.vstack((indicies_array, init_coor))

        return indicies_array

    ##########################################
    # Multi Map Compare functions (things related to the 2nd window)
    def MultiMapCompare_Display_Callback(self):
        self.mmcWindow.Map_listbox_mmcWindow.clear()
        listboxitems = list(self.Maps.keys())
        for i in range(len(listboxitems)):
            self.mmcWindow.Map_listbox_mmcWindow.addItem(listboxitems[i])
        self.mmcWindow.display_mmcWindow()

    def MultiMapCompare_exportMapData_Callback(self):
        pickeditem = self.mmcWindow.pickitem
        if pickeditem:
            df = pd.DataFrame(self.Maps[pickeditem.text()][0])
            name = self.mmcWindow.exportMapData_filename.text()
            dirpath = QFileDialog.getExistingDirectory(self, 'Select a directory to export')
            if dirpath != '':
                pathfile = os.path.join(dirpath, name + '.csv')
                df.to_csv(pathfile, index=False)
            else:
                print('Please a directory')

    # executes when a load button is clicked
    def MultiMapCompare_LoadMap_Callback(self, button):
        if self.mmcWindow.pickitem:
            if button is self.mmcWindow.slot1_load:
                self.MultiMapCompare_LoadMap_func(self.mmcWindow.pickitem, 1)
            elif button is self.mmcWindow.slot2_load:
                self.MultiMapCompare_LoadMap_func(self.mmcWindow.pickitem, 2)
            elif button is self.mmcWindow.slot3_load:
                self.MultiMapCompare_LoadMap_func(self.mmcWindow.pickitem, 3)
            elif button is self.mmcWindow.slot4_load:
                self.MultiMapCompare_LoadMap_func(self.mmcWindow.pickitem, 4)
            elif button is self.mmcWindow.slot5_load:
                self.MultiMapCompare_LoadMap_func(self.mmcWindow.pickitem, 5)
            elif button is self.mmcWindow.slot6_load:
                self.MultiMapCompare_LoadMap_func(self.mmcWindow.pickitem, 6)
        else:
            print("Choose an item from the listbox")

    def MultiMapCompare_LoadMap_func(self, pickeditem, num):
        self.mmcWindow.map_packet[num][6].setText(pickeditem.text())
        x_end = self.Maps[pickeditem.text()][1]
        y_end = self.Maps[pickeditem.text()][2]
        if self.mmcWindow.map_packet[num][0]:
            self.mmcWindow.map_packet[num][1].cla()
            self.mmcWindow.map_packet[num][4].remove()
            if (len(self.Maps[pickeditem.text()])) == 3:
                self.mmcWindow.map_packet[num][2] = self.mmcWindow.map_packet[num][1].imshow(
                    self.Maps[pickeditem.text()][0], cmap='jet',
                    aspect='auto', extent=[0, x_end, 0, y_end])
            elif (len(self.Maps[pickeditem.text()])) == 5:
                zmin = self.Maps[pickeditem.text()][3]
                tmax = self.Maps[pickeditem.text()][4]
                self.mmcWindow.map_packet[num][2] = self.mmcWindow.map_packet[num][1].imshow(
                    self.Maps[pickeditem.text()][0], cmap='jet',
                    aspect='auto', vmin=zmin, vmax=tmax, extent=[0, x_end, 0, y_end])
            self.mmcWindow.map_packet[num][1].set_xlabel('x, mm')
            self.mmcWindow.map_packet[num][1].set_ylabel('y, mm')
            self.mmcWindow.map_packet[num][4] = self.mmcWindow.map_packet[num][0].figure.colorbar(
                self.mmcWindow.map_packet[num][2])
            self.mmcWindow.map_packet[num][1].invert_yaxis()
            self.mmcWindow.map_packet[num][1].set_aspect('equal')
            self.mmcWindow.map_packet[num][0].draw()
        else:
            self.mmcWindow.map_packet[num][0] = FigureCanvas(plt.figure(tight_layout=True))
            self.mmcWindow.map_packet[num][3] = NavigationToolbar(self.mmcWindow.map_packet[num][0], self)
            self.mmcWindow.map_packet[num][5].addWidget(self.mmcWindow.map_packet[num][3])
            self.mmcWindow.map_packet[num][5].addWidget(self.mmcWindow.map_packet[num][0])
            self.mmcWindow.map_packet[num][1] = self.mmcWindow.map_packet[num][0].figure.subplots()
            if (len(self.Maps[pickeditem.text()])) == 3:
                self.mmcWindow.map_packet[num][2] = self.mmcWindow.map_packet[num][1].imshow(
                    self.Maps[pickeditem.text()][0], cmap='jet',
                    aspect='auto', extent=[0, x_end, 0, y_end])
            elif (len(self.Maps[pickeditem.text()])) == 5:
                zmin = self.Maps[pickeditem.text()][3]
                tmax = self.Maps[pickeditem.text()][4]
                self.mmcWindow.map_packet[num][2] = self.mmcWindow.map_packet[num][1].imshow(
                    self.Maps[pickeditem.text()][0], cmap='jet',
                    aspect='auto', vmin=zmin, vmax=tmax, extent=[0, x_end, 0, y_end])
            self.mmcWindow.map_packet[num][1].set_xlabel('x, mm')
            self.mmcWindow.map_packet[num][1].set_ylabel('y, mm')
            self.mmcWindow.map_packet[num][4] = self.mmcWindow.map_packet[num][0].figure.colorbar(
                self.mmcWindow.map_packet[num][2])
            self.mmcWindow.map_packet[num][1].invert_yaxis()
            self.mmcWindow.map_packet[num][1].set_aspect('equal')

    def setIM(self):
        global isIM
        isIM = True

    def setMS(self):
        global isIM
        isIM = False


if __name__ == "__main__":
    app = QApplication(sys.argv)
    mainwindow = MainGUIobject()
    # widget = mainwindow(None)
    widget = QtWidgets.QStackedWidget()
    widget.addWidget(mainwindow)
    widget.setFixedWidth(MW_width)
    widget.setFixedHeight(MW_height)
    widget.setWindowTitle("image_inspector")
    widget.show()
    sys.exit(app.exec_())

# Check if there is a ID file
# allow browsing and loading an ID file
# if there isn't an ID file, only annotate the m/z values and intensity
# do I really need to do the local max?

# New ID file may give multiple IDs, display all of them -> maybe loop them?
