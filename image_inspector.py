# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 15:58:09 2021

@author: wtd14
"""

# TODO: for 4th of July week:
# 1. Make it crash proof, so that the program will not terminate early under any circumstances.
# 2. Write good, python best practices comments for every function
# 3. Write a user manual for Image Inspector
# 4. Figure out if noise button should be implemented or not
# 7. (JC)Ask JC if the standard dev increasing when the image is flipped is okay?
# 15. Make ML heuristic scorer.
# 16. (JC)What color should the spectra plot colormap be?
# 17. (Esteban)Should the multimap maps disappear when it is destroyed?


import os
import sys

import numpy
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import (QApplication, QFileDialog, QListWidgetItem, QMainWindow)
import PyQt5.uic as uic
from PyQt5 import QtCore
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import math
import csv
from matplotlib.figure import Figure

# imported file
import ROI_class as roi  # ROI_class.py

MW_width = 1591
MW_height = 1051
isIM = None

button_style_sheet = ("QRadioButton{border:None}QRadioButton::indicator:unchecked{"
                      "border : 1px solid black;width : 25px;height : 12px;border-radius : 7px;}"
                      "QRadioButton::indicator:checked{border : 1px solid black;width : 25px;"
                      "height : 12px;border-radius : 7px;background-color : #598392}")


def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    return os.path.join(os.path.abspath("."), relative_path)


mainWindow_ui_path = resource_path("image_inspector_layout.ui")
mmcWindow_ui_path = resource_path("image_inspector_multicomp.ui")
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

        self.pickitem = None

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
        self.binI = None
        self.ROI_outline = {}
        self.ROI_Map_Data = {}
        self.pickedPointData = None
        self.mapData = None
        self.label = None
        self.scalefact = None
        self.spectra_toolbar = None
        self.annotation = None
        self.mzVals = None
        self.intensity = None
        self.drifts = None
        self.chosenDataIso = None
        self.view = None
        self.iso_view = None
        self.con_cbar = None
        self.original_image = None
        self.ROIData = None
        self.setupUi(self)
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowSystemMenuHint)

        # plots
        self._spectra_ax = None
        self._con_ax = None
        self._kin_ax = None

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
        self.zmax.sliderMoved.connect(self.zmax_Callback)
        self.zmax.valueChanged.connect(self.zmax_Callback)
        self.zmin.sliderMoved.connect(self.zmin_Callback)
        self.zmin.valueChanged.connect(self.zmin_Callback)
        self.temp_max.returnPressed.connect(self.temp_max_Callback)
        self.temp_min.returnPressed.connect(self.temp_min_Callback)
        self.max_iso.returnPressed.connect(self.temp_max_iso_Callback)
        self.min_iso.returnPressed.connect(self.temp_min_iso_Callback)
        self.zmax_isotope.sliderMoved.connect(self.zmax_isotope_Callback)
        self.zmax_isotope.valueChanged.connect(self.zmax_isotope_Callback)
        self.zmin_isotope.sliderMoved.connect(self.zmin_isotope_Callback)
        self.zmin_isotope.valueChanged.connect(self.zmin_isotope_Callback)
        self.ROI_select.clicked.connect(self.ROI_select_Callback_mask)
        self.exportROI.clicked.connect(self.exportROI_Callback)
        self.ROI_listbox.itemDoubleClicked.connect(self.ROI_listbox_Callback)
        self.importROI.clicked.connect(self.importROI_Callback)
        self.deleteROIbutton.clicked.connect(self.deleteROIbutton_Callback)
        self.clearROIbutton.clicked.connect(self.clearROI_Callback)
        self.exportROI_val.clicked.connect(self.exportROI_spectra_val_Callback)
        self.find_IDlist.clicked.connect(self.find_IDlist_Callback)
        self.Map_listbox.itemDoubleClicked.connect(self.Map_listbox_Callback)
        self.exportConcMap.clicked.connect(self.export_ConcMap_Callback)
        self.exportIsotopeMap.clicked.connect(self.export_IsotopeMap_Callback)
        self.deleteMapbutton.clicked.connect(self.deleteMapbutton_Callback)
        self.clearMapListboxbutton.clicked.connect(self.clearMapbutton_Callback)
        self.extract_Map_mzOI.clicked.connect(self.mzOI_extractMap_Callback)
        self.drift_scrollbar.sliderMoved.connect(self.drift_scrollbar_callback)
        self.drift_time.valueChanged.connect(self.drift_time_callback)
        self.one_drift_time.toggled.connect(self.button_changed_callback)
        self.all_drift_times.toggled.connect(self.button_changed_callback)
        self.reset_scatter.clicked.connect(self.reset_scatter_callback)
        self.reset_image.clicked.connect(self.reset_orig_image)
        self.set_mz_minmax.clicked.connect(self.change_mz)
        self.flipButton.clicked.connect(self.flip_figure)
        self.rotate_right_button.clicked.connect(self.rightRotate)
        self.rotate_left_button.clicked.connect(self.leftRotate)

        self.IMDataButton.clicked.connect(self.setIM)
        self.MSDataButton.clicked.connect(self.setMS)

        self.multiMW.clicked.connect(self.MultiMapCompare_Display_Callback)
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
        self.mmcWindow.importMapData.clicked.connect(self.MultiMapCompare_importMapData_Callback)

        # function image_inspector_OpeningFcn
        self.has_data = 0
        self.millimeter.setChecked(True)
        self.massplusone.setChecked(True)
        self.all_drift_times.setChecked(True)
        self.ideal_ratio.setText('1')
        QListWidgetItem('ROI list appears here', self.ROI_listbox)
        QListWidgetItem('Map list appears here', self.Map_listbox)

        self.fName = ""
        self.fName_mzOI = ""
        self.cubefilename = ""
        self.ROI = {}
        self.ROIplots = {}
        self.ROI_img_mean = {}
        self.ROI_Mass = {}
        self.ROIcount = 0
        self.Maps = {}
        self.Mapcount = 0
        self.Map_listselect_text = ""
        self.img_mean = 0

        # canvas
        self.spectra_canvas = None

        # ROI
        self.h = None
        self.ROI_listselect_text = ""
        self.ROI_listselect_array = []

        # multipmap compare variables
        self.ConcMapData = None  # map data, x_end, y_end, ()
        self.IsotopeMapData = None  # # map data, x_end, y_end, (iso_min, iso_max)

        # Style sheets
        self.micrometer.setStyleSheet(button_style_sheet)
        self.millimeter.setStyleSheet(button_style_sheet)
        self.centimeter.setStyleSheet(button_style_sheet)
        self.IMDataButton.setStyleSheet(button_style_sheet)
        self.MSDataButton.setStyleSheet(button_style_sheet)
        self.massplusone.setStyleSheet(button_style_sheet)
        self.massplustwo.setStyleSheet(button_style_sheet)
        self.one_drift_time.setStyleSheet(button_style_sheet)
        self.all_drift_times.setStyleSheet(button_style_sheet)

    def flip_figure(self):
        try:
            if self.mapData is not None:
                self.mapData.reverse()
            if self.ConcMapData is not None:
                self.ConcMapData.reverse()
                self.displayImage(self.ConcMapData, self.pixelSizeX, self.pixelSizeY)
            if self.IsotopeMapData is not None:
                self.IsotopeMapData.reverse()
                self.displayIsoImage(self.ConcMapData, self.IsotopeMapData, self.pixelSizeX, self.pixelSizeY)
        except AttributeError:
            print("Error: Flipping an array with no values.\n"
                  "Please return to an image with data if you wish to flip the image.")
            return

    def rightRotate(self):
        self.rotate(True)

    def leftRotate(self):
        self.rotate(False)

    def rotate(self, isRight):
        if self.mapData is None:
            return 0

        temp = self.pixelSizeX
        self.pixelSizeX = self.pixelSizeY
        self.pixelSizeY = temp

        if isRight:
            self.mapData = self.rotateRight(self.mapData)
        else:
            self.mapData = self.rotateLeft(self.mapData)

        if self.ConcMapData:
            if isRight:
                self.ConcMapData = self.rotateRight(self.ConcMapData)
            else:
                self.ConcMapData = self.rotateLeft(self.ConcMapData)
            self.displayImage(self.ConcMapData, self.pixelSizeX, self.pixelSizeY)
        if self.IsotopeMapData:
            if isRight:
                self.IsotopeMapData = self.rotateRight(self.IsotopeMapData)
            else:
                self.IsotopeMapData = self.rotateLeft(self.IsotopeMapData)
            self.displayIsoImage(self.ConcMapData, self.IsotopeMapData, self.pixelSizeX, self.pixelSizeY)
        return 0

    def rotateRight(self, origMap):
        newMap = []
        for i in range(len(origMap[0])):
            newLine = []
            for j in range(len(origMap) - 1, -1, -1):
                newLine.append(origMap[j][i])
            newMap.append(newLine)
        return newMap

    def rotateLeft(self, origMap):
        newMap = []
        for i in range(len(origMap[0]) - 1, -1, -1):
            newLine = []
            for j in range(len(origMap)):
                newLine.append(origMap[j][i])
            newMap.append(newLine)
        return newMap

    def reset_orig_image(self):
        if self.original_image:
            self.pickedPointData = None
            self.displayImage(self.original_image, self.pixelSizeX, self.pixelSizeY)
            self.start.clear()
            self.massbox.clear()
            while self.plot_kin.count():
                child = self.plot_kin.takeAt(0)
                if child.widget():
                    child.widget().deleteLater()

    def reset_scatter_callback(self):
        if isIM:
            self.one_drift_time.setChecked(False)
            self.all_drift_times.setChecked(True)
        if self.mzVals is None:
            print("Error: There is no original plot. Please select a .bin file and press 'GO' ")
            return
        self.displayScatter(self.mzVals, self.intensity, self.drifts)
        self.set_min_max_mz(self.mzVals)

    def change_mz(self):
        max_mz = self.max_mz.value()
        min_mz = self.min_mz.value()
        mzVals = self.mzVals
        intensities = self.intensity

        if mzVals is None:
            return 0

        processed_mz = []
        processed_intens = []

        if isIM:
            drifts = self.drifts
            processed_drift = []
            for i in range(len(mzVals)):
                if max_mz >= mzVals[i] >= min_mz:
                    processed_mz.append(mzVals[i])
                    processed_intens.append(intensities[i])
                    processed_drift.append(drifts[i])
            self.displayScatter(processed_mz, processed_intens, processed_drift)
        else:
            for i in range(len(mzVals)):
                if max_mz >= mzVals[i] >= min_mz:
                    processed_mz.append(mzVals[i])
                    processed_intens.append(intensities[i])
            self.displayScatter(processed_mz, processed_intens, None)
        try:
            self.set_min_max_mz(processed_mz)
        except ValueError:
            print("Error: There is no spectra to plot")
            return

    # --- Executes on button press in find_IDlist.
    # can load a new ID list the consists of two columns
    # must be .csv file
    # the 1st column must be named "m/z"
    # the 2nd column must be named "Lipid ID"
    def find_IDlist_Callback(self):
        # But this code does not handle .h5 or .mat files
        self.fName_IDlist = QFileDialog.getOpenFileName(self, 'Pick ID List', filter='*.csv')
        self.IDlist_name.setText(self.fName_IDlist[0])
        if self.fName_IDlist[0] == '':
            print("Please select a file to import as the ID list.")
            return
        self.ids_pd = pd.read_csv(self.fName_IDlist[0])

    def button_changed_callback(self):
        if not self.drifts:
            return 0
        if self.one_drift_time.isChecked():
            val = self.drift_time.value()
            self.show_mz_map(val)
        else:
            self.change_mz()
        return 0

    def drift_time_callback(self):
        val = self.drift_time.value()
        self.drift_scrollbar.setValue(val)
        if self.drifts is None:
            return
        if self.one_drift_time.isChecked():
            self.show_mz_map(val)

    def drift_scrollbar_callback(self):
        val = self.drift_scrollbar.sliderPosition()
        self.drift_time.setValue(val)
        if self.drifts is None:
            return
        if self.one_drift_time.isChecked():
            self.show_mz_map(val)

    def show_mz_map(self, val):
        theDrifts = np.asarray(self.drifts)
        drift_vals = np.where(val == theDrifts, self.drifts, 0)
        drift_vals = drift_vals.nonzero()

        mzVals = np.asarray(self.mzVals)[drift_vals]
        intensity = np.asarray(self.intensity)[drift_vals]

        theMax = self.max_mz.value()
        theMin = self.min_mz.value()

        in1 = np.where(theMax >= mzVals, mzVals, 0)
        in1 = in1.nonzero()

        mzVals = np.asarray(mzVals)[in1]
        intensity = np.asarray(intensity)[in1]

        in1 = np.where(mzVals >= theMin, mzVals, 0)
        in1 = in1.nonzero()

        mzVals = np.asarray(mzVals)[in1]
        intensity = np.asarray(intensity)[in1]

        self.displayScatter(mzVals, intensity, None, .3)
        return 0

    def mass_up_Callback(self):
        if not self.view:
            return 0
        if isIM:
            val = float(self.start.text()) + 1.0
            self.start.setText(str(val))
            self.im_point()
            return 0

    # --- Executes on button press in mass_down.
    # When the button is clicked, the mass index is decremented by 1
    # and the image for the new mass is displayed
    def mass_down_Callback(self):
        if not self.view:
            return 0
        if isIM:
            val = float(self.start.text()) - 1.0
            self.start.setText(str(val))
            self.im_point()
            return 0

    def temp_max_Callback(self):
        newMax = self.plot_con_temp_pressed_callback(self.temp_max.text())
        self.zmax.setValue(math.ceil(newMax))
        self.scale_image()

    def temp_min_Callback(self):
        newMin = self.plot_con_temp_pressed_callback(self.temp_min.text())
        self.zmin.setValue(math.ceil(newMin))
        self.scale_image()

    def plot_con_temp_pressed_callback(self, temp):
        try:
            newVal = float(temp)
        except ValueError:
            print("Please enter a numeric value.")
            newVal = float(self.min_int.text())
            return newVal
        if newVal > float(self.max_int.text()):
            print("Error: Trying to set the value above the max")
            newVal = float(self.max_int.text())
        elif newVal < float(self.min_int.text()):
            print("Error: Trying to set the value below the min")
            newVal = float(self.min_int.text())
        return newVal

    def temp_max_iso_Callback(self):
        newMax = self.plot_kin_temp_pressed_callback(self.max_iso.text())
        self.zmax_isotope.setValue(math.ceil(newMax))
        self.scale_iso_image()

    def temp_min_iso_Callback(self):
        newMin = self.plot_kin_temp_pressed_callback(self.min_iso.text())
        self.zmin_isotope.setValue(math.ceil(newMin))
        self.scale_iso_image()

    def plot_kin_temp_pressed_callback(self, temp):
        try:
            newVal = float(temp)
        except ValueError:
            print("Please enter a numeric value.")
            newVal = float(self.min_int_iso.text())
            return newVal
        if newVal > float(self.max_int_iso.text()):
            print("Error: Trying to set the value above the max")
            newVal = float(self.max_int_iso.text())
        elif newVal < float(self.min_int_iso.text()):
            print("Error: Trying to set the value below the min")
            newVal = float(self.min_int_iso.text())
        return newVal

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
        if self.start.text() == '':
            print("You must choose or input a point to the 'Selected Mass' box first.")
            return
        if self.massbox.text() == '':
            print("There is no plot to select")
            return
        else:
            if self.h:
                self.h.disconnect()
            self.h = roi.new_ROI(self.con_img)

    # This function plots a mass spectrum corresponing to a selected point
    # on the displayed image, or an image corresponding to a selected point
    # on the mass spectrum
    # --- Executes on button press in pick_point.
    def pick_point_Callback(self):
        if self.micrometer.isChecked():
            self.scalefact = 1e3
            self.label = 'μm'
        elif self.millimeter.isChecked():
            self.scalefact = 1
            self.label = 'mm'
        elif self.centimeter.isChecked():
            self.scalefact = 0.1
            self.label = 'cm'
        if self.start.text() == '':
            print("You must choose or input a point to the 'Selected Mass' box first.")
            return 0
        elif isIM:
            self.im_point()
            return 0
        else:
            self.ms_point()
            return 0

    def im_point(self):
        if self.start.text() == '':
            return
        try:
            picked_point = float(self.start.text())
        except ValueError:
            print("Please enter a number for the selected mass.")
            return
        max_diff = self.ppm_calc(picked_point)
        try:
            ideal_ratio = float(self.ideal_ratio.text())
        except ValueError:
            print("Please enter a number above 0 for the m/z spacing.")
            return
        if ideal_ratio <= 0:
            print("Please enter a spacing greater than 0")
            return
        self.massbox.setText(self.start.text())
        mapData = self.mapData

        chosenData = []
        chosenDataPlusOne = []
        chosenDataPlusTwo = []
        maxIntensity = 0
        maxIntensityPlusOne = 0
        maxIntensityPlusTwo = 0

        for line in mapData:
            lineData = []
            lineDataPlusOne = []
            lineDataPlusTwo = []
            for frame in line:
                theVal = 0
                theValPlusOne = 0
                theValPlusTwo = 0
                for scan in frame:
                    if picked_point - max_diff <= scan[0] < picked_point + max_diff:
                        theVal += scan[1]
                        if theVal > maxIntensity:
                            maxIntensity = theVal
                    elif picked_point + (1 / ideal_ratio) - max_diff <= scan[0] < picked_point + (
                            1 / ideal_ratio) + max_diff:
                        theValPlusOne += scan[1]
                        if theValPlusOne > maxIntensityPlusOne:
                            maxIntensityPlusOne = theValPlusOne
                    elif picked_point + (2 / ideal_ratio) - max_diff <= scan[0] < picked_point + (
                            2 / ideal_ratio) + max_diff:
                        theValPlusTwo += scan[1]
                        if theValPlusTwo > maxIntensityPlusTwo:
                            maxIntensityPlusTwo = theValPlusTwo
                lineData.append(theVal)
                lineDataPlusOne.append(theValPlusOne)
                lineDataPlusTwo.append(theValPlusTwo)
            chosenData.append(lineData)
            chosenDataPlusOne.append(lineDataPlusOne)
            chosenDataPlusTwo.append(lineDataPlusTwo)

        m_zero_sum = 0
        m_one_sum = 0
        m_two_sum = 0
        if maxIntensity != 0:
            denom = maxIntensity + maxIntensityPlusOne + maxIntensityPlusTwo
            m_zero_sum = round(maxIntensity / denom, 4)
            m_one_sum = round(maxIntensityPlusOne / denom, 4)
            m_two_sum = round(maxIntensityPlusTwo / denom, 4)

        self.Msumratio.setText(str(m_zero_sum))
        self.Mplusonesumratio.setText(str(m_one_sum))
        self.Mplustwosumratio.setText(str(m_two_sum))

        self.pickedPointData = None
        self.displayImage(chosenData, self.pixelSizeX, self.pixelSizeY)
        self.pickedPointData = chosenData

        if self.massplusone.isChecked():
            self.chosenDataIso = None
            self.displayIsoImage(chosenData, chosenDataPlusOne, self.pixelSizeX, self.pixelSizeY)
            self.chosenDataIso = chosenDataPlusOne
        elif self.massplustwo.isChecked():
            self.chosenDataIso = None
            self.displayIsoImage(chosenData, chosenDataPlusTwo, self.pixelSizeX, self.pixelSizeY)
            self.chosenDataIso = chosenDataPlusTwo

    def ms_point(self):
        if self.start.text() == '':
            return
        try:
            picked_point = float(self.start.text())
        except ValueError:
            print("Please enter a number for the selected mass.")
            return
        max_diff = self.ppm_calc(picked_point)
        try:
            ideal_ratio = float(self.ideal_ratio.text())
        except ValueError:
            print("Please enter a number above 0 for the m/z spacing.")
            return
        if ideal_ratio <= 0:
            print("Please enter a spacing greater than 0")
            return

        self.massbox.setText(self.start.text())

        mapData = self.mapData

        imageData = []
        image_plus_one = []
        image_plus_two = []

        for line in mapData:
            lineVals = []
            line_plus_one = []
            line_plus_two = []
            for scan in line:
                intensity = 0
                intensity_plus_one = 0
                intensity_plus_two = 0
                for val in scan:
                    if val[0] + max_diff >= picked_point >= val[0] - max_diff:
                        intensity += val[1]
                    elif val[0] + (1 / ideal_ratio) + max_diff >= picked_point >= val[0] + (1 / ideal_ratio) - max_diff:
                        intensity_plus_one += val[1]
                    elif val[0] + (2 / ideal_ratio) + max_diff >= picked_point >= val[0] + (2 / ideal_ratio) - max_diff:
                        intensity_plus_two += val[1]
                lineVals.append(intensity)
                line_plus_one.append(intensity_plus_one)
                line_plus_two.append(intensity_plus_two)
            imageData.append(lineVals)
            image_plus_one.append(line_plus_one)
            image_plus_two.append(line_plus_two)

        self.pickedPointData = None
        self.displayImage(imageData, self.pixelSizeX, self.pixelSizeY)
        self.pickedPointData = imageData

        not_used, m_zero_max_intensity = self.find_min_max_image(imageData)
        not_used, m_one_max_intensity = self.find_min_max_image(image_plus_one)
        not_used, m_two_max_intensity = self.find_min_max_image(image_plus_two)
        if m_zero_max_intensity != 0:
            denom = m_zero_max_intensity + m_one_max_intensity + m_two_max_intensity
            self.Msumratio.setText(str(round(m_zero_max_intensity / denom, 4)))
            self.Mplusonesumratio.setText(str(round(m_one_max_intensity / denom, 4)))
            self.Mplustwosumratio.setText(str(round(m_two_max_intensity / denom, 4)))

        if self.massplusone.isChecked():
            self.chosenDataIso = None
            self.displayIsoImage(imageData, image_plus_one, self.pixelSizeX, self.pixelSizeY)
            self.chosenDataIso = image_plus_one
        if self.massplustwo.isChecked():
            self.chosenDataIso = None
            self.displayIsoImage(imageData, image_plus_two, self.pixelSizeX, self.pixelSizeY)
            self.chosenDataIso = image_plus_two
        return 0

    def isotope_scalar(self, m_zero_intensity, isotope_intensity):
        new_data = []
        for i in range(len(m_zero_intensity)):
            theLine = []
            orig_line = m_zero_intensity[i]
            iso_line = isotope_intensity[i]
            for j in range(len(orig_line)):
                if orig_line[j] == 0:
                    ratio = 0
                else:
                    ratio = iso_line[j] / orig_line[j] + iso_line[j]
                theLine.append(ratio)
            new_data.append(theLine)
        return new_data

    # --- Executes on button press in find_file.
    def find_file_Callback(self):
        # But this code does not handle .h5 or .mat files
        self.fName = QFileDialog.getOpenFileName(self, 'Pick Data Cube', filter='*.mat, *.h5, *.bin *.csv')
        self.wspc_name.setText(self.fName[0])

    # --- Executes on button press in start_cube.
    def start_cube_Callback(self):
        self.cubefilename = self.fName[0]
        filename = self.cubefilename
        self.functionsCommonToAll()
        print("Working to read datacube")
        if filename.endswith('.mat'):
            print("This code can't process .mat files")
            return
        elif filename.endswith('.h5'):
            print("This code can't process .h5 files")
            return
        elif filename.endswith('.csv'):
            data = ((pd.read_csv(self.cubefilename, header=None)).to_numpy(numpy.float32)).flatten()
        elif filename.endswith('.bin'):
            file = open(self.cubefilename)
            data = np.fromfile(file, dtype=np.float32)
            file.close()
        else:
            print('Unexpected file extension')
            return
        try:
            if isIM:
                self.cubeAsIMData(data)
            elif self.MSDataButton.isChecked():
                self.cubeAsMSData(data)
            else:
                print("Please select whether the file is IM or MS Data")
        except IndexError:
            if isIM:
                print("Did you mean to select MS Data?")
            else:
                print("Did you mean to select IM Data?")

    def functionsCommonToAll(self):
        if self.micrometer.isChecked():
            self.scalefact = 1e3
            self.label = 'μm'
        elif self.millimeter.isChecked():
            self.scalefact = 1
            self.label = 'mm'
        elif self.centimeter.isChecked():
            self.scalefact = 0.1
            self.label = 'cm'

        while self.plot_kin.count():
            child = self.plot_kin.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        self.iso_view = None

        self.chosenDataIso = None

        self.max_int.clear()
        self.min_int.clear()
        self.max_int_iso.clear()
        self.min_int_iso.clear()
        self.max_iso.clear()
        self.min_iso.clear()
        self.temp_max.clear()
        self.temp_min.clear()
        self.Msumratio.clear()
        self.Msumstandard_error.clear()
        self.numberpoints.clear()
        self.Mplusonesumratio.clear()
        self.Mplusonesumstandard_error.clear()
        self.Mplustwosumratio.clear()
        self.Mplustwosumstandard_error.clear()
        self.start.clear()
        self.Noise_Output_Box.clear()
        self.exportConcMapname.setText("Conc. Map Name")
        self.exportIsotopeMapname.setText("Isotope Map Name")
        self.ID_Output_Box.clear()
        self.massbox.clear()

        self.ROIplots.clear()
        self.ROI.clear()
        self.ROI_img_mean.clear()
        self.ROIcount = 0
        self.ROIcountbox.setText("0")
        self.ROI_listbox.clear()
        QListWidgetItem('ROI list appears here', self.ROI_listbox)

        self.clearMapbutton_Callback()

        self.pickedPointData = None
        self.pick_IDthreshold.setValue(20)
        self.pick_IDthreshold.setMaximum(1000)

        self.zmax_isotope.setMaximum(99)
        self.zmax_isotope.setMinimum(0)
        self.zmin_isotope.setMaximum(99)
        self.zmin_isotope.setMinimum(0)
        self.mzVals = None
        self.drifts = None
        self.intensity = None

    def displayIsoImage(self, zero_image, imageData, pixelSizeX, pixelSizeY):
        # This needs to be a different function for both the original point and the scaling functions
        count = len(zero_image) * len(zero_image[0])
        self.numberpoints.setText(str(count))

        xend = len(imageData[0]) * (pixelSizeX / 1000)
        yend = len(imageData) * (pixelSizeY / 1000)

        while self.plot_kin.count():
            child = self.plot_kin.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

        zero_for_dev = np.asarray(zero_image).flatten()
        std_deviation = round(float(np.std(zero_for_dev)), 4)
        self.Msumstandard_error.setText(str(std_deviation))

        iso_data = self.isotope_scalar(zero_image, imageData)
        self.IsotopeMapData = iso_data
        iso_for_deviation = np.asarray(iso_data)
        iso_for_deviation = iso_for_deviation.flatten()

        std_deviation = round(float(np.std(iso_for_deviation)), 4)

        if self.massplusone.isChecked():
            self.Mplusonesumstandard_error.setText(str(std_deviation))
            self.Mplustwosumstandard_error.clear()
        elif self.massplustwo.isChecked():
            self.Mplustwosumstandard_error.setText(str(std_deviation))
            self.Mplusonesumstandard_error.clear()
        else:
            print("Please choose whether to display mass plus one or mass plus two.")
            return

        self.iso_view = FigureCanvas(Figure(figsize=(5, 3)))
        self.axes = self.iso_view.figure.subplots()
        self.toolbar = NavigationToolbar(self.iso_view, self)
        self.plot_kin.addWidget(self.iso_view)
        self.con_img2 = self.axes.imshow(iso_data, cmap='inferno', aspect=(yend / xend),
                                         extent=[0, xend * self.scalefact, 0, yend * self.scalefact])
        self.axes.set_xlabel("x, " + self.label)
        self.axes.set_ylabel("y, " + self.label)
        plt.colorbar(self.con_img2)
        self.iso_view.draw()
        if not self.chosenDataIso:
            theMin, theMax = self.find_min_max_image(imageData)
            # If this needs to be the scaled image just change the above line
            self.max_int_iso.setText(str(theMax))
            self.min_int_iso.setText(str(theMin))
            self.max_iso.setText(str(theMax))
            self.min_iso.setText(str(theMin))
            self.zmax_isotope.setMaximum(math.ceil(theMax))
            self.zmax_isotope.setMinimum(math.floor(theMin))
            self.zmin_isotope.setMaximum(math.ceil(theMax))
            self.zmin_isotope.setMinimum(math.floor(theMin))
            self.zmax_isotope.setValue(math.ceil(theMax))
            self.zmin_isotope.setValue(math.floor(theMin))

    def displayImage(self, imageData, pixelSizeX, pixelSizeY):
        xend = len(imageData[0]) * (pixelSizeX / 1000)
        yend = len(imageData) * (pixelSizeY / 1000)

        if self.view:
            self.plot_con.removeWidget(self.view)
        self.view = FigureCanvas(Figure(figsize=(5, 3)))
        self.axes = self.view.figure.subplots()
        self.toolbar = NavigationToolbar(self.view, self)
        self.plot_con.addWidget(self.view)
        self.con_img = self.axes.imshow(imageData, cmap='jet', aspect=(yend / xend),
                                        extent=[0, xend * self.scalefact, 0, yend * self.scalefact])
        plt.colorbar(self.con_img)
        self.axes.set_title('Points In Selected Region')
        self.axes.set_xlabel("x, " + self.label)
        self.axes.set_ylabel("y, " + self.label)
        self.view.draw()
        self.ConcMapData = imageData
        if not self.pickedPointData:
            theMin, theMax = self.find_min_max_image(imageData)
            self.max_int.setText(str(theMax))
            self.min_int.setText(str(theMin))
            self.temp_max.setText(str(theMax))
            self.temp_min.setText(str(theMin))
            self.zmax.setMaximum(math.ceil(theMax))
            self.zmax.setMinimum(math.floor(theMin))
            self.zmin.setMaximum(math.ceil(theMax))
            self.zmin.setMinimum(math.floor(theMin))
            self.zmax.setValue(math.ceil(theMax))
            self.zmin.setValue(math.floor(theMin))
            self.pickedPointData = imageData

    def find_min_max_image(self, imageData):
        theMin = sys.maxsize
        theMax = 0
        for line in imageData:
            lineMin = min(line)
            lineMax = max(line)
            if lineMax > theMax:
                theMax = lineMax
            if lineMin < theMin:
                theMin = lineMin
        return theMin, theMax

    def displayScatter(self, mzVals, intensity, drifts=None, pt_size=.01):
        while self.plot_spectra.count():  # This solved the double figure problem
            child = self.plot_spectra.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        plt.close('all')

        self.spectra_canvas = FigureCanvas(plt.figure(tight_layout=True))
        self.spectra_canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.spectra_canvas.setFocus()
        self.spectra_toolbar = NavigationToolbar(self.spectra_canvas, self)
        self.plot_spectra.addWidget(self.spectra_toolbar)
        self.plot_spectra.addWidget(self.spectra_canvas)
        self._spectra_ax = self.spectra_canvas.figure.subplots()
        if drifts is not None:
            x = self._spectra_ax.scatter(mzVals, intensity, s=pt_size, c=drifts, cmap="turbo", alpha=0.75, picker=True)
            plt.colorbar(x).set_label('Drift times')
            try:
                self.drift_scrollbar.setMinimum(int(min(drifts)))
                self.drift_scrollbar.setMaximum(int(max(drifts)))
                self.drift_time.setMinimum(int(min(drifts)))
                self.drift_time.setMaximum(int(max(drifts)))
            except ValueError:
                print("There are no values in the spectra")
                self.drift_scrollbar.setMinimum(0)
                self.drift_scrollbar.setMaximum(0)
                self.drift_time.setMinimum(0)
                self.drift_time.setMaximum(0)

        else:
            self._spectra_ax.scatter(mzVals, intensity, color='#393424', s=pt_size, alpha=0.75, picker=True)  # Can change to peaks here
        self._spectra_ax.set_title('Points In Selected Region')
        self._spectra_ax.set_xlabel('m/z')
        self._spectra_ax.set_ylabel('intensity')
        self.spectra_canvas.mpl_connect('pick_event', self.data_cursor_click)
        plt.ylabel('intensity')
        plt.xlabel('m/z')

    def set_min_max_mz(self, mzVals):
        try:
            self.min_mz.setRange(min(mzVals), max(mzVals))
            self.max_mz.setRange(min(mzVals), max(mzVals))
            self.min_mz.setValue(min(mzVals))
            self.max_mz.setValue(max(mzVals))
        except TypeError:
            print("Error: There is no plot. Please select a .bin file and press 'GO' ")
            return

    def cubeAsMSData(self, data):
        numLines = data[0]
        numScans = data[1]
        self.pixelSizeX = data[2]
        self.pixelSizeY = data[3]

        lineNum = 0
        i = 4

        mzVals = []
        intensities = []
        imageData = []
        mapData = []

        while lineNum < numLines:
            scanNum = 0
            line = []
            mapLine = []
            while scanNum < numScans:
                scan = 0
                numDataPoints = data[i]
                i += 1
                currDataPoint = 0
                scanVal = []
                while currDataPoint < numDataPoints:
                    currDataPoint += 1
                    mzVals.append(data[i])
                    intensities.append(data[i + 1])
                    scan += data[i + 1]
                    scanVal.append([data[i], data[i + 1]])
                    i += 2
                scanNum += 1
                line.append(scan)
                mapLine.append(scanVal)
            lineNum += 1
            imageData.append(line)
            mapData.append(mapLine)
        self.displayImage(imageData, self.pixelSizeX, self.pixelSizeY)
        self.displayScatter(mzVals, intensities, None)
        self.set_min_max_mz(mzVals)
        self.mapData = mapData
        self.mzVals = mzVals
        self.intensity = intensities
        self.original_image = imageData

    def cubeAsIMData(self, data):
        frameNum = data[0]
        fileNum = data[1]
        self.pixelSizeX = data[2]
        self.pixelSizeY = data[3]

        numFrames = 0
        numFiles = 0
        i = 4

        mzVals = []
        intensity = []
        drifts = []
        chosenData = []
        mapData = []
        lineData = []
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
            mapLine = []
            while numFrames < frameNum:
                totalDriftBins = data[i]
                frameDone = True
                currdriftBin = 0
                i += 1
                frameVal = []
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
                        frameVal.append([data[i], data[i + 1], driftTime])
                        i += 2
                    currdriftBin += 1
                mapLine.append(frameVal)
                if not valAdded:
                    lineData.append(0)
                elif valAdded:
                    lineData.append(theVal)
                valAdded = False
                theVal = 0
                numFrames += 1
            mapData.append(mapLine)

        self.displayScatter(mzVals, intensity, drifts)
        self.set_min_max_mz(mzVals)

        self.mzVals = mzVals
        self.intensity = intensity
        self.drifts = drifts
        self.mapData = mapData

        self.displayImage(chosenData, self.pixelSizeX, self.pixelSizeY)

        self.has_data = 1
        self.original_image = chosenData

    # --- Executes on slider movement.
    def zmax_Callback(self):
        self.temp_max.setText(str(self.zmax.sliderPosition()))
        self.scale_image()

    def zmin_Callback(self):
        self.temp_min.setText(str(self.zmin.sliderPosition()))
        self.scale_image()

    # This function updates the image to the current index
    def scale_image(self):
        if not self.pickedPointData:
            return
        x = self.pickedPointData
        data = []
        maximum = self.zmax.sliderPosition()
        minimum = self.zmin.sliderPosition()
        if minimum >= maximum:
            newMap = np.zeros((len(x), len(x[0])))
            self.displayImage(newMap, self.pixelSizeX, self.pixelSizeY)
            return 0
        if isIM:
            for line in x:
                newLine = []
                for frame in line:
                    newFrame = 0  # Changed this from a 0 to the minimum. Is it right??
                    if frame > maximum:
                        newFrame = maximum
                    elif maximum >= frame >= minimum:
                        newFrame = frame
                    newLine.append(newFrame)
                data.append(newLine)
            self.displayImage(data, self.pixelSizeX, self.pixelSizeY)
        else:
            for line in x:
                newLine = []
                for scan in line:  # Same as above line
                    newScan = 0
                    if scan > maximum:
                        newFrame = maximum
                    elif maximum >= scan >= minimum:
                        newScan = scan
                    newLine.append(newScan)
                data.append(newLine)
            self.displayImage(data, self.pixelSizeX, self.pixelSizeY)

    # --- Executes on slider movement.
    def zmax_isotope_Callback(self):
        self.max_iso.setText(str(self.zmax_isotope.sliderPosition()))
        self.scale_iso_image()

    # --- Executes on slider movement.
    def zmin_isotope_Callback(self):
        self.min_iso.setText(str(self.zmin_isotope.sliderPosition()))
        self.scale_iso_image()

    # This function updates the image to the current index
    # Figure out what this function is trying to do!
    def scale_iso_image(self):
        if self.chosenDataIso is None:  # This will be triggered only at the beginning
            return 0
        x = self.chosenDataIso
        data = []
        maximum = self.zmax_isotope.sliderPosition()
        minimum = self.zmin_isotope.sliderPosition()
        if minimum >= maximum:
            newMap = np.zeros((len(x), len(x[0])))
            self.displayIsoImage(newMap, newMap, self.pixelSizeX, self.pixelSizeY)
            return 0
        for line in x:
            newLine = []
            for frame in line:
                newFrame = 0
                if frame > maximum:
                    newFrame = maximum
                elif maximum >= frame >= minimum:
                    newFrame = frame
                newLine.append(newFrame)
            data.append(newLine)
        self.displayIsoImage(self.pickedPointData, data, self.pixelSizeX, self.pixelSizeY)

    def export_ConcMap_Callback(self):
        if self.ConcMapData is None:
            print("There is no image to export")
            return
        # when clicking on the exportConcMap button, it will save the filename and concentration of the map
        self.Maps[self.exportConcMapname.text()] = self.ConcMapData
        self.refreshMaplistbox()
        self.Mapcount += 1

    def export_IsotopeMap_Callback(self):
        if self.IsotopeMapData is None:
            print("There is no image to export")
            return
        # when clicking on the exportIsotopeMap button, it will save the filename and concentration of the map
        self.Maps[self.exportIsotopeMapname.text()] = self.IsotopeMapData
        self.Mapcount += 1
        self.refreshMaplistbox()

    def Map_listbox_Callback(self, item):
        self.Map_listselect_text = item.text()

    def deleteMapbutton_Callback(self):
        if self.Mapcount == 0:
            return
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
            return
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
        if self.h is None:
            print("Please choose an ROI before processing")
            return
        if self.view:
            self.binI = self.h.get_mask().astype(int)
            self.binI = np.flipud(self.binI)
            f = np.argwhere(np.ravel(self.binI, order='C'))[:, 0]
            self.ROI_outline[self.exportROIfilename.text()] = f
            self.numberpoints.setText(str(len(f)))
        self.ROIcount = self.ROIcount + 1
        self.ROIcountbox.setText(str(self.ROIcount))
        self.ROI[self.exportROIfilename.text()] = self.binI
        self.ROI_Mass[self.exportROIfilename.text()] = float(self.massbox.text())
        self.ROI_Map_Data[self.exportROIfilename.text()] = self.mapData
        self.refreshROIlistbox()
        return 0

    def refreshROIlistbox(self):
        if len(self.ROI) == 0:
            # set box with default text
            self.ROI_listbox.clear()
            # del self.ROI_listselect_text
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
        try:
            self.ROI_listselect_array = self.ROI[item.text()]
        except KeyError:
            print("This is not a ROI. Please choose a ROI. ")
            return
        self.ROI_listselect_text = item.text()
        self.ROI_img_mean[item.text()] = self.img_mean
        self._con_ax = 1
        self.ROIplots[self.ROI_listselect_text] = item  # What am I doing here??
        # self.ROI[self.ROI_listselect_text] =
        return 0

    # loads the spectra from the selected ROI in the ROI_listbox
    # --- Executes on button press in 'Load selected ROI spectra' button
    def importROI_Callback(self):
        if (self.ROI_listselect_text == ""):
            print("No item selected")
        else:
            try:
                chosenVal = self.ROI_Mass[self.ROI_listselect_text]
            except ValueError:
                print("Error: The chosen m/z must be numeric")
                return
            ppm = self.ppm_calc(chosenVal)

            self.ROIData = []
            mzVals = []
            intensities = []
            drifts = None

            outline = self.ROI_outline[self.ROI_listselect_text]
            mapData = self.ROI_Map_Data[self.ROI_listselect_text]
            if self.drifts is None:
                counter = 0
                for line in mapData:
                    for scan in line:
                        if counter not in outline:
                            counter += 1
                            continue
                        else:
                            for point in scan:
                                mz = point[0]
                                if chosenVal + ppm >= mz >= chosenVal - ppm:
                                    mzVals.append(mz)
                                    intensities.append(point[1])
                            counter += 1
                for i in range(len(mzVals)):
                    self.ROIData.append([mzVals[i], intensities[i]])
            else:
                drifts = []
                counter = 0
                for line in mapData:
                    for frame in line:
                        if counter not in outline:
                            counter += 1
                            continue
                        else:
                            for point in frame:
                                mz = point[0]
                                if chosenVal + ppm >= mz >= chosenVal - ppm:
                                    mzVals.append(mz)
                                    intensities.append(point[1])
                                    drifts.append(point[2])
                            counter += 1
                for i in range(len(mzVals)):
                    self.ROIData.append([mzVals[i], intensities[i], drifts[i]])

            self.displayScatter(mzVals, intensities, drifts, 100 / len(mzVals))
            self.set_min_max_mz(mzVals)

    def exportROI_spectra_val_Callback(self):
        if (self.ROI_listselect_text == ""):
            print("No item selected. Double click an item in the ROI listbox")
        else:
            path = QFileDialog.getSaveFileName(self, 'Save CSV', os.getenv('HOME'), 'CSV(*.csv)')
            if path[0] != '':
                with open(path[0], 'w', newline='') as csv_file:
                    writer = csv.writer(csv_file, dialect='excel', delimiter=',', lineterminator='\n')
                    if isIM:
                        writer.writerow(["M/Z Value", "Intensity", "Drift Time"])
                        for line in self.ROIData:
                            writer.writerow(line)
                    else:
                        writer.writerow(["M/Z Value", "Intensity"])
                        for line in self.ROIData:
                            writer.writerow(line)
            else:
                print("No location selected. Please select a location")

    def deleteROIbutton_Callback(self):
        if self.ROIcount == 0:
            self.ROI_listselect_text = ""
            return
        else:
            if self._con_ax:
                if self.ROI_listselect_text in self.ROIplots:
                    del self.ROIplots[self.ROI_listselect_text]
                    del self.ROI[self.ROI_listselect_text]
                    del self.ROI_img_mean[self.ROI_listselect_text]
                    self.ROIcount = self.ROIcount - 1
                    self.ROIcountbox.setText(str(self.ROIcount))
                    self.refreshROIlistbox()
                    print('plot removed')
                    if isIM:
                        self.im_point()
                    else:
                        self.ms_point()
                else:
                    print('Item does not exist. Please double click on another item')
        self.ROI_listselect_text = ""

    def clearROI_Callback(self):
        self.ROIplots.clear()
        self.ROI.clear()
        self.ROI_img_mean.clear()
        self.ROIcount = 0
        self.ROIcountbox.setText("0")
        self.refreshROIlistbox()
        self.reset_scatter_callback()
        self.ROI_listselect_text = ""
        if isIM:
            self.im_point()
        else:
            self.ms_point()

    # --- Executes on button press in find_file_mzOI.
    def find_file_mzOI_Callback(self):
        self.fName_mzOI = QFileDialog.getOpenFileName(self, 'Pick list: m/z of interest', filter='*.csv')
        self.mzOI_listname.setText(self.fName_mzOI[0])
        # self.fName_mzOI_flag = True

    def mzOI_extractMap_Callback(self):
        if self.mzOI_listname.text() == '':
            print("Please choose a .csv file with m/z values")
            return 0
        if self.mzVals is None:
            print("Please choose a .bin file to process")
            return 0
        mzOI = []
        try:
            mzOI = pd.read_csv(self.mzOI_listname.text())
        except Exception:
            print("Please choose a .csv file with m/z values")
            return 0
        mzOI = mzOI.to_numpy()

        mzVals = []
        intensity = []
        drifts = None
        if isIM:
            drifts = []

        for mz in mzOI:
            diff = self.ppm_calc(mz)
            for i in range(len(self.mzVals)):
                if mz + diff >= self.mzVals[i] >= mz - diff:
                    mzVals.append(self.mzVals[i])
                    intensity.append(self.intensity[i])
                    if isIM:
                        drifts.append(self.drifts[i])

        self.displayScatter(mzVals, intensity, drifts)

    def data_cursor_click(self, event):
        # indexes = event.ind
        # i = indexes[0]
        theXmOverZ = event.mouseevent.lastevent.xdata
        theY = event.mouseevent.lastevent.ydata

        self.start.setText("%.5f" % theXmOverZ)
        self.spectra_annotation(theXmOverZ, theY)

    # This is a function that calculates how far away a value must be to be defined as a different lipid.
    def ppm_calc(self, mzVal):
        if float(self.pick_IDthreshold.value()) == 0:
            return 0  # is this correct??
        return mzVal * float(self.pick_IDthreshold.value()) / 1e6

    def spectra_annotation(self, mz, intensity):
        diff = self.ppm_calc(mz)
        lipid_map = {}
        lipid_id = "not defined"
        if self.ids_pd is not None:
            mz_vals = self.ids_pd['m/z']
            lipid_ids = self.ids_pd['Lipid ID']
            for i in range(mz_vals.size):
                x = mz_vals[i]
                if (mz + diff) > x > (mz - diff):
                    difference = abs(mz - x)
                    lipid_map[difference] = lipid_ids[i]
        if len(lipid_map.keys()) > 0:
            key = min(lipid_map.keys())
            lipid_id = lipid_map[key]

        self.ID_Output_Box.setText(lipid_id)

        if self.annotation is not None:
            self.annotation.remove()

        if isIM:
            mz_vals = np.asarray(self.mzVals)
            drifts = np.asarray(self.drifts)
            vals = np.where(mz + diff > mz_vals, mz_vals, 0)
            vals = np.where(mz - diff < vals, drifts, 0)

            vals = vals.nonzero()
            finals = drifts[vals]
            if len(finals) != 0:
                range_low = min(finals)
                range_high = max(finals)
            else:
                range_low = "None."
                range_high = "Try a higher minimum ppm."
            self.annotation = self._spectra_ax.annotate(
                "X = {0:.4f}\nY = {1:.4f}\nID = {2}\nDrift Range = {3}-{4}".format
                (mz, intensity, lipid_id, range_low, range_high), xy=(mz, intensity), xycoords='data',
                va='bottom', ha='left',
                bbox=dict(boxstyle='square, pad=0.3', facecolor='white'))
        else:
            self.annotation = self._spectra_ax.annotate(
                "X = {0:.4f}\nY = {1:.4f}\nID = {2}".format
                (mz, intensity, lipid_id), xy=(mz, intensity), xycoords='data',
                va='bottom', ha='left',
                bbox=dict(boxstyle='square, pad=0.3', facecolor='white'))
        self.spectra_canvas.draw()

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
        if pickeditem is None or pickeditem == '':
            print('Please choose a map')
            return 0
        theMap = self.Maps[pickeditem.text()]
        df = pd.DataFrame(theMap)
        name = self.mmcWindow.exportMapData_filename.text()
        dirpath = QFileDialog.getExistingDirectory(self, 'Select a directory to export')
        if dirpath != '':
            pathfile = os.path.join(dirpath, name + '.csv')
            df.to_csv(pathfile, index=False)
        else:
            print('Please choose a directory')

    def MultiMapCompare_importMapData_Callback(self):
        path = QFileDialog.getOpenFileName(self, 'Select a map to import', filter='*.csv')
        if path[0] == '':
            print("No map to import.")
            return
        df = pd.read_csv(path[0])
        self.Maps[self.mmcWindow.importMapData_filename.text()] = df.to_numpy()
        self.MultiMapCompare_Display_Callback()
        return 0

    # executes when a load button is clicked
    def MultiMapCompare_LoadMap_Callback(self, button):
        if self.mmcWindow.pickitem:
            text = self.mmcWindow.pickitem.text()
            if text == '':
                print("Please select a map.")
                return 0
            if button is self.mmcWindow.slot1_load:
                self.MultiMapCompare_LoadMap_func(text, 1)
            elif button is self.mmcWindow.slot2_load:
                self.MultiMapCompare_LoadMap_func(text, 2)
            elif button is self.mmcWindow.slot3_load:
                self.MultiMapCompare_LoadMap_func(text, 3)
            elif button is self.mmcWindow.slot4_load:
                self.MultiMapCompare_LoadMap_func(text, 4)
            elif button is self.mmcWindow.slot5_load:
                self.MultiMapCompare_LoadMap_func(text, 5)
            elif button is self.mmcWindow.slot6_load:
                self.MultiMapCompare_LoadMap_func(text, 6)
        else:
            print("Choose an item from the listbox")

    def MultiMapCompare_LoadMap_func(self, item, num):
        if item == '':
            print("No item selected. Please double click on an item in the map listbox.")
            return
        self.mmcWindow.map_packet[num][6].setText(item)
        chosenMap = self.Maps[item]

        xend = len(chosenMap[0]) * (self.pixelSizeX / 1000)
        yend = len(chosenMap) * (self.pixelSizeY / 1000)

        if self.mmcWindow.map_packet[num][5]:
            self.mmcWindow.map_packet[num][5].removeWidget(self.mmcWindow.map_packet[num][3])
            self.mmcWindow.map_packet[num][5].removeWidget(self.mmcWindow.map_packet[num][0])

        self.mmcWindow.map_packet[num][0] = FigureCanvas(plt.figure(tight_layout=True))
        self.mmcWindow.map_packet[num][3] = NavigationToolbar(self.mmcWindow.map_packet[num][0], self)
        self.mmcWindow.map_packet[num][5].addWidget(self.mmcWindow.map_packet[num][3])
        self.mmcWindow.map_packet[num][5].addWidget(self.mmcWindow.map_packet[num][0])
        self.mmcWindow.map_packet[num][1] = self.mmcWindow.map_packet[num][0].figure.subplots()

        self.mmcWindow.map_packet[num][2] = self.mmcWindow.map_packet[num][1].imshow(chosenMap, cmap='jet',
                                                                                     interpolation='gaussian',
                                                                                     aspect=(yend / xend),
                                                                                     extent=[0, xend, 0, yend])
        self.mmcWindow.map_packet[num][1].set_xlabel('x, mm')
        self.mmcWindow.map_packet[num][1].set_ylabel('y, mm')
        self.mmcWindow.map_packet[num][4] = self.mmcWindow.map_packet[num][0].figure.colorbar(
            self.mmcWindow.map_packet[num][2])
        return 0

    def setIM(self):
        global isIM
        isIM = True
        self.all_drift_times.setCheckable(True)
        self.all_drift_times.setChecked(True)
        self.one_drift_time.setCheckable(True)
        self.drift_time.setDisabled(False)
        self.drift_scrollbar.setDisabled(False)

    def setMS(self):
        global isIM
        isIM = False
        self.all_drift_times.setChecked(False)
        self.all_drift_times.setCheckable(False)
        self.one_drift_time.setCheckable(False)
        self.drift_time.setDisabled(True)
        self.drift_scrollbar.setDisabled(True)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = MainGUIobject()
    widget = QtWidgets.QStackedWidget()
    widget.addWidget(main_window)
    widget.setFixedWidth(MW_width)
    widget.setFixedHeight(MW_height)
    widget.setWindowTitle("image_inspector")
    widget.show()
    sys.exit(app.exec_())
