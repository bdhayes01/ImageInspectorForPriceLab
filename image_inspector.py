# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 15:58:09 2021

@author: wtd14
"""

# TODO: for 10th of July week:
# 1. Group all functions based on where in the ui layout they are.
# 2. Make an outline of the different documents that you will write.
#       a. User manual (Quickstart guide and more complicated functions.
#       b. Input formatting notes and examples.
#       c. In-depth code comments (to be written in the python editor, using best practices)
# 3. Write Image Inspector documents
# 4. Make it crash proof, so that the program will not terminate early under any circumstances.
# 5. Make ML heuristic scorer.
#           Questions:
# 6. (JC)Should the noise button be implemented?
# 7. (JC)The standard dev increases when the image is flipped or rotated. Is that okay?
# 8. (JC)What color should the spectra plot colormap be?
#       a. Make examples for this in colab that can be shown to JC. Screenshot them!
# 9. Put the documents onto the lab website online.
# 10. (JC)Ask if the mass up/ down should add 1 or the m/z spacing?
# 11. Put a form for reporting errors on the lab website.

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
        self.pick_item = None

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

        # color_bars
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
        self.pick_item = item

    def display_mmcWindow(self):
        self.show()


loaded_ui_main = uic.loadUiType(mainWindow_ui_path)[0]


class MainGUIobject(QtWidgets.QMainWindow, loaded_ui_main):
    def __init__(self, parent=None):  # Initialization of the code
        QtWidgets.QMainWindow.__init__(self, parent)
        super(MainGUIobject, self).__init__()
        self.con_img = None
        self.toolbar = None
        self.axes = None
        self.pixel_size_x = None
        self.pixel_size_y = None
        self.binI = None
        self.ROI_outline = {}
        self.ROI_Map_Data = {}
        self.picked_point_data = None
        self.map_data = None
        self.label = None
        self.scale_factor = None
        self.spectra_toolbar = None
        self.annotation = None
        self.mz_vals = None
        self.intensity = None
        self.drifts = None
        self.chosen_data_iso = None
        self.view = None
        self.iso_view = None
        self.original_image = None
        self.ROIData = None
        self.setupUi(self)
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowSystemMenuHint)

        # plots
        self._spectra_ax = None
        self._con_ax = False

        # load in default ID file
        self.ids_pd = pd.read_csv(DESI_ID_path)

        # multi map compare instance
        self.mmcWindow = MultiMapCompareobject()

        # Function connections
        self.find_file.clicked.connect(self.find_data_file)
        self.find_file_mzOI.clicked.connect(self.find_file_mzOI_Callback)
        self.start_cube.clicked.connect(self.process_file)
        self.pick_point.clicked.connect(self.pick_point_Callback)
        self.mass_up.clicked.connect(self.mass_up_callback)
        self.mass_down.clicked.connect(self.mass_down_callback)
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
        self.ROI_listbox.itemDoubleClicked.connect(self.roi_listbox_callback)
        self.importROI.clicked.connect(self.import_roi_callback)
        self.deleteROIbutton.clicked.connect(self.delete_roi_callback)
        self.clearROIbutton.clicked.connect(self.clear_roi_callback)
        self.exportROI_val.clicked.connect(self.export_roi_spectra_callback)
        self.find_IDlist.clicked.connect(self.find_IDlist_Callback)
        self.Map_listbox.itemDoubleClicked.connect(self.map_listbox_callback)
        self.exportConcMap.clicked.connect(self.export_ConcMap_Callback)
        self.exportIsotopeMap.clicked.connect(self.export_IsotopeMap_Callback)
        self.deleteMapbutton.clicked.connect(self.delete_map_callback)
        self.clearMapListboxbutton.clicked.connect(self.clear_map_callback)
        self.extract_Map_mzOI.clicked.connect(self.mzOI_extractMap_Callback)
        self.drift_scrollbar.sliderMoved.connect(self.drift_scrollbar_callback)
        self.drift_time.valueChanged.connect(self.drift_time_callback)
        self.one_drift_time.toggled.connect(self.button_changed_callback)
        self.all_drift_times.toggled.connect(self.button_changed_callback)
        self.reset_scatter.clicked.connect(self.reset_scatter_callback)
        self.reset_image.clicked.connect(self.reset_orig_image)
        self.set_mz_minmax.clicked.connect(self.change_mz)
        self.flipButton.clicked.connect(self.flip_figure)
        self.rotate_right_button.clicked.connect(self.right_rotate)
        self.rotate_left_button.clicked.connect(self.left_rotate)

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
        self.millimeter.setChecked(True)
        self.massplusone.setChecked(True)
        self.all_drift_times.setChecked(True)
        self.ideal_ratio.setText('1')
        QListWidgetItem('ROI list appears here', self.ROI_listbox)
        QListWidgetItem('Map list appears here', self.Map_listbox)

        self.file_name = ""
        self.ROI = {}
        self.ROIplots = {}
        self.ROI_Mass = {}
        self.ROIcount = 0
        self.Maps = {}
        self.map_count = 0
        self.Map_listselect_text = ""

        # canvas
        self.spectra_canvas = None

        # ROI
        self.outline = None
        self.ROI_listselect_text = ""

        # multipmap compare variables
        self.conc_map_data = None
        self.isotope_map_data = None

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

    #####  Functions in the top boxes  #####
    def find_data_file(self):
        self.file_name = QFileDialog.getOpenFileName(self, 'Pick Data Cube', filter='*.mat, *.h5, *.bin *.csv')
        self.working_file_name.setText(self.file_name[0])

    # setIM sets the global isIM variable to true and enables all buttons needed to inspect IM data. 
    def setIM(self):
        global isIM
        isIM = True
        self.all_drift_times.setCheckable(True)
        self.all_drift_times.setChecked(True)
        self.one_drift_time.setCheckable(True)
        self.drift_time.setDisabled(False)
        self.drift_scrollbar.setDisabled(False)

    # setMS sets the global isIM variable to false and disables any buttons not needed to inspect MS data. 
    def setMS(self):
        global isIM
        isIM = False
        self.all_drift_times.setChecked(False)
        self.all_drift_times.setCheckable(False)
        self.one_drift_time.setCheckable(False)
        self.drift_time.setDisabled(True)
        self.drift_scrollbar.setDisabled(True)

    def process_file(self):
        file_name = self.file_name[0]
        self.reset_all()
        print("Working to read datacube")
        if file_name.endswith('.mat'):
            print("This code can't process .mat files")
            return
        elif file_name.endswith('.h5'):
            print("This code can't process .h5 files")
            return
        elif file_name.endswith('.csv'):
            data = ((pd.read_csv(file_name, header=None)).to_numpy(numpy.float32)).flatten()
            # Opens the input csv file and converts the file to the correct format.
        elif file_name.endswith('.bin'):
            file = open(file_name)
            data = np.fromfile(file, dtype=np.float32)
            file.close()
        else:
            print('Unexpected file extension')
            return
        try:
            if isIM:
                self.process_im_data(data)
            elif self.MSDataButton.isChecked():
                self.process_ms_data(data)
            else:
                print("Please select whether the file is IM or MS Data")
        except IndexError:
            # If an error is thrown it could be because the user selected the incorrect input format
            if isIM:
                print("Did you mean to select MS Data?")
            else:
                print("Did you mean to select IM Data?")

    # reset_all resets
    def reset_all(self):
        if self.micrometer.isChecked():
            self.scale_factor = 1e3
            self.label = 'μm'
        elif self.millimeter.isChecked():
            self.scale_factor = 1
            self.label = 'mm'
        elif self.centimeter.isChecked():
            self.scale_factor = 0.1
            self.label = 'cm'

        while self.iso_plot.count():
            child = self.iso_plot.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        self.iso_view = None

        self.chosen_data_iso = None

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
        self.ROIcount = 0
        self.ROIcountbox.setText("0")
        self.ROI_listbox.clear()
        QListWidgetItem('ROI list appears here', self.ROI_listbox)

        self.clear_map_callback()

        self.picked_point_data = None
        self.pick_IDthreshold.setValue(20)
        self.pick_IDthreshold.setMaximum(1000)

        self.zmax_isotope.setMaximum(99)
        self.zmax_isotope.setMinimum(0)
        self.zmin_isotope.setMaximum(99)
        self.zmin_isotope.setMinimum(0)
        self.mz_vals = None
        self.drifts = None
        self.intensity = None

    def process_ms_data(self, data):
        num_lines = data[0]
        num_scans = data[1]
        self.pixel_size_x = data[2]
        self.pixel_size_y = data[3]

        line_num = 0
        i = 4

        mz_vals = []
        intensities = []
        image_data = []
        map_data = []

        while line_num < num_lines:
            scan_num = 0
            line = []
            map_line = []
            while scan_num < num_scans:
                scan = 0
                num_data_points = data[i]
                i += 1
                curr_data_point = 0
                scan_val = []
                while curr_data_point < num_data_points:
                    curr_data_point += 1
                    mz_vals.append(data[i])
                    intensities.append(data[i + 1])
                    scan += data[i + 1]
                    scan_val.append([data[i], data[i + 1]])
                    i += 2
                scan_num += 1
                line.append(scan)
                map_line.append(scan_val)
            line_num += 1
            image_data.append(line)
            map_data.append(map_line)
        self.display_image(image_data, self.pixel_size_x, self.pixel_size_y)
        self.display_spectra(mz_vals, intensities, None)
        self.set_min_max_mz(mz_vals)
        self.map_data = map_data
        self.mz_vals = mz_vals
        self.intensity = intensities
        self.original_image = image_data

    def process_im_data(self, data):
        frame_num = data[0]
        file_num = data[1]
        self.pixel_size_x = data[2]
        self.pixel_size_y = data[3]

        num_frames = 0
        num_files = 0
        i = 4

        mz_vals = []
        intensity = []
        drifts = []
        chosen_data = []
        map_data = []

        while num_files < file_num:
            line_data = []
            map_line = []
            while num_frames < frame_num:
                total_drift_bins = data[i]
                curr_drift_bin = 0
                i += 1
                frame_val = []
                curr_val = 0
                while curr_drift_bin < total_drift_bins:
                    num_values = data[i]
                    i += 1
                    drift_time = data[i]
                    i += 1
                    for a in range(int(num_values)):
                        curr_val += data[i + 1]
                        intensity.append(data[i + 1])
                        drifts.append(drift_time)
                        mz_vals.append(data[i])
                        frame_val.append([data[i], data[i + 1], drift_time])
                        i += 2
                    curr_drift_bin += 1
                map_line.append(frame_val)
                line_data.append(curr_val)
                num_frames += 1
            chosen_data.append(line_data)
            num_files += 1
            num_frames = 0
            map_data.append(map_line)

        self.display_spectra(mz_vals, intensity, drifts)
        self.set_min_max_mz(mz_vals)
        self.display_image(chosen_data, self.pixel_size_x, self.pixel_size_y)

        self.mz_vals = mz_vals
        self.intensity = intensity
        self.drifts = drifts
        self.map_data = map_data
        self.original_image = chosen_data

    ##### Functions in the MS Image Controls box #####
    def mass_up_callback(self):
        if not self.view:
            return 0
        self.start.setText(
            str(float(self.start.text()) + 1.0))  # TODO: Ask if the mass up/ down should add 1 or the m/z spacing?
        if isIM:
            self.im_point()
        else:
            self.ms_point()
        return 0

    def mass_down_callback(self):
        if not self.view:
            return 0
        self.start.setText(str(float(self.start.text()) - 1.0))
        if isIM:
            self.im_point()
        else:
            self.ms_point()
        return 0

    def reset_orig_image(self):
        if self.original_image:
            self.picked_point_data = None
            self.display_image(self.original_image, self.pixel_size_x, self.pixel_size_y)
            self.start.clear()
            self.massbox.clear()
            self.isotope_map_data = None
            while self.iso_plot.count():
                child = self.iso_plot.takeAt(0)
                if child.widget():
                    child.widget().deleteLater()

    def export_ConcMap_Callback(self):
        if self.conc_map_data is None:
            print("There is no image to export")
            return
        # when clicking on the exportConcMap button, it will save the filename and concentration of the map
        self.Maps[self.exportConcMapname.text()] = self.conc_map_data
        self.refresh_map_listbox()
        self.map_count += 1

    def export_IsotopeMap_Callback(self):
        if self.isotope_map_data is None:
            print("There is no image to export")
            return
        # when clicking on the exportIsotopeMap button, it will save the filename and concentration of the map
        self.Maps[self.exportIsotopeMapname.text()] = self.isotope_map_data
        self.map_count += 1
        self.refresh_map_listbox()

    ##### Flipping functions #####
    def flip_figure(self):
        try:
            if self.original_image is not None:
                self.original_image.reverse()
            if self.map_data is not None:
                self.map_data.reverse()
            if self.conc_map_data is not None:
                self.conc_map_data.reverse()
                self.display_image(self.conc_map_data, self.pixel_size_x, self.pixel_size_y)
            if self.isotope_map_data is not None:
                self.isotope_map_data.reverse()
                self.display_iso_image(self.conc_map_data, self.isotope_map_data, self.pixel_size_x, self.pixel_size_y)
        except AttributeError:
            print("Error: Flipping an array with no values.\n"
                  "Please return to an image with data if you wish to flip the image.")
            return

    def right_rotate(self):
        self.rotate(True)

    def left_rotate(self):
        self.rotate(False)

    def rotate(self, isRight):
        if self.map_data is None:
            return
        self.pixel_size_x, self.pixel_size_y = self.pixel_size_y, self.pixel_size_x  # switch variables

        if isRight:
            self.map_data = self.rotate_right(self.map_data)
        else:
            self.map_data = self.rotate_left(self.map_data)

        if self.original_image:
            if isRight:
                self.original_image = self.rotate_right(self.original_image)
            else:
                self.original_image = self.rotate_left(self.original_image)

        if self.conc_map_data:
            if isRight:
                self.conc_map_data = self.rotate_right(self.conc_map_data)
            else:
                self.conc_map_data = self.rotate_left(self.conc_map_data)
            self.display_image(self.conc_map_data, self.pixel_size_x, self.pixel_size_y)

        if self.isotope_map_data:
            if isRight:
                self.isotope_map_data = self.rotate_right(self.isotope_map_data)
            else:
                self.isotope_map_data = self.rotate_left(self.isotope_map_data)
            self.display_iso_image(self.conc_map_data, self.isotope_map_data, self.pixel_size_x, self.pixel_size_y)
        return 0

    def rotate_right(self, orig_map):
        new_map = []
        for line in range(len(orig_map[0])):
            new_line = []
            for scan in range(len(orig_map) - 1, -1, -1):
                new_line.append(orig_map[scan][line])
            new_map.append(new_line)
        return new_map

    def rotate_left(self, orig_map):
        new_map = []
        for line in range(len(orig_map[0]) - 1, -1, -1):
            new_line = []
            for scan in range(len(orig_map)):
                new_line.append(orig_map[scan][line])
            new_map.append(new_line)
        return new_map

    ##### Original map scaling  #####
    def zmax_Callback(self):
        self.temp_max.setText(str(self.zmax.sliderPosition()))
        self.scale_image()

    def zmin_Callback(self):
        self.temp_min.setText(str(self.zmin.sliderPosition()))
        self.scale_image()

    def scale_image(self):
        if not self.picked_point_data:
            return 0
        image = self.picked_point_data
        data = []
        maximum = self.zmax.sliderPosition()
        minimum = self.zmin.sliderPosition()
        if minimum >= maximum:
            new_map = np.zeros((len(image), len(image[0])))
            self.display_image(new_map, self.pixel_size_x, self.pixel_size_y)
            return 0
        for line in image:
            newLine = []
            for scan in line:
                newFrame = 0
                if scan > maximum:
                    newFrame = maximum
                elif maximum >= scan >= minimum:
                    newFrame = scan
                newLine.append(newFrame)
            data.append(newLine)
        self.display_image(data, self.pixel_size_x, self.pixel_size_y)

    def temp_max_Callback(self):
        new_max = self.plot_con_temp_pressed_callback(self.temp_max.text())
        self.zmax.setValue(math.ceil(new_max))
        self.scale_image()

    def temp_min_Callback(self):
        new_min = self.plot_con_temp_pressed_callback(self.temp_min.text())
        self.zmin.setValue(math.ceil(new_min))
        self.scale_image()

    def plot_con_temp_pressed_callback(self, temp):
        try:
            new_val = float(temp)
        except ValueError:
            print("Please enter a numeric value.")
            new_val = float(self.min_int.text())
            return new_val
        if new_val > float(self.max_int.text()):
            print("Error: Trying to set the value above the max")
            new_val = float(self.max_int.text())
        elif new_val < float(self.min_int.text()):
            print("Error: Trying to set the value below the min")
            new_val = float(self.min_int.text())
        return new_val

    ##### Isotope map scaling  #####
    def zmax_isotope_Callback(self):
        self.max_iso.setText(str(self.zmax_isotope.sliderPosition()))
        self.scale_iso_image()

    def zmin_isotope_Callback(self):
        self.min_iso.setText(str(self.zmin_isotope.sliderPosition()))
        self.scale_iso_image()

    def scale_iso_image(self):
        if self.chosen_data_iso is None:  # This will be triggered only at the beginning
            return 0
        image = self.chosen_data_iso
        data = []
        maximum = self.zmax_isotope.sliderPosition()
        minimum = self.zmin_isotope.sliderPosition()
        if minimum >= maximum:
            new_map = np.zeros((len(image), len(image[0])))
            self.display_iso_image(new_map, new_map, self.pixel_size_x, self.pixel_size_y)
            return 0
        for line in image:
            new_line = []
            for frame in line:
                new_frame = 0
                if frame > maximum:
                    new_frame = maximum
                elif maximum >= frame >= minimum:
                    new_frame = frame
                new_line.append(new_frame)
            data.append(new_line)
        self.display_iso_image(self.picked_point_data, data, self.pixel_size_x, self.pixel_size_y)

    def temp_max_iso_Callback(self):
        new_max = self.iso_plot_temp_pressed_callback(self.max_iso.text())
        self.zmax_isotope.setValue(math.ceil(new_max))
        self.scale_iso_image()

    def temp_min_iso_Callback(self):
        new_min = self.iso_plot_temp_pressed_callback(self.min_iso.text())
        self.zmin_isotope.setValue(math.ceil(new_min))
        self.scale_iso_image()

    def iso_plot_temp_pressed_callback(self, temp):
        try:
            new_val = float(temp)
        except ValueError:
            print("Please enter a numeric value.")
            new_val = float(self.min_int_iso.text())
            return new_val
        if new_val > float(self.max_int_iso.text()):
            print("Error: Trying to set the value above the max")
            new_val = float(self.max_int_iso.text())
        elif new_val < float(self.min_int_iso.text()):
            print("Error: Trying to set the value below the min")
            new_val = float(self.min_int_iso.text())
        return new_val

    ##### Map Controls #####
    def pick_point_Callback(self):
        if self.micrometer.isChecked():
            self.scale_factor = 1e3
            self.label = 'μm'
        elif self.millimeter.isChecked():
            self.scale_factor = 1
            self.label = 'mm'
        elif self.centimeter.isChecked():
            self.scale_factor = 0.1
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
        map_data = self.map_data

        chosen_data = []
        chosen_data_plus_one = []
        chosen_data_plus_two = []
        max_intensity = 0
        max_intensity_plus_one = 0
        max_intensity_plus_two = 0

        for line in map_data:
            line_data = []
            line_data_plus_one = []
            line_data_plus_two = []
            for frame in line:
                the_val = 0
                the_val_plus_one = 0
                the_val_plus_two = 0
                for scan in frame:
                    if picked_point - max_diff <= scan[0] < picked_point + max_diff:
                        the_val += scan[1]
                        if the_val > max_intensity:
                            max_intensity = the_val
                    elif picked_point + (1 / ideal_ratio) - max_diff <= scan[0] < picked_point + (
                            1 / ideal_ratio) + max_diff:
                        the_val_plus_one += scan[1]
                        if the_val_plus_one > max_intensity_plus_one:
                            max_intensity_plus_one = the_val_plus_one
                    elif picked_point + (2 / ideal_ratio) - max_diff <= scan[0] < picked_point + (
                            2 / ideal_ratio) + max_diff:
                        the_val_plus_two += scan[1]
                        if the_val_plus_two > max_intensity_plus_two:
                            max_intensity_plus_two = the_val_plus_two
                line_data.append(the_val)
                line_data_plus_one.append(the_val_plus_one)
                line_data_plus_two.append(the_val_plus_two)
            chosen_data.append(line_data)
            chosen_data_plus_one.append(line_data_plus_one)
            chosen_data_plus_two.append(line_data_plus_two)

        m_zero_sum = 0
        m_one_sum = 0
        m_two_sum = 0
        if max_intensity != 0:
            denom = max_intensity + max_intensity_plus_one + max_intensity_plus_two
            m_zero_sum = round(max_intensity / denom, 4)
            m_one_sum = round(max_intensity_plus_one / denom, 4)
            m_two_sum = round(max_intensity_plus_two / denom, 4)

        self.Msumratio.setText(str(m_zero_sum))
        self.Mplusonesumratio.setText(str(m_one_sum))
        self.Mplustwosumratio.setText(str(m_two_sum))

        self.picked_point_data = None
        self.display_image(chosen_data, self.pixel_size_x, self.pixel_size_y)
        self.picked_point_data = chosen_data

        if self.massplusone.isChecked():
            self.chosen_data_iso = None
            self.display_iso_image(chosen_data, chosen_data_plus_one, self.pixel_size_x, self.pixel_size_y)
            self.chosen_data_iso = chosen_data_plus_one
        elif self.massplustwo.isChecked():
            self.chosen_data_iso = None
            self.display_iso_image(chosen_data, chosen_data_plus_two, self.pixel_size_x, self.pixel_size_y)
            self.chosen_data_iso = chosen_data_plus_two

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

        map_data = self.map_data

        image_data = []
        image_plus_one = []
        image_plus_two = []

        for line in map_data:
            line_vals = []
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
                line_vals.append(intensity)
                line_plus_one.append(intensity_plus_one)
                line_plus_two.append(intensity_plus_two)
            image_data.append(line_vals)
            image_plus_one.append(line_plus_one)
            image_plus_two.append(line_plus_two)

        self.picked_point_data = None
        self.display_image(image_data, self.pixel_size_x, self.pixel_size_y)
        self.picked_point_data = image_data

        not_used, m_zero_max_intensity = self.find_min_max_image(image_data)
        not_used, m_one_max_intensity = self.find_min_max_image(image_plus_one)
        not_used, m_two_max_intensity = self.find_min_max_image(image_plus_two)
        if m_zero_max_intensity != 0:
            denom = m_zero_max_intensity + m_one_max_intensity + m_two_max_intensity
            self.Msumratio.setText(str(round(m_zero_max_intensity / denom, 4)))
            self.Mplusonesumratio.setText(str(round(m_one_max_intensity / denom, 4)))
            self.Mplustwosumratio.setText(str(round(m_two_max_intensity / denom, 4)))

        if self.massplusone.isChecked():
            self.chosen_data_iso = None
            self.display_iso_image(image_data, image_plus_one, self.pixel_size_x, self.pixel_size_y)
            self.chosen_data_iso = image_plus_one
        if self.massplustwo.isChecked():
            self.chosen_data_iso = None
            self.display_iso_image(image_data, image_plus_two, self.pixel_size_x, self.pixel_size_y)
            self.chosen_data_iso = image_plus_two
        return 0

    def ROI_select_Callback_mask(self):
        if self.start.text() == '':
            print("You must choose or input a point to the 'Selected Mass' box first.")
            return
        if self.massbox.text() == '':
            print("There is no plot to select")
            return
        else:
            if self.outline:
                self.outline.disconnect()
            self.outline = roi.new_ROI(self.con_img)

    def exportROI_Callback(self):
        if self.outline is None:
            print("Please choose an ROI before processing")
            return
        if self.view:
            self.binI = self.outline.get_mask().astype(int)
            self.binI = np.flipud(self.binI)
            f = np.argwhere(np.ravel(self.binI, order='C'))[:, 0]
            self.ROI_outline[self.exportROIfilename.text()] = f
            self.numberpoints.setText(str(len(f)))
        self.ROIcount = self.ROIcount + 1
        self.ROIcountbox.setText(str(self.ROIcount))
        self.ROI[self.exportROIfilename.text()] = self.binI
        self.ROI_Mass[self.exportROIfilename.text()] = float(self.massbox.text())
        self.ROI_Map_Data[self.exportROIfilename.text()] = self.map_data
        self.refresh_ROI_listbox()
        return 0

    ##### Drift time box functions #####
    def drift_scrollbar_callback(self):
        val = self.drift_scrollbar.sliderPosition()
        self.drift_time.setValue(val)
        if self.drifts is None:
            return
        if self.one_drift_time.isChecked():
            self.show_mz_map(val)

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

    def show_mz_map(self, val):
        drift_vals = np.where(val == np.asarray(self.drifts), self.drifts, 0)
        drift_vals = drift_vals.nonzero()

        mz_vals = np.asarray(self.mz_vals)[drift_vals]
        intensity = np.asarray(self.intensity)[drift_vals]

        the_max = self.max_mz.value()
        the_min = self.min_mz.value()

        index = np.where(the_max >= mz_vals, mz_vals, 0)
        index = index.nonzero()

        mz_vals = np.asarray(mz_vals)[index]
        intensity = np.asarray(intensity)[index]

        index = np.where(mz_vals >= the_min, mz_vals, 0)
        index = index.nonzero()

        mz_vals = np.asarray(mz_vals)[index]
        intensity = np.asarray(intensity)[index]

        self.display_spectra(mz_vals, intensity, None)
        return 0

    def change_mz(self):
        max_mz = self.max_mz.value()
        min_mz = self.min_mz.value()
        mz_vals = self.mz_vals
        intensities = self.intensity

        if mz_vals is None:
            return 0

        processed_mz = []
        processed_intens = []

        if isIM:
            drifts = self.drifts
            processed_drift = []
            for i in range(len(mz_vals)):
                if max_mz >= mz_vals[i] >= min_mz:
                    processed_mz.append(mz_vals[i])
                    processed_intens.append(intensities[i])
                    processed_drift.append(drifts[i])
            self.display_spectra(processed_mz, processed_intens, processed_drift)
        else:
            for i in range(len(mz_vals)):
                if max_mz >= mz_vals[i] >= min_mz:
                    processed_mz.append(mz_vals[i])
                    processed_intens.append(intensities[i])
            self.display_spectra(processed_mz, processed_intens, None)
        try:
            self.set_min_max_mz(processed_mz)
        except ValueError:
            print("Error: There is no spectra to plot")
            return

    def reset_scatter_callback(self):
        if isIM:
            self.one_drift_time.setChecked(False)
            self.all_drift_times.setChecked(True)
        if self.mz_vals is None:
            print("Error: There is no original plot. Please select a .bin file and press 'GO' ")
            return
        self.display_spectra(self.mz_vals, self.intensity, self.drifts)
        self.set_min_max_mz(self.mz_vals)

    ##### ROI Listbox functions #####
    def roi_listbox_callback(self, item):
        # the outline of the ROI should appear
        try:
            self.ROI[item.text()]
        except KeyError:
            print("This is not a ROI. Please choose a ROI. ")
            return
        self.ROI_listselect_text = item.text()
        self._con_ax = True
        self.ROIplots[self.ROI_listselect_text] = item
        return 0

    def import_roi_callback(self):
        if self.ROI_listselect_text == "":
            print("No item selected")
        else:
            try:
                chosen_val = self.ROI_Mass[self.ROI_listselect_text]
            except ValueError:
                print("Error: The chosen m/z must be numeric")
                return
            ppm = self.ppm_calc(chosen_val)

            self.ROIData = []
            mz_vals = []
            intensities = []
            drifts = None

            outline = self.ROI_outline[self.ROI_listselect_text]
            map_data = self.ROI_Map_Data[self.ROI_listselect_text]
            if self.drifts is None:
                counter = 0
                for line in map_data:
                    for scan in line:
                        if counter not in outline:
                            counter += 1
                            continue
                        else:
                            for point in scan:
                                mz = point[0]
                                if chosen_val + ppm >= mz >= chosen_val - ppm:
                                    mz_vals.append(mz)
                                    intensities.append(point[1])
                            counter += 1
                for i in range(len(mz_vals)):
                    self.ROIData.append([mz_vals[i], intensities[i]])
            else:
                drifts = []
                counter = 0
                for line in map_data:
                    for frame in line:
                        if counter not in outline:
                            counter += 1
                            continue
                        else:
                            for point in frame:
                                mz = point[0]
                                if chosen_val + ppm >= mz >= chosen_val - ppm:
                                    mz_vals.append(mz)
                                    intensities.append(point[1])
                                    drifts.append(point[2])
                            counter += 1
                for i in range(len(mz_vals)):
                    self.ROIData.append([mz_vals[i], intensities[i], drifts[i]])

            self.display_spectra(mz_vals, intensities, drifts)
            self.set_min_max_mz(mz_vals)

    def export_roi_spectra_callback(self):
        if self.ROI_listselect_text == "":
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

    def delete_roi_callback(self):
        if self.ROIcount == 0:
            self.ROI_listselect_text = ""
            return
        else:
            if self._con_ax:
                if self.ROI_listselect_text in self.ROIplots:
                    del self.ROIplots[self.ROI_listselect_text]
                    del self.ROI[self.ROI_listselect_text]
                    self.ROIcount = self.ROIcount - 1
                    self.ROIcountbox.setText(str(self.ROIcount))
                    self.refresh_ROI_listbox()
                    print('plot removed')
                    if isIM:
                        self.im_point()
                    else:
                        self.ms_point()
                else:
                    print('Item does not exist. Please double click on another item')
        self.ROI_listselect_text = ""

    def clear_roi_callback(self):
        self.ROIplots.clear()
        self.ROI.clear()
        self.ROIcount = 0
        self.ROIcountbox.setText("0")
        self.refresh_ROI_listbox()
        self.reset_scatter_callback()
        self.ROI_listselect_text = ""
        if isIM:
            self.im_point()
        else:
            self.ms_point()

    def refresh_ROI_listbox(self):
        if len(self.ROI) == 0:
            self.ROI_listbox.clear()
            QListWidgetItem('ROI list appears here', self.ROI_listbox)
        else:
            self.ROI_listbox.clear()
            list_box_items = list(self.ROI.keys())
            for i in range(len(list_box_items)):
                self.ROI_listbox.addItem(list_box_items[i])

    ##### Map Listbox functions #####
    def delete_map_callback(self):
        if self.map_count == 0:
            return
        else:
            if self.Map_listselect_text:
                del self.Maps[self.Map_listselect_text]
                self.refresh_map_listbox()
                self.map_count -= 1
                print('Map removed')
            else:
                print('Item does not exist. Please double click on another item')

    def clear_map_callback(self):
        if self.map_count == 0:
            return
        else:
            self.Maps.clear()
            self.map_count = 0
            self.refresh_map_listbox()

    def map_listbox_callback(self, item):
        self.Map_listselect_text = item.text()

    def refresh_map_listbox(self):
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
            list_box_items = list(self.Maps.keys())
            for i in range(len(list_box_items)):
                self.Map_listbox.addItem(list_box_items[i])

    ##### m/z of interest box functions #####
    def find_file_mzOI_Callback(self):
        self.mzOI_listname.setText(QFileDialog.getOpenFileName(self, 'Pick list: m/z of interest', filter='*.csv')[0])

    def mzOI_extractMap_Callback(self):
        if self.mzOI_listname.text() == '':
            print("Please choose a .csv file with m/z values")
            return 0
        if self.mz_vals is None:
            print("Please choose a .bin file to process")
            return 0
        try:
            mzOI = pd.read_csv(self.mzOI_listname.text())
            mzOI = mzOI.to_numpy()
        except IOError:
            print("Please choose a .csv file with m/z values")
            return 0

        mz_vals = []
        intensity = []
        drifts = None
        if isIM:
            drifts = []
        for mz in mzOI:
            diff = self.ppm_calc(mz)
            for i in range(len(self.mz_vals)):
                if mz + diff >= self.mz_vals[i] >= mz - diff:
                    mz_vals.append(self.mz_vals[i])
                    intensity.append(self.intensity[i])
                    if isIM:
                        drifts.append(self.drifts[i])
        self.display_spectra(mz_vals, intensity, drifts)

    ##### m/z ID box functions #####
    def find_IDlist_Callback(self):
        file_name = QFileDialog.getOpenFileName(self, 'Pick ID List', filter='*.csv')
        self.IDlist_name.setText(file_name[0])
        if file_name[0] == '':
            print("Please select a file to import as the ID list.")
            return
        self.ids_pd = pd.read_csv(file_name[0])

    ##### Image display functions #####
    def display_iso_image(self, zero_image, imageData, pixelSizeX, pixelSizeY):
        # This needs to be a different function for both the original point and the scaling functions
        count = len(zero_image) * len(zero_image[0])
        self.numberpoints.setText(str(count))

        x_end = len(imageData[0]) * (pixelSizeX / 1000)
        y_end = len(imageData) * (pixelSizeY / 1000)

        while self.iso_plot.count():
            child = self.iso_plot.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

        zero_for_dev = np.asarray(zero_image).flatten()
        std_deviation = round(float(np.std(zero_for_dev)), 4)
        self.Msumstandard_error.setText(str(std_deviation))

        iso_data = self.isotope_scalar(zero_image, imageData)
        self.isotope_map_data = iso_data
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
        self.iso_plot.addWidget(self.iso_view)
        con_img2 = self.axes.imshow(iso_data, cmap='inferno', aspect=(y_end / x_end),
                                    extent=[0, x_end * self.scale_factor, 0, y_end * self.scale_factor])
        self.axes.set_xlabel("x, " + self.label)
        self.axes.set_ylabel("y, " + self.label)
        plt.colorbar(con_img2)
        self.iso_view.draw()
        if not self.chosen_data_iso:
            the_min, the_max = self.find_min_max_image(imageData)
            # If this needs to be the scaled image just change the above line
            self.max_int_iso.setText(str(the_max))
            self.min_int_iso.setText(str(the_min))
            self.max_iso.setText(str(the_max))
            self.min_iso.setText(str(the_min))
            self.zmax_isotope.setMaximum(math.ceil(the_max))
            self.zmax_isotope.setMinimum(math.floor(the_min))
            self.zmin_isotope.setMaximum(math.ceil(the_max))
            self.zmin_isotope.setMinimum(math.floor(the_min))
            self.zmax_isotope.setValue(math.ceil(the_max))
            self.zmin_isotope.setValue(math.floor(the_min))

    def isotope_scalar(self, m_zero_intensity, isotope_intensity):
        new_data = []
        for line in range(len(m_zero_intensity)):
            the_line = []
            orig_line = m_zero_intensity[line]
            iso_line = isotope_intensity[line]
            for scan in range(len(orig_line)):
                if orig_line[scan] == 0:
                    ratio = 0
                else:
                    ratio = iso_line[scan] / orig_line[scan] + iso_line[scan]
                the_line.append(ratio)
            new_data.append(the_line)
        return new_data

    def display_image(self, image_data, pixel_size_x, pixel_size_y):
        x_end = len(image_data[0]) * (pixel_size_x / 1000)
        y_end = len(image_data) * (pixel_size_y / 1000)

        if self.view:
            self.plot_con.removeWidget(self.view)
        self.view = FigureCanvas(Figure(figsize=(5, 3)))
        self.axes = self.view.figure.subplots()
        self.toolbar = NavigationToolbar(self.view, self)
        self.plot_con.addWidget(self.view)
        self.con_img = self.axes.imshow(image_data, cmap='jet', aspect=(y_end / x_end),
                                        extent=[0, x_end * self.scale_factor, 0, y_end * self.scale_factor])
        plt.colorbar(self.con_img)
        self.axes.set_title('Points In Selected Region')
        self.axes.set_xlabel("x, " + self.label)
        self.axes.set_ylabel("y, " + self.label)
        self.view.draw()
        self.conc_map_data = image_data
        if not self.picked_point_data:
            the_min, the_max = self.find_min_max_image(image_data)
            self.max_int.setText(str(the_max))
            self.min_int.setText(str(the_min))
            self.temp_max.setText(str(the_max))
            self.temp_min.setText(str(the_min))
            self.zmax.setMaximum(math.ceil(the_max))
            self.zmax.setMinimum(math.floor(the_min))
            self.zmin.setMaximum(math.ceil(the_max))
            self.zmin.setMinimum(math.floor(the_min))
            self.zmax.setValue(math.ceil(the_max))
            self.zmin.setValue(math.floor(the_min))
            self.picked_point_data = image_data

    ##### Spectra display functions #####
    def display_spectra(self, mz_vals, intensity, drifts=None, pt_size=.01):
        while self.plot_spectra.count():  # This solved the double figure problem
            child = self.plot_spectra.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        plt.close('all')

        if 100 / len(mz_vals) > pt_size:
            pt_size = 100 / len(mz_vals)

        self.spectra_canvas = FigureCanvas(plt.figure(tight_layout=True))
        self.spectra_canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.spectra_canvas.setFocus()
        self.spectra_toolbar = NavigationToolbar(self.spectra_canvas, self)
        self.plot_spectra.addWidget(self.spectra_toolbar)
        self.plot_spectra.addWidget(self.spectra_canvas)
        self._spectra_ax = self.spectra_canvas.figure.subplots()
        if drifts is not None:
            spectra = self._spectra_ax.scatter(mz_vals, intensity, s=pt_size, c=drifts, cmap="turbo", alpha=0.75, picker=True)
            plt.colorbar(spectra).set_label('Drift times')
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
            self._spectra_ax.scatter(mz_vals, intensity, color='#393424', s=pt_size, alpha=0.75, picker=True)
            # Can change to peaks here
        self._spectra_ax.set_title('Points In Selected Region')
        self._spectra_ax.set_xlabel('m/z')
        self._spectra_ax.set_ylabel('intensity')
        self.spectra_canvas.mpl_connect('pick_event', self.data_cursor_click)
        plt.ylabel('intensity')
        plt.xlabel('m/z')

    def data_cursor_click(self, event):
        theXmOverZ = event.mouseevent.lastevent.xdata
        theY = event.mouseevent.lastevent.ydata
        self.start.setText("%.5f" % theXmOverZ)
        self.spectra_annotation(theXmOverZ, theY)

    def spectra_annotation(self, mz, intensity):
        diff = self.ppm_calc(mz)
        lipid_map = {}
        lipid_id = "not defined"
        if self.ids_pd is not None:
            mz_vals = self.ids_pd['m/z']
            lipid_ids = self.ids_pd['Lipid ID']
            for i in range(mz_vals.size):
                mz_val = mz_vals[i]
                if (mz + diff) > mz_val > (mz - diff):
                    difference = abs(mz - mz_val)
                    lipid_map[difference] = lipid_ids[i]
        if len(lipid_map.keys()) > 0:
            key = min(lipid_map.keys())
            lipid_id = lipid_map[key]

        self.ID_Output_Box.setText(lipid_id)

        if self.annotation is not None:
            self.annotation.remove()

        if isIM:
            mz_vals = np.asarray(self.mz_vals)
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

    def set_min_max_mz(self, mz_vals):
        try:
            self.min_mz.setRange(min(mz_vals), max(mz_vals))
            self.max_mz.setRange(min(mz_vals), max(mz_vals))
            self.min_mz.setValue(min(mz_vals))
            self.max_mz.setValue(max(mz_vals))
        except TypeError:
            print("Error: There is no plot. Please select a .bin file and press 'GO' ")
            return

    ##### Misc. functions #####
    def find_min_max_image(self, image_data):
        the_min = sys.maxsize
        the_max = 0
        for line in image_data:
            line_min = min(line)
            line_max = max(line)
            if line_max > the_max:
                the_max = line_max
            if line_min < the_min:
                the_min = line_min
        return the_min, the_max

    def ppm_calc(self, mz_val):
        if float(self.pick_IDthreshold.value()) == 0:
            return 0
        return mz_val * float(self.pick_IDthreshold.value()) / 1e6

    ##########################################
    # Multi Map Compare functions (things related to the 2nd window)
    def MultiMapCompare_Display_Callback(self):
        self.mmcWindow.pick_item = None
        self.mmcWindow.Map_listbox_mmcWindow.clear()
        list_box_items = list(self.Maps.keys())
        for i in range(len(list_box_items)):
            self.mmcWindow.Map_listbox_mmcWindow.addItem(list_box_items[i])
        self.mmcWindow.display_mmcWindow()

    def MultiMapCompare_exportMapData_Callback(self):
        picked_item = self.mmcWindow.pick_item
        if picked_item is None or picked_item == '':
            print('Please choose a map')
            return 0
        df = pd.DataFrame(self.Maps[picked_item.text()])
        name = self.mmcWindow.exportMapData_filename.text()
        dir_path = QFileDialog.getExistingDirectory(self, 'Select a directory to export')
        if dir_path != '':
            pathfile = os.path.join(dir_path, name + '.csv')
            df.to_csv(pathfile, index=False)
        else:
            print('Please choose a directory')

    def MultiMapCompare_importMapData_Callback(self):
        path = QFileDialog.getOpenFileName(self, 'Select a map to import', filter='*.csv')
        if path[0] == '':
            print("No map to import.")
            return
        self.Maps[self.mmcWindow.importMapData_filename.text()] = pd.read_csv(path[0]).to_numpy()
        self.MultiMapCompare_Display_Callback()
        return 0

    # executes when a load button is clicked
    def MultiMapCompare_LoadMap_Callback(self, button):
        if self.mmcWindow.pick_item:
            text = self.mmcWindow.pick_item.text()
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
        chosen_map = self.Maps[item]

        x_end = len(chosen_map[0]) * (self.pixel_size_x / 1000)
        y_end = len(chosen_map) * (self.pixel_size_y / 1000)

        if self.mmcWindow.map_packet[num][5]:
            self.mmcWindow.map_packet[num][5].removeWidget(self.mmcWindow.map_packet[num][3])
            self.mmcWindow.map_packet[num][5].removeWidget(self.mmcWindow.map_packet[num][0])

        self.mmcWindow.map_packet[num][0] = FigureCanvas(plt.figure(tight_layout=True))
        self.mmcWindow.map_packet[num][3] = NavigationToolbar(self.mmcWindow.map_packet[num][0], self)
        self.mmcWindow.map_packet[num][5].addWidget(self.mmcWindow.map_packet[num][3])
        self.mmcWindow.map_packet[num][5].addWidget(self.mmcWindow.map_packet[num][0])
        self.mmcWindow.map_packet[num][1] = self.mmcWindow.map_packet[num][0].figure.subplots()

        self.mmcWindow.map_packet[num][2] = self.mmcWindow.map_packet[num][1].imshow(chosen_map, cmap='jet',
                                                                                     interpolation='gaussian',
                                                                                     aspect=(y_end / x_end),
                                                                                     extent=[0, x_end, 0, y_end])
        self.mmcWindow.map_packet[num][1].set_xlabel('x, mm')
        self.mmcWindow.map_packet[num][1].set_ylabel('y, mm')
        self.mmcWindow.map_packet[num][4] = self.mmcWindow.map_packet[num][0].figure.colorbar(
            self.mmcWindow.map_packet[num][2])
        return 0


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
