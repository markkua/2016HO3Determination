# python3
# -*- coding: utf-8 -*-

import math
import sys
import time
from typing import List

import numpy as np
from vtk import *
from jplephem.spk import SPK
from PyQt5.QtCore import *
from PyQt5.QtGui import *
#from PyQt5.QtWidgets import QMainWindow, QApplication, QSplitter, QFrame, QVBoxLayout, QStatusBar, QLabel, QCalendarWidget
from PyQt5.QtWidgets import *
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

# TODO 减小体积


Colors = vtk.vtkNamedColors()


class MainWindow(QMainWindow):
    displaying: bool = False  # 是否在循环刷新vtk控件
    frame_rate = 25  # 帧率
    delta_tdb = 0.2  # 帧之间的时间差  TODO 调整速率

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        # 设置标题，大小，图标
        self.setWindowTitle("2016HO3")
        self.resize(1200, 700)
        self.setWindowIcon(QIcon("data/icon.jpg"))
        # 初始化数据
        self.dataProvider = DataProvider("data/de430.bsp", "data/status.eph")
        self.dataProvider.current_tdb = date2TDB(QDate.currentDate())  # todo  当前日期
        # 初始化VTK控件
        self.vtkWidget = MyVTKWidget()  # VTK控件对象
        self._insertPlanets()  # 插入行星显示主体
        # 初始化 UI
        self._initUI()
        # connections
        self.okButton.clicked.connect(self._on_OKButton_clicked)
        self.calendar.selectionChanged.connect(self._onCalChanged)
        log("MainWindow Ready.", "info")
        self.window_status_label.setText("Ready")
        return

    def _initUI(self):
        # 分割窗口，并添加控件
        # 右半部分
        rSplitter = QSplitter(Qt.Vertical)
        rFrame = QFrame()
        rLayout = QVBoxLayout()
        rLayout.addWidget(self.vtkWidget)
        rLayout.setContentsMargins(0, 0, 0, 0)  # ?
        rFrame.setLayout(rLayout)
        rSplitter.addWidget(rFrame)
        # status bar
        self.window_status_label = QLabel()
        self.window_status_label.setMinimumWidth(200)
        self.window_status_label.setAlignment(Qt.AlignRight)
        self.statusBar = QStatusBar()
        self.statusBar.setMinimumWidth(300)
        self.statusBar.addPermanentWidget(self.window_status_label)
        rSplitter.addWidget(self.statusBar)
        # 拉伸比例
        rSplitter.setStretchFactor(0, 9)
        rSplitter.setStretchFactor(1, 1)

        # 左边
        self.lSplitter = QSplitter(Qt.Vertical)
        lFrame = QFrame()
        lLayout = QVBoxLayout()
        # 大标题2016HO3
        lLabel0 = QLabel()
        lLabel0.setText("2016HO3")
        lLabel0.setFont(QFont("Microsoft YaHei", 40, QFont.Bold))
        lLabel0.setAlignment(Qt.AlignCenter)
        lLayout.addWidget(lLabel0)  # 0
        lLayout.setStretch(0, 1.5)
        # Select start date
        lLabel = QLabel()
        lLabel.setText("Select start date: ")
        lLabel.setFont(QFont("Microsoft YaHei", 20, QFont.Bold))
        lLabel.setAlignment(Qt.AlignCenter)
        lLayout.addWidget(lLabel)  # 1
        lLayout.setStretch(1, 1)
        # 日历
        # self.calendar = QCalendarWidget()
        # self.calendar.setGridVisible(True)
        # self.calendar.setDateRange(QDate(1900, 1, 1), QDate(2195, 12, 24))
        # self.calendar.setMaximumHeight(300)
        # self.calendar.setMinimumWidth(30)
        # lLayout.addWidget(self.calendar)  # 2
        self.speed_layout = QHBoxLayout()
        self.speedButton_origion = QPushButton("▶")
        self.speedButton_faster = QPushButton(">>>")
        self.speedButton_slower = QPushButton("<<<")
        self.speed_layout.addWidget(self.speedButton_origion)
        self.speed_layout.addWidget(self.speedButton_faster)
        self.speed_layout.addWidget(self.speedButton_slower)
        lLayout.addWidget(self.speed_layout)
        lLayout.setStretch(2, 4)
        # 日历标签
        self.calendar_label = QLabel()
        self.calendar_label.setFont(QFont("Microsoft YaHei", 10, QFont.Bold))
        self.calendar_label.setAlignment(Qt.AlignCenter)
        self._onCalChanged()
        lLayout.addWidget(self.calendar_label)  # 3
        lLayout.setStretch(3, 1)

        self.okButton = QPushButton("OK")
        self.okButton.setMinimumHeight(50)
        self.okButton.setFont(QFont("Microsoft YaHei", 20, QFont.Bold))
        lLayout.addWidget(self.okButton)  # 4
        lLayout.setStretch(4, 1)

        tempLabel = QLabel()
        tempLabel.setMaximumHeight(50)
        lLayout.addWidget(tempLabel)  # 5
        lLayout.setStretch(5, 1)

        lFrame.setLayout(lLayout)
        self.lSplitter.addWidget(lFrame)

        # 总布局
        splitter_main = QSplitter(Qt.Horizontal)
        splitter_main.addWidget(self.lSplitter)  # 左边
        splitter_main.addWidget(rSplitter)  # 右边
        splitter_main.setStretchFactor(0, 1)
        splitter_main.setStretchFactor(1, 5)
        self.setCentralWidget(splitter_main)
        return

    def _insertPlanets(self):
        # 太阳
        self.vtkWidget.renderer.AddActor(
            self.dataProvider.buildSphereSource('Sun', np.array([0, 0, 0]), 0.2, "Yellow"))
        # 地球
        self.vtkWidget.renderer.AddActor(
            self.dataProvider.buildSphereSource('Earth', np.array([100, 0, 0]), 0.09, "SkyBlue"))
        # 2016HO3
        self.vtkWidget.renderer.AddActor(
            self.dataProvider.buildSphereSource('2016HO3', np.array([100, 1, 0]), 0.03, "Pink"))
        # Mercury
        self.vtkWidget.renderer.AddActor(
            self.dataProvider.buildSphereSource('Mercury', np.array([100, 1, 0]), 0.05, "Gold"))
        # Venus
        self.vtkWidget.renderer.AddActor(
            self.dataProvider.buildSphereSource('Venus', np.array([100, 1, 0]), 0.04, "Green"))
        # TODO 改半径，改为实际大小
        return None

    def _on_OKButton_clicked(self):
        log("OK clicked", "info")
        if not self.displaying:
            self.displaying = True
            log("go, displaying=" + self.displaying.__str__(), "debug")
            self._display_circle()
        else:  # displaying
            self.displaying = False
            log("pause, displaying=" + self.displaying.__str__(), "debug")
        QApplication.processEvents()
        return

    def _display_circle(self):
        log("start displaying", "info")
        self.window_status_label.setText("Displaying")
        while self.displaying:
            frame_start = time.time()  # 帧开始时间
            frame_time = 1 / self.frame_rate  # 帧时长
            # update data (takes time)
            self._update_data()
            # 控制帧率
            time_spent = time.time() - frame_start
            time_left = frame_time - time_spent
            if time_left > 0:
                time.sleep(time_left)
            self._refresh_vtkDisplay()
            QApplication.processEvents()
        self.window_status_label.setText("Pause")
        log("displaying circle quit", "info")
        return

    def _update_data(self):
        current_tdb = self.dataProvider.current_tdb
        new_tdb = current_tdb + self.delta_tdb
        self.dataProvider.update_to_tdb(new_tdb)
        # todo
        return

    def _refresh_vtkDisplay(self):
        self.window_status_label.setText("refresh")
        self.vtkWidget.GetRenderWindow().Render()  # 刷新
        # Show current tdb
        self.statusBar.showMessage("TDB=%.4f" % self.dataProvider.current_tdb)
        QApplication.processEvents()
        return

    def _onCalChanged(self):
        date = self.calendar.selectedDate()
        tdb = date2TDB(date)
        self.calendar_label.setText("%s  TDB = %.2f" % (str(date.toPyDate()), tdb))
        # todo


# VTK控件对象
class MyVTKWidget(QVTKRenderWindowInteractor):
    def __init__(self):
        QVTKRenderWindowInteractor.__init__(self)
        self.renderer = vtk.vtkRenderer()
        self.GetRenderWindow().AddRenderer(self.renderer)
        self.interactor = self.GetRenderWindow().GetInteractor()
        # self.interactor.SetInteractorStyle(None)  # 禁用交互
        # todo 相机
        camera = vtk.vtkCamera()
        camera.SetViewAngle(40)
        camera.SetFocalPoint(0, 0, 0)
        camera.SetPosition(4, 0, 5)
        camera.SetViewUp(0, 0, 1)
        # print(camera)
        self.renderer.SetActiveCamera(camera)


# 读取2016HO3的文本星历
class EphemerisReader:
    """ 读取星历文件 """
    def __init__(self, filename: str):
        self.filename = filename
        print("ephReader filename = %s" % self.filename)  # TODO test
        with open(self.filename, 'r') as file:
            self.content = file.readlines()
        self.min_tdb, self.max_tdb = self._get_tdb_range()

    def _get_tdb_range(self):
        index = 0
        while self.content[index] != "$$SOE\n":
            if not index < self.content.__len__():
                print("indexStatus search error: $$SOE out of range")
            index += 1
        index += 1
        earliest_TAB, Position, Velocity = self._readStatus(index)
        while self.content[index] != "$$EOE\n":
            index += 3
        latest_TDB, Position, Velocity = self._readStatus(index - 3)
        return earliest_TAB, latest_TDB

    def indexStatus(self, targetTDB: float):
        # 搜索该日期对应记录的首行行号（如果不相等则向前取最近的记录）
        index = 0
        # 跳过Head
        while self.content[index] != "$$SOE\n":
            if not index < self.content.__len__():
                print("indexStatus search error: $$SOE out of range")
            index += 1
        index += 1
        # 数据开始
        while self.content[index] != "$$EOE\n":
            line = self.content[index]
            TDB1 = float(line[0:17])
            while index + 3 < self.content.__len__():
                line = self.content[index + 3]
                TDB2 = float(line[0:17])
                if TDB1 <= targetTDB < TDB2:
                    break
                else:
                    index += 3
            return index

    def _readStatus(self, lineIndex: int):
        # lineIndex: 3行为一个记录，Index指向记录首行
        if not lineIndex < self.content.__len__():
            print("_readStatus error: index out of range")
        try:
            line = self.content[lineIndex]
            TDB = float(line[0:17])
            Position = np.zeros(3)
            Velocity = np.zeros(3)
            line = self.content[lineIndex+1]
            j = 9
            for i in range(3):
                temp = (line[j:j + 22]).lower()
                Position[i] = float(temp)
                j += 24
            for i in range(3):
                temp = (line[j:j + 22]).lower()
                Velocity[i] = float(temp)
                j += 24
            return TDB, Position, Velocity
        except ValueError:
            print("ERROR: _readStatus | lineIndex=" + lineIndex.__str__())
            return None, None, None

    def getPosition(self, tdb: float):
        if not self.min_tdb <= tdb < self.max_tdb:
            print(self.__class__(), " tdb out of range")
        index = self.indexStatus(tdb)
        tdb1, Position1, Velocity = self._readStatus(index)
        tdb2, Position2, Velocity = self._readStatus(index+3)
        Position1 = (Position2 - Position1) / (tdb2 - tdb1) * (tdb - tdb1) + Position1
        return Position1


# 读取de430.bsp星历文件
class BspReader:
    def __init__(self, filename='de430.bsp'):
        self._kernel = SPK.open(filename)
        self._PlanetDict = {"Mercury": 1, "Venus": 2, "Earth": 3, "Mars": 4, "Jupiter": 5, "Saturn": 6, "Uranus": 7,
                            "Neptune": 8, "Pluto": 9, "Sun": 10}

    def getPosition(self, key: str, tdb: float):
        target = self._PlanetDict[key]
        position = self._kernel[0, target].compute(tdb)
        return position


class DataProvider:
    def __init__(self, bsp_file="data/de430.bsp", eph_file="data/status.eph"):
        # data
        self._sphereSourceDic = {}  # 所有的行星的sphereSource
        self._bspReader = BspReader(bsp_file)
        self._ephReader = EphemerisReader(eph_file)
        self.scale = 1e-8  # 缩放尺度
        self.min_tdb, self.max_tdb = self._ephReader.min_tdb, self._ephReader.max_tdb
        self.current_tdb = 0  # todo

    def buildSphereSource(self, key: str, Center: np.array(3), radius: float, color: str) -> vtk.vtkActor:
        if key in self._sphereSourceDic.keys():
            print("sphereSource already exists, can not insert " + key)
            return vtk.vtkActor()
        sphereSource = vtk.vtkSphereSource()
        # Make the surface smooth.
        sphereSource.SetPhiResolution(100)
        sphereSource.SetThetaResolution(100)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphereSource.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(Colors.GetColor3d(color))
        sphereSource.SetCenter(Center[0], Center[1], Center[2])
        sphereSource.SetRadius(radius)
        # 加入数据字典
        self._sphereSourceDic[key] = sphereSource
        return actor

    def update_to_tdb(self, tdb: float):
        if not self.min_tdb <= tdb < self.max_tdb:
            print("dataProvider error: tdb out of range")
            return
        for key in self._sphereSourceDic.keys():
            if key == '2016HO3':
                Center = self._ephReader.getPosition(tdb)
            else:
                Center = self._bspReader.getPosition(key, tdb)
            Center *= self.scale
            self._sphereSourceDic[key].SetCenter(Center[0], Center[1], Center[2])
        self.current_tdb = tdb
        return


def date2TDB(date: QDate):
        # 2458551.000000000 = A.D. 2019-Mar-08
        defaultTDB = 2458551
        defaultDate = QDate(2019, 3, 8)
        tdb = defaultTDB + defaultDate.daysTo(date)
        return tdb

# 输出日志的函数
def log(log: str, level: str):
    level = level.lower()
    if level == "debug":
        print("[Debug] %s" % log)
    elif level == "info":
        print("[Info]  %s" % log)
    elif level == "warn":
        print("[Warn]  %s" % log)
    elif level == "error":
        print("[error] %s" % log)
    else:
        print("[error] Log level error")


if __name__ == "__main__":
    app = QApplication(sys.argv)

    win = MainWindow()

    win.show()

    win.vtkWidget.interactor.Initialize()

    sys.exit(app.exec_())
