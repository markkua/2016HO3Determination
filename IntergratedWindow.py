# python3
# -*- coding: utf-8 -*-

from math import asin, acos, pi
import sys
import time
from typing import List

import numpy as np
from vtk import *  # TODO 轻量化
from jplephem.spk import SPK

import pandas as pd
import spiceypy as spice
from numpy.linalg import norm
from scipy.integrate import solve_ivp

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import QMainWindow, QApplication, QSplitter, QFrame, QWidget,\
     QVBoxLayout, QHBoxLayout, QStatusBar, QLabel, QCalendarWidget, QPushButton
# from PyQt5.QtWidgets import *
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
# TODO 减小体积


Colors = vtk.vtkNamedColors()

# TODO 同步显示TDB和日期


class MainWindow(QMainWindow):
    displaying: bool = False  # 是否在循环刷新vtk控件
    frame_rate = 25  # 帧率
    origion_delta_tdb = 0.05
    delta_tdb = origion_delta_tdb

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        # 设置标题，大小，图标
        self.setWindowTitle("2016HO3")
        self.resize(1200, 700)
        self.setWindowIcon(QIcon("data/icon.jpg"))
        # 初始化数据
        self.dataProvider = DataProvider("data/de430.bsp", "data/status.eph")
        self.dataProvider.current_tdb = date2TDB(QDate.currentDate())  # 设置为当前日期
        # 初始化VTK控件
        self.vtkWidget = MyVTKWidget()  # VTK控件对象
        self._init_planets()  # 插入行星显示主体
        self._update_data()  # 更新球体坐标到当前日期
        # 初始化 UI
        self._initUI()
        # connections
        self.okButton.clicked.connect(self._on_OKButton_clicked)
        self.calendar.selectionChanged.connect(self._onCalChanged)
        self.speedButton_moreSlower.clicked.connect(self._on_speedButton_moreSlower)
        self.speedButton_slower.clicked.connect(self._on_speedButton_slower)
        self.speedButton_origion.clicked.connect(self._on_speedButton_origion)
        self.speedButton_faster.clicked.connect(self._on_speedButton_faster)
        self.speedButton_moreFaster.clicked.connect(self._on_speedButton_moreFaster)
        # log
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
        lLabel0.setFont(QFont("Microsoft YaHei", 30, QFont.Bold))
        lLabel0.setAlignment(Qt.AlignCenter)
        lLayout.addWidget(lLabel0)  # 0
        lLayout.setStretch(0, 1.5)
        # Select start date
        lLabel = QLabel()
        # lLabel.setText("Select start date: ")
        # lLabel.setFont(QFont("Microsoft YaHei", 20, QFont.Bold))
        lLabel.setAlignment(Qt.AlignCenter)
        lLayout.addWidget(lLabel)  # 1
        lLayout.setStretch(1, 1)
        # 日历
        self.calendar = QCalendarWidget()
        self.calendar.setGridVisible(True)
        self.calendar.setDateRange(QDate(1900, 1, 1), QDate(2195, 12, 24))
        self.calendar.setMaximumHeight(300)
        # lLayout.addWidget(self.calendar)  # 2
        # 速度控制控件
        self.speed_box = QWidget()
        self.speed_layout = QHBoxLayout()
        self.speed_box.setLayout(self.speed_layout)
        speed_button_width = 40  # 按钮宽度
        self.speedButton_moreSlower = QPushButton("<<<")
        self.speedButton_moreSlower.setMaximumWidth(speed_button_width)
        self.speedButton_slower = QPushButton("<")
        self.speedButton_slower.setMaximumWidth(speed_button_width)
        self.speedButton_origion = QPushButton("▶")
        self.speedButton_origion.setMaximumWidth(speed_button_width)
        self.speedButton_faster = QPushButton(">")
        self.speedButton_faster.setMaximumWidth(speed_button_width)
        self.speedButton_moreFaster = QPushButton(">>>")
        self.speedButton_moreFaster.setMaximumWidth(speed_button_width)
        self.speed_layout.addWidget(self.speedButton_moreSlower)
        self.speed_layout.addWidget(self.speedButton_slower)
        self.speed_layout.addWidget(self.speedButton_origion)
        self.speed_layout.addWidget(self.speedButton_faster)
        self.speed_layout.addWidget(self.speedButton_moreFaster)
        lLayout.addWidget(self.speed_box)
        lLayout.setStretch(2, 4)
        # 日历标签
        self.calendar_label = QLabel()
        self.calendar_label.setFont(QFont("Microsoft YaHei", 10, QFont.Bold))
        self.calendar_label.setAlignment(Qt.AlignCenter)
        self._onCalChanged()
        lLayout.addWidget(self.calendar_label)  # 3
        lLayout.setStretch(3, 1)
        # OK按钮
        self.okButton = QPushButton("▶")
        self.okButton.setMinimumHeight(50)
        self.okButton.setFont(QFont("Microsoft YaHei", 20, QFont.Bold))
        lLayout.addWidget(self.okButton)  # 4
        lLayout.setStretch(4, 1)
        # 空标签(为了占地方)
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

    def _init_planets(self):
        # 太阳
        sun_actor = self.dataProvider.build_sphere_actor('Sun', np.array([0, 0, 0]), 0.2, "Yellow")
        
        self.vtkWidget.renderer.AddActor(sun_actor)
        # 光照
        sun_light = vtkLight()
        sun_light.SetColor(1, 1, 1)
        sun_light.SetPosition(0, 0, 0)
        # self.vtkWidget.renderer.AddLight(sun_light)
        
        
        # 地球
        self.vtkWidget.renderer.AddActor(
            self.dataProvider.build_sphere_actor('Earth', np.array([0, 0, 0]), 0.09, "SkyBlue"))
        # 2016HO3
        self.vtkWidget.renderer.AddActor(
            self.dataProvider.build_sphere_actor('2016HO3', np.array([0, 1, 0]), 0.03, "Pink"))
        # Mercury
        self.vtkWidget.renderer.AddActor(
            self.dataProvider.build_sphere_actor('Mercury', np.array([0, 1, 0]), 0.05, "Gold"))
        # Venus
        self.vtkWidget.renderer.AddActor(
            self.dataProvider.build_sphere_actor('Venus', np.array([0, 1, 0]), 0.04, "Green"))
        # TODO 改半径，改为实际大小
        return None

    def _on_OKButton_clicked(self):
        log("OK clicked", "info")
        if not self.displaying:
            self.displaying = True
            self.okButton.setText("||")
            log("go, displaying=" + self.displaying.__str__(), "debug")
            self._display_loop()
        else:  # displaying
            self.displaying = False
            self.okButton.setText("▶")
            log("pause, displaying=" + self.displaying.__str__(), "debug")
        QApplication.processEvents()
        return

    def _display_loop(self):
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
        return

    def _refresh_vtkDisplay(self):
        self.vtkWidget.GetRenderWindow().Render()  # 刷新
        # Show current tdb
        self.vtkWidget.commentTextActor.SetInput("TDB=%.4f" % self.dataProvider.current_tdb)
        QApplication.processEvents()
        return

    def _onCalChanged(self):
        date = self.calendar.selectedDate()
        tdb = date2TDB(date)
        self.calendar_label.setText("%s  TDB = %.2f" % (str(date.toPyDate()), tdb))
        # todo

    # TODO 重新设置速率步长
    large_step = 0.3
    small_step = 0.1
    
    def _on_speedButton_moreSlower(self):
        step = self.large_step
        if self.delta_tdb - step > 0:
            self.delta_tdb -= step
        log("speed more slower, delta_tdb=%f" % self.delta_tdb, "info")

    def _on_speedButton_slower(self):
        step = self.small_step
        if self.delta_tdb - step > 0:
            self.delta_tdb -= step
        log("speed slower, delta_tdb=%f" % self.delta_tdb, "info")
    
    def _on_speedButton_origion(self):
        self.delta_tdb = self.origion_delta_tdb
        log("speed origion", "info")
    
    def _on_speedButton_faster(self):
        step = self.small_step
        self.delta_tdb += step
        log("speed faster, delta_tdb=%f" % self.delta_tdb, "info")

    def _on_speedButton_moreFaster(self):
        step = self.large_step
        self.delta_tdb += step
        log("speed more faster, delta_tdb=%f" % self.delta_tdb, "info")
    
    def closeEvent(self, *args, **kwargs):
        log("Window Closed.", "Info")
        exit()


# VTK控件对象
class MyVTKWidget(QVTKRenderWindowInteractor):
    def __init__(self):
        QVTKRenderWindowInteractor.__init__(self)
        self.renderer = vtk.vtkRenderer()
        self.GetRenderWindow().AddRenderer(self.renderer)
        self.interactor = self.GetRenderWindow().GetInteractor()
        # self.interactor.SetInteractorStyle(None)  # 禁用交互
        
        # 相机
        camera = vtk.vtkCamera()
        camera.SetViewAngle(40)
        camera.SetFocalPoint(0, 0, 0)
        camera.SetPosition(4, 0, 3)
        camera.SetViewUp(0, 0, 1)
        # print(camera)
        self.renderer.SetActiveCamera(camera)
        
        # 文字注记
        self.commentTextActor = vtkTextActor()
        self.commentTextActor.SetDisplayPosition(0, 0)
        self.commentTextActor.SetInput("Commont")
        self.renderer.AddActor(self.commentTextActor)


# 读取2016HO3的文本星历
class EphemerisReader:
    """ 读取星历文件 """
    def __init__(self, filename: str):
        self.filename = filename
        log("ephReader filename = %s" % self.filename, "debug")
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


# 为vtk管理、更新数据的部件
class DataProvider:
    def __init__(self, bsp_file="data/de430.bsp", eph_file="data/status.eph"):
        # data
        self._sphereSourceDic = {}  # 所有的行星的sphereSource
        self._bspReader = BspReader(bsp_file)
        self._ephReader = EphemerisReader(eph_file)
        self.scale = 1e-8  # 缩放尺度
        self.min_tdb, self.max_tdb = self._ephReader.min_tdb, self._ephReader.max_tdb
        self.current_tdb = 0

    def build_sphere_actor(self, key: str, Center: np.array(3), radius: float, color: str) -> vtk.vtkActor:
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


# 轨道积分器
class EquationBuilder:
    def __init__(self):
        self.P_Solar = 4.56e-6
        self.debug_info = []
        
        self.PARAM = {
            "Ref_Body": "Sun",
            "Ref_Frame": "J2000",
            "Surface_Area": 11.42,
            "Satellite_Mass": 666,
            "CR": 1.21,
            "Radiation": True,
            "Central_Body": "Sun",
            "Occulting_Bodies": [
                "Earth",
                "Moon",
                "Mercury",
                "Venus"
            ],
            "Perturbation_Bodies": [
                "Mercury Barycenter",
                "Venus Barycenter",
                "Earth Moon Barycenter",
                "Mars Barycenter",
                "Jupiter Barycenter",
                "Saturn Barycenter",
                "Uranus Barycenter",
                "Neptune Barycenter",
                "Pluto Barycenter",
            ]
            
        }
        
        # spice.furnsh(r"data\latest_leapseconds.tls.pc")
        # spice.furnsh(r"data\ORHM_______________00038.BSP")
        spice.furnsh(r"data\de430.bsp")
        # spice.furnsh(r"data\pck00010.tpc")
        # spice.furnsh(r"data\Gravity.tpc")
        # spice.furnsh(r"data\mar097.bsp")
        
    def perturbation(self, name, r_mex, t):
        """
        计算摄动天体加速度
        :self.PARAM name: 摄动天体名称
        :self.PARAM r_mex: 人造卫星（MEX）的位置向量
        :self.PARAM t: TDB时刻
        :return: 对应天体摄动加速度加速度[ax, ay, az]
        """
        r_body = spice.spkezr(name, t, self.PARAM["Ref_Frame"], 'None', self.PARAM["Ref_Body"])[0][:3]
        r_cent = spice.spkezr(self.PARAM["Central_Body"], t, self.PARAM["Ref_Frame"], 'None', self.PARAM["Ref_Body"])[0][:3]
        
        GM_body = spice.bodvrd(name, "GM", 1)[1][0]
        
        r_body_mex = r_mex - r_body
        r_cent_body = r_body - r_cent
        
        a = - GM_body * (r_cent_body / norm(r_cent_body) ** 3 + r_body_mex / norm(r_body_mex) ** 3)
        return a
    
    def agl_between(self, vec1, vec2):
        """
        计算两向量之间的夹角
        :self.PARAM vec1: 向量1
        :self.PARAM vec2: 向量2
        :return: 夹角弧度值
        """
        return acos(vec1 @ vec2 / norm(vec1) / norm(vec2))
    
    def body_shadow_function(self, r_mex, name, t):
        """
        计算阴影函数（shadow function）见英文教材P81 3.4.2
        :self.PARAM r_mex: 人造卫星（MEX）的位置向量
        :self.PARAM name: 遮挡天体的名称
        :self.PARAM t: TDB时刻
        :return: 阴影函数v
        """
        r_body = spice.spkezr(name, t, self.PARAM["Ref_Frame"], 'None', self.PARAM["Ref_Body"])[0][:3]
        r_sun = spice.spkezr("Sun", t, self.PARAM["Ref_Frame"], 'None', self.PARAM["Ref_Body"])[0][:3]
        r_body_mex = r_mex - r_body
        r_mex_sun = r_sun - r_mex
        
        R_body = spice.bodvrd(name, "RADII", 3)[1][0]
        RS = spice.bodvrd("Sun", "RADII", 3)[1][0]
        a = asin(RS / norm(r_mex_sun))
        b = asin(R_body / norm(r_body_mex))
        c = self.agl_between(-1 * r_body_mex, r_mex_sun)
        v = self.shadow_function(a, b, c)
        return v
    
    def shadow_function(self, a, b, c):
        """
        由天体的apparent radius（a、b）和两者间的apparent separation计算阴影函数
        见英文教材P82
        :self.PARAM a: 被遮挡天体（太阳）的apparent radius
        :self.PARAM b: 遮挡天体的apparent radius
        :self.PARAM c: 被遮挡天体和遮挡天体的 apparent separation
        :return: 阴影函数v
        """
        if abs(a - b) < c < a + b:
            x = (c ** 2 + a ** 2 - b ** 2) / 2 / c
            y = (a ** 2 - x ** 2) ** .5
            A = a ** 2 * acos(x / a) + b ** 2 * acos((c - x) / b) - c * y
            v = 1 - A / pi / a ** 2
        elif c >= a + b:
            # no occultation
            v = 1
        else:  # c <= abs(a - b)
            if a < b:
                # total occultation
                v = 0
            else:
                # partial but maximum
                v = 1 - b ** 2 / a ** 2
        return v
    
    def solar_radiation_pressure(self, t, y):
        """
        计算太阳辐射光压 见英文版教材P79 3.75
        :self.PARAM t: TDB时刻
        :self.PARAM y: 对应时刻人造卫星的状态向量[x, y, z, vx, vy, vz]
        :return: 太阳辐射光压加速度
        """
        r_mex = y[:3]
        r_sun = spice.spkezr("Sun", t, self.PARAM["Ref_Frame"], 'None', self.PARAM["Ref_Body"])[0][:3]
        r_sun_mex = r_mex - r_sun
        
        AU_km = spice.gdpool("AU", 0, 1)[0]
        v = min([self.body_shadow_function(r_mex, name, t) for name in self.PARAM["Occulting_Bodies"]])
        CR, A, m = [self.PARAM[i] for i in ["CR", "Surface_Area", "Satellite_Mass"]]
        a = v * self.P_Solar * CR * A / m * AU_km ** 2 * r_sun_mex / norm(r_sun_mex) ** 3
        
        return a * 1e-3  # km * s^-2
        pass
    
    def n_body_equation(self, t, y):
        """
        n体问题微分方程 y_dot = f(t, y)
        见英文版教材P117
        :self.PARAM t: 对应tdb时刻
        :self.PARAM y: 对应时刻人造卫星的状态向量[x, y, z, vx, vy, vz]
        :return: y_dot
        """
        y_dot = np.empty((6,))
        y_dot[:3] = y[3:]
        
        r_mex = y[:3]
        r_central = spice.spkezr(self.PARAM["Central_Body"], t, self.PARAM["Ref_Frame"], 'None', self.PARAM["Ref_Body"])[0][:3]
        r_central_mex = r_mex - r_central
        GM_central = spice.bodvrd(self.PARAM["Central_Body"], "GM", 1)[1][0]
        a_central = - GM_central * r_central_mex / norm(r_central_mex) ** 3
        perturbations = sum([self.perturbation(name, r_mex, t) for name in self.PARAM["Perturbation_Bodies"]])
        y_dot[3:] = a_central + perturbations
        
        if self.PARAM["Radiation"]:
            radiation_pressure = self.solar_radiation_pressure(t, y)
            y_dot[3:] += radiation_pressure
        
        return y_dot  # v, a


def date2TDB(date):
        # 2458551.000000000 = A.D. 2019-Mar-08
        defaultTDB = 2458551
        defaultDate = QDate(2019, 3, 8)
        tdb = defaultTDB + defaultDate.daysTo(date)
        return tdb


# 输出日志的函数
def log(message: str, level: str):
    level = level.lower()
    if level == "debug":
        print("[Debug] %s" % message)
    elif level == "info":
        print("\033[0;36m[Info]  %s\033[0m" % message)
    elif level == "warn":
        print("\033[0;34m[Warn]  %s\033[0m" % message)
    elif level == "error":
        print("\033[0;31m[error] %s\033[0m" % message)
    else:
        print("\033[0;31m[error] Log level error\033[0m")


if __name__ == "__main__":
    app = QApplication(sys.argv)

    win = MainWindow()

    win.show()

    log("warn", "warn")
    log("error", "error")
    
    # 开始运行
    win.okButton.click()
    
    win.vtkWidget.interactor.Initialize()

    sys.exit(app.exec_())
