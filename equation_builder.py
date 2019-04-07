from math import asin, acos, pi, atan, sin, cos, atan2
import numpy as np
import pandas as pd
import spiceypy as spice
from numpy.linalg import norm
from scipy.integrate import solve_ivp


P_Solar = 4.56e-6
debug_info = []

PARAM = {
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


def perturbation(name, r_mex, t):
    """
    计算摄动天体加速度
    :param name: 摄动天体名称
    :param r_mex: 人造卫星（MEX）的位置向量
    :param t: TDB时刻
    :return: 对应天体摄动加速度加速度[ax, ay, az]
    """
    r_body = spice.spkezr(name, t, PARAM["Ref_Frame"], 'None', PARAM["Ref_Body"])[0][:3]
    r_cent = spice.spkezr(PARAM["Central_Body"], t, PARAM["Ref_Frame"], 'None', PARAM["Ref_Body"])[0][:3]

    GM_body = spice.bodvrd(name, "GM", 1)[1][0]

    r_body_mex = r_mex - r_body
    r_cent_body = r_body - r_cent

    a = - GM_body * (r_cent_body / norm(r_cent_body) ** 3 + r_body_mex / norm(r_body_mex) ** 3)
    return a


def agl_between(vec1, vec2):
    """
    计算两向量之间的夹角
    :param vec1: 向量1
    :param vec2: 向量2
    :return: 夹角弧度值
    """
    return acos(vec1 @ vec2 / norm(vec1) / norm(vec2))


def body_shadow_function(r_mex, name, t):
    """
    计算阴影函数（shadow function）见英文教材P81 3.4.2
    :param r_mex: 人造卫星（MEX）的位置向量
    :param name: 遮挡天体的名称
    :param t: TDB时刻
    :return: 阴影函数v
    """
    r_body = spice.spkezr(name, t, PARAM["Ref_Frame"], 'None', PARAM["Ref_Body"])[0][:3]
    r_sun = spice.spkezr("Sun", t, PARAM["Ref_Frame"], 'None', PARAM["Ref_Body"])[0][:3]
    r_body_mex = r_mex - r_body
    r_mex_sun = r_sun - r_mex

    R_body = spice.bodvrd(name, "RADII", 3)[1][0]
    RS = spice.bodvrd("Sun", "RADII", 3)[1][0]
    a = asin(RS / norm(r_mex_sun))
    b = asin(R_body / norm(r_body_mex))
    c = agl_between(-1 * r_body_mex, r_mex_sun)
    v = shadow_function(a, b, c)
    return v


def shadow_function(a, b, c):
    """
    由天体的apparent radius（a、b）和两者间的apparent separation计算阴影函数
    见英文教材P82
    :param a: 被遮挡天体（太阳）的apparent radius
    :param b: 遮挡天体的apparent radius
    :param c: 被遮挡天体和遮挡天体的 apparent separation
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


def solar_radiation_pressure(t, y):
    """
    计算太阳辐射光压 见英文版教材P79 3.75
    :param t: TDB时刻
    :param y: 对应时刻人造卫星的状态向量[x, y, z, vx, vy, vz]
    :return: 太阳辐射光压加速度
    """
    r_mex = y[:3]
    r_sun = spice.spkezr("Sun", t, PARAM["Ref_Frame"], 'None', PARAM["Ref_Body"])[0][:3]
    r_sun_mex = r_mex - r_sun

    AU_km = spice.gdpool("AU", 0, 1)[0]
    v = min([body_shadow_function(r_mex, name, t) for name in PARAM["Occulting_Bodies"]])
    CR, A, m = [PARAM[i] for i in ["CR", "Surface_Area", "Satellite_Mass"]]
    a = v * P_Solar * CR * A / m * AU_km ** 2 * r_sun_mex / norm(r_sun_mex) ** 3

    return a * 1e-3  # km * s^-2
    pass


def n_body_equation(t, y):
    """
    n体问题微分方程 y_dot = f(t, y)
    见英文版教材P117
    :param t: 对应tdb时刻
    :param y: 对应时刻人造卫星的状态向量[x, y, z, vx, vy, vz]
    :return: y_dot
    """
    y_dot = np.empty((6,))
    y_dot[:3] = y[3:]

    r_mex = y[:3]
    r_central = spice.spkezr(PARAM["Central_Body"], t, PARAM["Ref_Frame"], 'None', PARAM["Ref_Body"])[0][:3]
    r_central_mex = r_mex - r_central
    GM_central = spice.bodvrd(PARAM["Central_Body"], "GM", 1)[1][0]
    a_central = - GM_central * r_central_mex / norm(r_central_mex) ** 3
    perturbations = sum([perturbation(name, r_mex, t) for name in PARAM["Perturbation_Bodies"]])
    y_dot[3:] = a_central + perturbations

    if PARAM["Radiation"]:
        radiation_pressure = solar_radiation_pressure(t, y)
        y_dot[3:] += radiation_pressure

    return y_dot  # v, a


if __name__ == '__main__':
    # print(debug_info)
    # TODO: Different Ref_Body of MEX lead to different error ?
    # PARAM['Ref_Body'] = 'SSB'
    start = '2003-06-29 19:00:00'
    end = '2003-07-02 12:00:00'
    freq = '1H'

    t0 = spice.str2et(start)
    y0 = spice.spkezr("MEX", t0, PARAM["Ref_Frame"], 'None', PARAM["Ref_Body"])[0]
    time_series = pd.date_range(start, end, freq=freq).to_pydatetime()
    ts = np.array([spice.str2et(str(t)) for t in time_series])
    result = solve_ivp(n_body_equation, (ts[0], ts[-1]), y0, method='BDF',
                       t_eval=ts, rtol=1e-13, atol=1e-16)
    solution = result.y.T
    print(solution)
    pass
