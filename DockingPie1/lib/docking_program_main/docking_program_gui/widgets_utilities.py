from pymol.Qt import QtWidgets, QtCore, QtGui
from pymol import Qt

def grid_dimensions(on_value_changed_function):

    spacing = QtWidgets.QLabel("All:")
    spacing_scroll_vis = QtWidgets.QDoubleSpinBox()
    spacing_scroll_vis.setRange(0, 10)
    spacing_scroll_vis.setSingleStep(0.01)
    spacing_scroll_vis.setDecimals(3)
    spacing_scroll_vis.setValue(0.375)
    spacing_scroll_vis.valueChanged.connect(on_value_changed_function)
    #
    x = QtWidgets.QLabel("X:")
    x_scroll_vis = QtWidgets.QSpinBox()
    x_scroll_vis.setRange(1, 50)
    x_scroll_vis.setSingleStep(1)
    x_scroll_vis.setValue(4)
    x_scroll_vis.valueChanged.connect(on_value_changed_function)

    y = QtWidgets.QLabel("Y:")
    y_scroll_vis = QtWidgets.QSpinBox()
    y_scroll_vis.setRange(1, 50)
    y_scroll_vis.setSingleStep(1)
    y_scroll_vis.setValue(4)
    y_scroll_vis.valueChanged.connect(on_value_changed_function)

    z = QtWidgets.QLabel("Z:")
    z_scroll_vis = QtWidgets.QSpinBox()
    z_scroll_vis.setRange(1, 50)
    z_scroll_vis.setSingleStep(1)
    z_scroll_vis.setValue(4)
    z_scroll_vis.valueChanged.connect(on_value_changed_function)

    return spacing, spacing_scroll_vis, x, x_scroll_vis, y, y_scroll_vis, z, z_scroll_vis


def grid_position(on_value_changed_function):

    x = QtWidgets.QLabel("X:")
    x_scroll = QtWidgets.QDoubleSpinBox()
    x_scroll.setRange(-1000, 1000)
    x_scroll.setSingleStep(00.10)
    x_scroll.valueChanged.connect(on_value_changed_function)

    y = QtWidgets.QLabel("Y:")
    y_scroll = QtWidgets.QDoubleSpinBox()
    y_scroll.setRange(-1000, 1000)
    y_scroll.setSingleStep(00.10)
    y_scroll.valueChanged.connect(on_value_changed_function)

    z = QtWidgets.QLabel("Z:")
    z_scroll = QtWidgets.QDoubleSpinBox()
    z.setBuddy(z_scroll)
    z_scroll.setRange(-1000, 1000)
    z_scroll.setSingleStep(00.10)
    z_scroll.valueChanged.connect(on_value_changed_function)

    return x, x_scroll, y, y_scroll, z, z_scroll
