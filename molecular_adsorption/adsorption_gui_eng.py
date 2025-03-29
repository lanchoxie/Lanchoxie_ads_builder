#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
分子吸附模拟可视化界面
"""

import os
import sys
import json
import time
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                            QLabel, QPushButton, QFileDialog, QComboBox, QSpinBox, QDoubleSpinBox,
                            QTabWidget, QMessageBox, QGroupBox, QRadioButton, QSlider, QComboBox,
                            QLineEdit, QCheckBox, QSplitter, QFrame, QDialog, QSpinBox,
                            QSizePolicy, QListWidget, QAbstractItemView, QButtonGroup,
                            QDialogButtonBox, QInputDialog, QProgressDialog, QAction)

from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QFont, QIcon

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.path import Path
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation
from mpl_toolkits.mplot3d import proj3d
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

from .languages import TRANSLATIONS, LanguageManager

from ase import Atoms
from ase.io import read, write
from ase.visualize.plot import plot_atoms
from ase.data import covalent_radii
from ase.data.colors import jmol_colors

# 导入我们的模块
# 当直接运行时
if __name__ == "__main__":
    # 添加父目录到模块搜索路径
    import os
    import sys
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, parent_dir)
    from molecular_adsorption.position_molecule import (find_surface_atoms,  
                               is_molecule_inside_substrate_xy, find_bridge_sites,
                               find_hollow_sites, validate_system, position_molecule_on_site)
# 当作为模块导入时
else:
    try:
        # 当作为独立脚本运行时
        from molecular_adsorption.position_molecule import (find_surface_atoms,  
                                   is_molecule_inside_substrate_xy, find_bridge_sites,
                                   find_hollow_sites, validate_system, position_molecule_on_site)
    # 当作为模块导入时
    except ImportError:
        # 尝试绝对导入
        try:
            from position_molecule import (find_surface_atoms, 
                                       is_molecule_inside_substrate_xy, find_bridge_sites,
                                       find_hollow_sites, validate_system, position_molecule_on_site)
        # 如果绝对导入失败，回退到相对导入
        except ImportError:
            from .position_molecule import (find_surface_atoms, 
                               is_molecule_inside_substrate_xy, find_bridge_sites,
                               find_hollow_sites, validate_system, position_molecule_on_site)
        
# 定义标准元素颜色表（全局变量）
ELEMENT_COLORS = {
    'H': '#FFFFFF',  # 白色
    'He': '#D9FFFF',
    'Li': '#CC80FF',
    'Be': '#C2FF00',
    'B': '#FFB5B5',
    'C': '#909090',  # 灰色
    'N': '#3050F8',  # 蓝色
    'O': '#FF0D0D',  # 红色
    'F': '#90E050',
    'Ne': '#B3E3F5',
    'Na': '#AB5CF2',
    'Mg': '#8AFF00',
    'Al': '#BFA6A6',
    'Si': '#F0C8A0',
    'P': '#FF8000',  # 橙色
    'S': '#FFFF30',  # 黄色
    'Cl': '#1FF01F',  # 绿色
    'Ar': '#80D1E3',
    'K': '#8F40D4',
    'Ca': '#3DFF00',
    'Sc': '#E6E6E6',
    'Ti': '#BFC2C7',
    'V': '#A6A6AB',
    'Cr': '#8A99C7',
    'Mn': '#9C7AC7',
    'Fe': '#E06633',  # 棕色
    'Co': '#F090A0',
    'Ni': '#50D050',
    'Cu': '#C88033',  # 铜色
    'Zn': '#7D80B0',
    'Ga': '#C28F8F',
    'Ge': '#668F8F',
    'As': '#BD80E3',
    'Se': '#FFA100',
    'Br': '#A62929',
    'Kr': '#5CB8D1',
    'Rb': '#702EB0',
    'Sr': '#00FF00',
    'Y': '#94FFFF',
    'Zr': '#94E0E0',
    'Nb': '#73C2C9',
    'Mo': '#54B5B5',
    'Tc': '#3B9E9E',
    'Ru': '#248F8F',
    'Rh': '#0A7D8C',
    'Pd': '#006985',
    'Ag': '#C0C0C0',  # 银色
    'Cd': '#FFD98F',
    'In': '#A67573',
    'Sn': '#668080',
    'Sb': '#9E63B5',
    'Te': '#D47A00',
    'I': '#940094',
    'Xe': '#429EB0',
    'Cs': '#57178F',
    'Ba': '#00C900',
    'La': '#70D4FF',
    'Ce': '#FFFFC7',
    'Pr': '#D9FFC7',
    'Nd': '#C7FFC7',
    'Pm': '#A3FFC7',
    'Sm': '#8FFFC7',
    'Eu': '#61FFC7',
    'Gd': '#45FFC7',
    'Tb': '#30FFC7',
    'Dy': '#1FFFC7',
    'Ho': '#00FF9C',
    'Er': '#00E675',
    'Tm': '#00D452',
    'Yb': '#00BF38',
    'Lu': '#00AB24',
    'Hf': '#4DC2FF',
    'Ta': '#4DA6FF',
    'W': '#2194D6',
    'Re': '#267DAB',
    'Os': '#266696',
    'Ir': '#175487',
    'Pt': '#D0D0E0',
    'Au': '#FFD123',  # 金色
    'Hg': '#B8B8D0',
    'Tl': '#A6544D',
    'Pb': '#575961',
    'Bi': '#9E4FB5',
    'Po': '#AB5C00',
    'At': '#754F45',
    'Rn': '#428296',
    'Fr': '#420066',
    'Ra': '#007D00',
    'Ac': '#70ABFA',
    'Th': '#00BAFF',
    'Pa': '#00A1FF',
    'U': '#008FFF',
    'Np': '#0080FF',
    'Pu': '#006BFF',
    'Am': '#545CF2',
    'Cm': '#785CE3',
    'Bk': '#8A4FE3',
    'Cf': '#A136D4',
    'Es': '#B31FD4',
    'Fm': '#B31FBA',
    'Md': '#B30DA6',
    'No': '#BD0D87',
    'Lr': '#C70066',
    'Rf': '#CC0059',
    'Db': '#D1004F',
    'Sg': '#D90045',
    'Bh': '#E00038',
    'Hs': '#E6002E',
    'Mt': '#EB0026'
}

def is_molecule_pbc_valid(molecule, substrate, projection_type='xy'):
    """
    判断分子在周期性边界条件下的状态
    
    Args:
        molecule: 分子对象
        substrate: 基底对象
        projection_type: 投影类型 ('xy', 'xz', 'yz')
    
    Returns:
        tuple: (状态, 颜色, 镜像透明度)
        状态: 'valid' (绿色), 'pbc' (蓝色), 'overlap' (红色)
    """
    if substrate is None:
        return 'valid', 'green', 0.3
    
    # 获取投影方向的索引
    if projection_type == 'xy':
        x_idx, y_idx = 0, 1
        cell_vectors = [substrate.get_cell()[0, [x_idx, y_idx]], 
                       substrate.get_cell()[1, [x_idx, y_idx]]]
    elif projection_type == 'xz':
        x_idx, y_idx = 0, 2
        cell_vectors = [substrate.get_cell()[0, [x_idx, y_idx]], 
                       substrate.get_cell()[2, [x_idx, y_idx]]]
    else:  # yz
        x_idx, y_idx = 1, 2
        cell_vectors = [substrate.get_cell()[1, [x_idx, y_idx]], 
                       substrate.get_cell()[2, [x_idx, y_idx]]]
    
    # 获取分子在投影平面的点
    molecule_points = molecule.get_positions()[:, [x_idx, y_idx]]
    molecule_min_x, molecule_min_y = molecule_points.min(axis=0)
    molecule_max_x, molecule_max_y = molecule_points.max(axis=0)
    
    # 只在 XZ 和 YZ 投影中检查分子和基底的重叠
    if projection_type in ['xz', 'yz']:
        # 获取基底在投影平面的点
        substrate_points = substrate.get_positions()[:, [x_idx, y_idx]]
        substrate_min_x, substrate_min_y = substrate_points.min(axis=0)
        substrate_max_x, substrate_max_y = substrate_points.max(axis=0)
        
        # 检查分子和基底是否重叠
        if (molecule_min_x < substrate_max_x and molecule_max_x > substrate_min_x and
            molecule_min_y < substrate_max_y and molecule_max_y > substrate_min_y):
            return 'overlap', 'red', 0.3
    
    # 检查是否跨越边界
    directions = []
    
    # 检查第一个方向（X或Y）
    if molecule_min_x < 0:
        directions.append(np.array([1, 0]))
    if molecule_max_x > cell_vectors[0][0]:
        directions.append(np.array([-1, 0]))
    
    # 检查第二个方向（Y或Z）
    if molecule_min_y < 0:
        directions.append(np.array([0, 1]))
    if molecule_max_y > cell_vectors[1][1]:
        directions.append(np.array([0, -1]))
    
    # 检查对角方向
    if (molecule_min_x < 0 and molecule_min_y < 0) or (molecule_max_x > cell_vectors[0][0] and molecule_max_y > cell_vectors[1][1]):
        directions.append(np.array([1, 1]))
    if (molecule_min_x < 0 and molecule_max_y > cell_vectors[1][1]) or (molecule_max_x > cell_vectors[0][0] and molecule_min_y < 0):
        directions.append(np.array([1, -1]))
    
    # 如果没有跨越边界，返回绿色
    if len(directions) == 0:
        return 'valid', 'green', 0.3
    
    # 检查是否有重叠
    has_overlap = False
    for direction in directions:
        offset = (direction[0] * cell_vectors[0] + 
                 direction[1] * cell_vectors[1])
        mirror_points = molecule_points + offset
        
        # 检查原始分子和镜像分子是否重叠
        if (molecule_min_x < molecule_max_x and molecule_min_y < molecule_max_y and
            np.any((mirror_points[:, 0] >= molecule_min_x) & (mirror_points[:, 0] <= molecule_max_x) &
                  (mirror_points[:, 1] >= molecule_min_y) & (mirror_points[:, 1] <= molecule_max_y))):
            has_overlap = True
            break
    
    if has_overlap:
        return 'overlap', 'red', 0.3
    else:
        return 'pbc', 'blue', 0.3


# 默认颜色（用于未知元素）
DEFAULT_COLOR = '#777777'  # 中灰色

class MplCanvas(FigureCanvas):
    """Matplotlib画布用于3D显示原子结构"""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        plt.rcParams['axes.facecolor'] = 'white'  # 设置全局背景颜色
        
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.patch.set_facecolor('white')  # 设置图表背景为白色
        self.axes = self.fig.add_subplot(111, projection='3d')
        self.axes.set_facecolor('white')  # 设置坐标轴背景为白色
        
        # 确保坐标轴有标签
        self.axes.set_xlabel('X (Å)')
        self.axes.set_ylabel('Y (Å)')
        self.axes.set_zlabel('Z (Å)')
        
        # 设置初始视角
        self.axes.view_init(elev=30, azim=30)
        
        super(MplCanvas, self).__init__(self.fig)
        self.setParent(parent)  # 确保设置了父组件
        
        FigureCanvas.setSizePolicy(self,
                                  QSizePolicy.Expanding,
                                  QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        
        # 立即绘制以显示白色背景
        self.draw()


class MplCanvas2D(FigureCanvas):
    """2D Matplotlib canvas"""
    
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi, tight_layout=True)
        super().__init__(self.fig)
        self.setParent(parent)
        self.axes = self.fig.add_subplot(111)


    def plot_projection(self, substrate, molecule, check_position=True):
        """Plot 2D projection of substrate and molecule"""
        if not substrate and not molecule:
            self.clear_plot()
            return
        
        try:
            self.axes.clear()
            
            # 根据投影类型设置坐标轴标签
            if self.projection_type == 'xy':
                x_label, y_label = 'X (Å)', 'Y (Å)'
                x_idx, y_idx = 0, 1
            elif self.projection_type == 'xz':
                x_label, y_label = 'X (Å)', 'Z (Å)'
                x_idx, y_idx = 0, 2
            else:  # yz
                x_label, y_label = 'Y (Å)', 'Z (Å)'
                x_idx, y_idx = 1, 2
            
            # 显示坐标轴和标签
            self.axes.set_frame_on(True)
            self.axes.spines['top'].set_visible(False)
            self.axes.spines['right'].set_visible(False)
            self.axes.set_xlabel(x_label, fontsize=8)
            self.axes.set_ylabel(y_label, fontsize=8)
            
            # 设置等比例显示
            self.axes.set_aspect('equal')
            
            # 绘制基底投影（考虑周期性）
            if substrate is not None:
                cell = substrate.get_cell()
                
                # 获取周期性边界
                if self.projection_type == 'xy':
                    # 绘制晶格边界（虚线）
                    x_vec = cell[0, [x_idx, y_idx]]
                    y_vec = cell[1, [x_idx, y_idx]]
                    corners = np.array([
                        [0, 0],
                        x_vec,
                        x_vec + y_vec,
                        y_vec,
                        [0, 0]
                    ])
                elif self.projection_type == 'xz':
                    # 绘制晶格边界（虚线）
                    x_vec = cell[0, [x_idx, y_idx]]
                    z_vec = cell[2, [x_idx, y_idx]]
                    corners = np.array([
                        [0, 0],
                        x_vec,
                        x_vec + z_vec,
                        z_vec,
                        [0, 0]
                    ])
                else:  # yz
                    # 绘制晶格边界（虚线）
                    y_vec = cell[1, [x_idx, y_idx]]
                    z_vec = cell[2, [x_idx, y_idx]]
                    corners = np.array([
                        [0, 0],
                        y_vec,
                        y_vec + z_vec,
                        z_vec,
                        [0, 0]
                    ])
                
                # 绘制晶格边界
                self.axes.plot(corners[:, 0], corners[:, 1],
                             linestyle='--', color='gray', alpha=0.5)
                
                # 绘制基底
                positions = substrate.get_positions()
                points = positions[:, [x_idx, y_idx]]
                try:
                    hull = ConvexHull(points)
                    hull_points = points[hull.vertices]
                    hull_points = np.vstack([hull_points, hull_points[0]])
                    self.axes.fill(hull_points[:, 0], hull_points[:, 1],
                                 alpha=0.2, color='blue', label=self.lang.get_text('substrate'))
                except:
                    min_x, min_y = points.min(axis=0)
                    max_x, max_y = points.max(axis=0)
                    box = np.array([[min_x, min_y],
                                  [max_x, min_y],
                                  [max_x, max_y],
                                  [min_x, max_y],
                                  [min_x, min_y]])
                    self.axes.fill(box[:, 0], box[:, 1],
                                 alpha=0.2, color='blue', label=self.lang.get_text('substrate'))
            
            # 绘制分子投影（考虑周期性）
            if molecule is not None:
                positions = molecule.get_positions()
                points = positions[:, [x_idx, y_idx]]
                
                # 更新分子状态和颜色
                if self.projection_type == 'xy':
                    # 获取 XY 平面的状态
                    xy_status, xy_color, xy_alpha = is_molecule_pbc_valid(molecule, substrate, 'xy')
                    ProjectionCanvas.molecule_status = xy_status
                    ProjectionCanvas.molecule_color = xy_color
                    ProjectionCanvas.molecule_mirror_alpha = xy_alpha
                elif self.projection_type == 'xz':
                    # 获取 XZ 平面的状态
                    xz_status, xz_color, xz_alpha = is_molecule_pbc_valid(molecule, substrate, 'xz')
                    # 如果 XZ 平面显示红色（重叠），则覆盖 XY 平面的状态
                    if xz_status == 'overlap':
                        ProjectionCanvas.molecule_status = xz_status
                        ProjectionCanvas.molecule_color = xz_color
                        ProjectionCanvas.molecule_mirror_alpha = xz_alpha
                else:  # yz
                    # 获取 YZ 平面的状态
                    yz_status, yz_color, yz_alpha = is_molecule_pbc_valid(molecule, substrate, 'yz')
                    # 如果 YZ 平面显示红色（重叠），则覆盖 XY 平面的状态
                    if yz_status == 'overlap':
                        ProjectionCanvas.molecule_status = yz_status
                        ProjectionCanvas.molecule_color = yz_color
                        ProjectionCanvas.molecule_mirror_alpha = yz_alpha
                
                # 使用类变量中的颜色
                color = ProjectionCanvas.molecule_color
                mirror_alpha = ProjectionCanvas.molecule_mirror_alpha
                
                # 绘制原始分子
                try:
                    hull = ConvexHull(points)
                    hull_points = points[hull.vertices]
                    hull_points = np.vstack([hull_points, hull_points[0]])
                    self.axes.fill(hull_points[:, 0], hull_points[:, 1],
                                 alpha=0.5, color=color, label=self.lang.get_text('molecule'))
                except:
                    min_x, min_y = points.min(axis=0)
                    max_x, max_y = points.max(axis=0)
                    box = np.array([[min_x, min_y],
                                  [max_x, min_y],
                                  [max_x, max_y],
                                  [min_x, max_y],
                                  [min_x, min_y]])
                    self.axes.fill(box[:, 0], box[:, 1],
                                 alpha=0.5, color=color, label=self.lang.get_text('molecule'))
                
                # 处理周期性边界条件
                if substrate is not None:
                    cell = substrate.get_cell()
                    if self.projection_type == 'xy':
                        cell_vectors = [cell[0, [x_idx, y_idx]], cell[1, [x_idx, y_idx]]]
                    elif self.projection_type == 'xz':
                        cell_vectors = [cell[0, [x_idx, y_idx]], cell[2, [x_idx, y_idx]]]
                    else:  # yz
                        cell_vectors = [cell[1, [x_idx, y_idx]], cell[2, [x_idx, y_idx]]]
                    
                    # 检查分子是否跨越边界
                    min_x, min_y = points.min(axis=0)
                    max_x, max_y = points.max(axis=0)
                    
                    # 定义需要复制的方向
                    directions = []
                    
                    # 检查第一个方向（X或Y）
                    if min_x < 0:
                        directions.append(np.array([1, 0]))
                    if max_x > cell_vectors[0][0]:
                        directions.append(np.array([-1, 0]))
                    
                    # 检查第二个方向（Y或Z）
                    if min_y < 0:
                        directions.append(np.array([0, 1]))
                    if max_y > cell_vectors[1][1]:
                        directions.append(np.array([0, -1]))
                    
                    # 检查对角方向
                    if (min_x < 0 and min_y < 0) or (max_x > cell_vectors[0][0] and max_y > cell_vectors[1][1]):
                        directions.append(np.array([1, 1]))
                    if (min_x < 0 and max_y > cell_vectors[1][1]) or (max_x > cell_vectors[0][0] and min_y < 0):
                        directions.append(np.array([1, -1]))
                    
                    # 绘制镜像分子
                    for direction in directions:
                        # 计算偏移量
                        offset = (direction[0] * cell_vectors[0] + 
                                direction[1] * cell_vectors[1])
                        
                        # 复制并平移分子
                        mirror_points = points + offset
                        
                        try:
                            hull = ConvexHull(mirror_points)
                            hull_points = mirror_points[hull.vertices]
                            hull_points = np.vstack([hull_points, hull_points[0]])
                            self.axes.fill(hull_points[:, 0], hull_points[:, 1],
                                         alpha=mirror_alpha, color=color,
                                         label=f'{self.lang.get_text("molecule")}({self.lang.get_text("mirror")})')
                        except:
                            min_x, min_y = mirror_points.min(axis=0)
                            max_x, max_y = mirror_points.max(axis=0)
                            box = np.array([[min_x, min_y],
                                          [max_x, min_y],
                                          [max_x, max_y],
                                          [min_x, max_y],
                                          [min_x, min_y]])
                            self.axes.fill(box[:, 0], box[:, 1],
                                         alpha=mirror_alpha, color=color,
                                         label=f'{self.lang.get_text("molecule")}({self.lang.get_text("mirror")})')
            
            # 标记锚点
            if hasattr(self, 'anchor_atom_index') and self.anchor_atom_index is not None and molecule is not None:
                anchor_pos = molecule.positions[self.anchor_atom_index]
                self.axes.scatter(anchor_pos[x_idx], anchor_pos[y_idx], 
                                color='yellow', marker='o', s=100, label=self.lang.get_text('anchor_point'))
                self.axes.annotate(f'({anchor_pos[x_idx]:.1f}, {anchor_pos[y_idx]:.1f})',
                                xy=(anchor_pos[x_idx], anchor_pos[y_idx]),
                                xytext=(5, 5), textcoords='offset points', fontsize=8)
            
            # 标记吸附位点
            if hasattr(self, 'site_position') and self.site_position is not None:
                self.axes.scatter(self.site_position[x_idx], self.site_position[y_idx],
                                color='red', marker='*', s=100, label=self.lang.get_text('adsorption_site'))
                self.axes.annotate(f'({self.site_position[x_idx]:.1f}, {self.site_position[y_idx]:.1f})',
                                xy=(self.site_position[x_idx], self.site_position[y_idx]),
                                xytext=(5, 5), textcoords='offset points', fontsize=8)
            
            # 只在顶部画布显示图例
            if self.with_legend:
                # 获取高度信息
                height_text = ""
                if substrate is not None:
                    substrate_height = substrate.get_cell()[2, 2]
                    height_text = f"{self.lang.get_text('lattice_height')}: {substrate_height:.1f} Å"
                if molecule is not None:
                    molecule_height = np.max(molecule.positions[:, 2])
                    if height_text:
                        height_text += f"\n{self.lang.get_text('molecule_highest_point')}: {molecule_height:.1f} Å"
                    else:
                        height_text = f"{self.lang.get_text('molecule_highest_point')}: {molecule_height:.1f} Å"
                
                # 添加高度信息到图例
                if height_text:
                    height_patch = mpatches.Patch(color='none', label=height_text)
                    handles, labels = self.axes.get_legend_handles_labels()
                    handles.append(height_patch)
                    labels.append(height_text)
                    self.axes.legend(handles, labels, loc='upper center',
                                bbox_to_anchor=(0.5, 1.2), ncol=2,
                                fontsize='small', framealpha=0.8)
            
            # 设置相等的纵横比
            self.axes.set_aspect('equal')
            
            # 调整布局
            self.figure.tight_layout()
            if self.with_legend:
                self.figure.subplots_adjust(top=0.8)  # 为图例留出空间
            
            self.draw()
            
        except Exception as e:
            print(f"Error plotting {self.projection_type} projection: {str(e)}")
            import traceback
            traceback.print_exc()

class ProjectionCanvas(MplCanvas2D):
    # 添加类变量来存储分子的状态和颜色
    molecule_status = None
    molecule_color = None
    molecule_mirror_alpha = None
    
    def __init__(self, parent=None, projection_type="XY", width=10, height=4, dpi=100, with_legend=False):
        super().__init__(parent, width, height, dpi)
        self.projection_type = projection_type.lower()
        self.with_legend = with_legend
        self.anchor_atom_index = None
        self.site_position = None
        # 获取语言管理器实例
        self.lang = LanguageManager()
    
    def plot_projection(self, substrate, molecule, check_position=True):
        """Plot 2D projection of substrate and molecule"""
        if not substrate and not molecule:
            self.clear_plot()
            return
        
        try:
            self.axes.clear()
            
            # 根据投影类型设置坐标轴标签
            if self.projection_type == 'xy':
                x_label, y_label = 'X (Å)', 'Y (Å)'
                x_idx, y_idx = 0, 1
            elif self.projection_type == 'xz':
                x_label, y_label = 'X (Å)', 'Z (Å)'
                x_idx, y_idx = 0, 2
            else:  # yz
                x_label, y_label = 'Y (Å)', 'Z (Å)'
                x_idx, y_idx = 1, 2
            
            # 显示坐标轴和标签
            self.axes.set_frame_on(True)
            self.axes.spines['top'].set_visible(False)
            self.axes.spines['right'].set_visible(False)
            self.axes.set_xlabel(x_label, fontsize=8)
            self.axes.set_ylabel(y_label, fontsize=8)
            
            # 设置等比例显示
            self.axes.set_aspect('equal')
            
            # 绘制基底投影（考虑周期性）
            if substrate is not None:
                cell = substrate.get_cell()
                
                # 获取周期性边界
                if self.projection_type == 'xy':
                    # 绘制晶格边界（虚线）
                    x_vec = cell[0, [x_idx, y_idx]]
                    y_vec = cell[1, [x_idx, y_idx]]
                    corners = np.array([
                        [0, 0],
                        x_vec,
                        x_vec + y_vec,
                        y_vec,
                        [0, 0]
                    ])
                elif self.projection_type == 'xz':
                    # 绘制晶格边界（虚线）
                    x_vec = cell[0, [x_idx, y_idx]]
                    z_vec = cell[2, [x_idx, y_idx]]
                    corners = np.array([
                        [0, 0],
                        x_vec,
                        x_vec + z_vec,
                        z_vec,
                        [0, 0]
                    ])
                else:  # yz
                    # 绘制晶格边界（虚线）
                    y_vec = cell[1, [x_idx, y_idx]]
                    z_vec = cell[2, [x_idx, y_idx]]
                    corners = np.array([
                        [0, 0],
                        y_vec,
                        y_vec + z_vec,
                        z_vec,
                        [0, 0]
                    ])
                
                # 绘制晶格边界
                self.axes.plot(corners[:, 0], corners[:, 1],
                             linestyle='--', color='gray', alpha=0.5)
                
                # 绘制基底
                positions = substrate.get_positions()
                points = positions[:, [x_idx, y_idx]]
                try:
                    hull = ConvexHull(points)
                    hull_points = points[hull.vertices]
                    hull_points = np.vstack([hull_points, hull_points[0]])
                    self.axes.fill(hull_points[:, 0], hull_points[:, 1],
                                 alpha=0.2, color='blue', label=self.lang.get_text('substrate'))
                except:
                    min_x, min_y = points.min(axis=0)
                    max_x, max_y = points.max(axis=0)
                    box = np.array([[min_x, min_y],
                                  [max_x, min_y],
                                  [max_x, max_y],
                                  [min_x, max_y],
                                  [min_x, min_y]])
                    self.axes.fill(box[:, 0], box[:, 1],
                                 alpha=0.2, color='blue', label=self.lang.get_text('substrate'))
            
            # 绘制分子投影（考虑周期性）
            if molecule is not None:
                positions = molecule.get_positions()
                points = positions[:, [x_idx, y_idx]]
                
                # 更新分子状态和颜色
                if self.projection_type == 'xy':
                    # 获取 XY 平面的状态
                    xy_status, xy_color, xy_alpha = is_molecule_pbc_valid(molecule, substrate, 'xy')
                    # 获取 XZ 平面的状态
                    xz_status, xz_color, xz_alpha = is_molecule_pbc_valid(molecule, substrate, 'xz')
                    # 获取 YZ 平面的状态
                    yz_status, yz_color, yz_alpha = is_molecule_pbc_valid(molecule, substrate, 'yz')
                    
                    # 如果任何一个平面显示红色（重叠），则使用红色
                    if xz_status == 'overlap' or yz_status == 'overlap':
                        ProjectionCanvas.molecule_status = 'overlap'
                        ProjectionCanvas.molecule_color = 'red'
                        ProjectionCanvas.molecule_mirror_alpha = 0.3
                    else:
                        ProjectionCanvas.molecule_status = xy_status
                        ProjectionCanvas.molecule_color = xy_color
                        ProjectionCanvas.molecule_mirror_alpha = xy_alpha
                elif self.projection_type == 'xz':
                    # 获取 XZ 平面的状态
                    xz_status, xz_color, xz_alpha = is_molecule_pbc_valid(molecule, substrate, 'xz')
                    ProjectionCanvas.molecule_status = xz_status
                    ProjectionCanvas.molecule_color = xz_color
                    ProjectionCanvas.molecule_mirror_alpha = xz_alpha
                else:  # yz
                    # 获取 YZ 平面的状态
                    yz_status, yz_color, yz_alpha = is_molecule_pbc_valid(molecule, substrate, 'yz')
                    ProjectionCanvas.molecule_status = yz_status
                    ProjectionCanvas.molecule_color = yz_color
                    ProjectionCanvas.molecule_mirror_alpha = yz_alpha
                
                # 使用类变量中的颜色
                color = ProjectionCanvas.molecule_color
                mirror_alpha = ProjectionCanvas.molecule_mirror_alpha
                
                # 绘制原始分子
                try:
                    hull = ConvexHull(points)
                    hull_points = points[hull.vertices]
                    hull_points = np.vstack([hull_points, hull_points[0]])
                    self.axes.fill(hull_points[:, 0], hull_points[:, 1],
                                 alpha=0.5, color=color, label=self.lang.get_text('molecule'))
                except:
                    min_x, min_y = points.min(axis=0)
                    max_x, max_y = points.max(axis=0)
                    box = np.array([[min_x, min_y],
                                  [max_x, min_y],
                                  [max_x, max_y],
                                  [min_x, max_y],
                                  [min_x, min_y]])
                    self.axes.fill(box[:, 0], box[:, 1],
                                 alpha=0.5, color=color, label=self.lang.get_text('molecule'))
                
                # 处理周期性边界条件
                if substrate is not None:
                    cell = substrate.get_cell()
                    if self.projection_type == 'xy':
                        cell_vectors = [cell[0, [x_idx, y_idx]], cell[1, [x_idx, y_idx]]]
                    elif self.projection_type == 'xz':
                        cell_vectors = [cell[0, [x_idx, y_idx]], cell[2, [x_idx, y_idx]]]
                    else:  # yz
                        cell_vectors = [cell[1, [x_idx, y_idx]], cell[2, [x_idx, y_idx]]]
                    
                    # 检查分子是否跨越边界
                    min_x, min_y = points.min(axis=0)
                    max_x, max_y = points.max(axis=0)
                    
                    # 定义需要复制的方向
                    directions = []
                    
                    # 检查第一个方向（X或Y）
                    if min_x < 0:
                        directions.append(np.array([1, 0]))
                    if max_x > cell_vectors[0][0]:
                        directions.append(np.array([-1, 0]))
                    
                    # 检查第二个方向（Y或Z）
                    if min_y < 0:
                        directions.append(np.array([0, 1]))
                    if max_y > cell_vectors[1][1]:
                        directions.append(np.array([0, -1]))
                    
                    # 检查对角方向
                    if (min_x < 0 and min_y < 0) or (max_x > cell_vectors[0][0] and max_y > cell_vectors[1][1]):
                        directions.append(np.array([1, 1]))
                    if (min_x < 0 and max_y > cell_vectors[1][1]) or (max_x > cell_vectors[0][0] and min_y < 0):
                        directions.append(np.array([1, -1]))
                    
                    # 绘制镜像分子
                    for direction in directions:
                        # 计算偏移量
                        offset = (direction[0] * cell_vectors[0] + 
                                direction[1] * cell_vectors[1])
                        
                        # 复制并平移分子
                        mirror_points = points + offset
                        
                        try:
                            hull = ConvexHull(mirror_points)
                            hull_points = mirror_points[hull.vertices]
                            hull_points = np.vstack([hull_points, hull_points[0]])
                            self.axes.fill(hull_points[:, 0], hull_points[:, 1],
                                         alpha=mirror_alpha, color=color,
                                         label=f'{self.lang.get_text("molecule")}({self.lang.get_text("mirror")})')
                        except:
                            min_x, min_y = mirror_points.min(axis=0)
                            max_x, max_y = mirror_points.max(axis=0)
                            box = np.array([[min_x, min_y],
                                          [max_x, min_y],
                                          [max_x, max_y],
                                          [min_x, max_y],
                                          [min_x, min_y]])
                            self.axes.fill(box[:, 0], box[:, 1],
                                         alpha=mirror_alpha, color=color,
                                         label=f'{self.lang.get_text("molecule")}({self.lang.get_text("mirror")})')
            
            # 标记锚点
            if hasattr(self, 'anchor_atom_index') and self.anchor_atom_index is not None and molecule is not None:
                anchor_pos = molecule.positions[self.anchor_atom_index]
                self.axes.scatter(anchor_pos[x_idx], anchor_pos[y_idx], 
                                color='yellow', marker='o', s=100, label=self.lang.get_text('anchor_point'))
                self.axes.annotate(f'({anchor_pos[x_idx]:.1f}, {anchor_pos[y_idx]:.1f})',
                                xy=(anchor_pos[x_idx], anchor_pos[y_idx]),
                                xytext=(5, 5), textcoords='offset points', fontsize=8)
            
            # 标记吸附位点
            if hasattr(self, 'site_position') and self.site_position is not None:
                self.axes.scatter(self.site_position[x_idx], self.site_position[y_idx],
                                color='red', marker='*', s=100, label=self.lang.get_text('adsorption_site'))
                self.axes.annotate(f'({self.site_position[x_idx]:.1f}, {self.site_position[y_idx]:.1f})',
                                xy=(self.site_position[x_idx], self.site_position[y_idx]),
                                xytext=(5, 5), textcoords='offset points', fontsize=8)
            
            # 只在顶部画布显示图例
            if self.with_legend:
                # 获取高度信息
                height_text = ""
                if substrate is not None:
                    substrate_height = substrate.get_cell()[2, 2]
                    height_text = f"{self.lang.get_text('lattice_height')}: {substrate_height:.1f} Å"
                if molecule is not None:
                    molecule_height = np.max(molecule.positions[:, 2])
                    if height_text:
                        height_text += f"\n{self.lang.get_text('molecule_highest_point')}: {molecule_height:.1f} Å"
                    else:
                        height_text = f"{self.lang.get_text('molecule_highest_point')}: {molecule_height:.1f} Å"
                
                # 添加高度信息到图例
                if height_text:
                    height_patch = mpatches.Patch(color='none', label=height_text)
                    handles, labels = self.axes.get_legend_handles_labels()
                    handles.append(height_patch)
                    labels.append(height_text)
                    self.axes.legend(handles, labels, loc='upper center',
                                bbox_to_anchor=(0.5, 1.2), ncol=2,
                                fontsize='small', framealpha=0.8)
            
            self.draw()
        except Exception as e:
            print(f"绘制投影时出错: {e}")
            import traceback
            traceback.print_exc()
            self.axes.clear()
            self.draw()
    
    def clear_plot(self):
        """Clear the projection plot"""
        self.axes.clear()
        self.draw()

class ExportOptionsDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        # 获取语言管理器实例
        self.lang = LanguageManager()
        self.setWindowTitle(self.lang.get_text('export_options'))
        self.initUI()
        
    def initUI(self):
        layout = QVBoxLayout()
        
        # 基础名称输入
        name_layout = QHBoxLayout()
        name_label = QLabel(self.lang.get_text('base_name'))
        self.name_input = QLineEdit()
        name_layout.addWidget(name_label)
        name_layout.addWidget(self.name_input)
        layout.addLayout(name_layout)
        
        # 导出内容选择
        content_group = QGroupBox(self.lang.get_text('export_content'))
        content_layout = QVBoxLayout()
        
        self.export_adsorption = QCheckBox(self.lang.get_text('export_adsorption_system'))
        self.export_adsorption.setChecked(True)
        content_layout.addWidget(self.export_adsorption)
        
        self.export_structures = QCheckBox(self.lang.get_text('export_substrate_and_molecule_structures'))
        self.export_structures.setChecked(True)
        content_layout.addWidget(self.export_structures)
        
        self.export_figure = QCheckBox(self.lang.get_text('export_picture'))
        self.export_figure.setChecked(True)
        content_layout.addWidget(self.export_figure)
        
        self.export_json = QCheckBox(self.lang.get_text('export_current_state'))
        self.export_json.setChecked(True)
        content_layout.addWidget(self.export_json)
        
        content_group.setLayout(content_layout)
        layout.addWidget(content_group)
        
        # 文件格式选择
        format_group = QGroupBox(self.lang.get_text('file_format'))
        format_layout = QVBoxLayout()
        
        self.format_combo = QComboBox()
        self.format_combo.addItems([
            "VASP Cartesian",
            "VASP Fractional",
            "XSF",
            "XYZ",
            "CIF"
        ])
        format_layout.addWidget(self.format_combo)
        
        format_group.setLayout(format_layout)
        layout.addWidget(format_group)
        
        # 按钮
        button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
            Qt.Horizontal, self)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
        
        self.setLayout(layout)
        
    def get_export_options(self):
        format_map = {
            "VASP Cartesian": "vasp_cartesian",
            "VASP Fractional": "vasp_fractional",
            "XSF": "xsf",
            "XYZ": "xyz",
            "CIF": "cif"
        }
        
        return {
            'base_name': self.name_input.text(),
            'export_adsorption': self.export_adsorption.isChecked(),
            'export_structures': self.export_structures.isChecked(),
            'export_figure': self.export_figure.isChecked(),
            'export_json': self.export_json.isChecked(),
            'format': format_map[self.format_combo.currentText()]
        }

class SiteSelectionDialog(QDialog):
    """Dialog for selecting adsorption sites on the substrate"""
    def __init__(self, substrate, site_type, element=None, parent=None):
        super().__init__(parent)
        self.substrate = substrate
        self.site_type = site_type
        self.element_filter = element
        self.parent_window = parent  # 保存父窗口引用，用于获取anchor信息
        
        # 获取语言管理器实例
        self.lang = LanguageManager()
        
        # 初始化变量
        self.sites = []
        self.site_elements = []  # 每个位点对应的元素列表
        self.selected_site = None
        self.filtered_sites = []  # 筛选后的位点
        self.filtered_site_elements = []  # 筛选后的位点元素
        self.filter_elements = []  # 用户选择的筛选元素
        
        self.initUI()
        self.plot_substrate_and_sites()
        
    def initUI(self):
        """Initialize the UI components"""
        # 使用site_type设置标题
        site_type_text = ""
        if self.site_type == "top":
            site_type_text = self.lang.get_text('site_types')[0]  # 顶位
        elif self.site_type == "bridge":
            site_type_text = self.lang.get_text('site_types')[1]  # 桥位
        elif self.site_type == "hollow":
            site_type_text = self.lang.get_text('site_types')[2]  # 空位
            
        self.setWindowTitle(f"{self.lang.get_text('select')}{site_type_text}{self.lang.get_text('adsorption_site')}")
        self.setMinimumSize(800, 600)
        
        # 创建主布局
        main_layout = QVBoxLayout()
        
        # 创建上部控制面板
        control_panel = QWidget()
        control_layout = QVBoxLayout(control_panel)
        
        # 创建筛选控件部分
        filter_group = QGroupBox(self.lang.get_text('site_filter'))
        filter_layout = QVBoxLayout()
        
        # 添加"是否筛选"复选框
        self.filter_checkbox = QCheckBox(self.lang.get_text('enable_element_filter'))
        self.filter_checkbox.stateChanged.connect(self.toggle_filter_controls)
        filter_layout.addWidget(self.filter_checkbox)
        
        # 添加表面层数选择
        surface_layer_layout = QHBoxLayout()
        surface_layer_layout.addWidget(QLabel(self.lang.get_text('surface_layer')))
        self.surface_layer_spin = QSpinBox()
        self.surface_layer_spin.setMinimum(1)
        self.surface_layer_spin.setMaximum(10)  # 最多允许10层
        self.surface_layer_spin.setValue(1)
        self.surface_layer_spin.valueChanged.connect(self.surface_layers_changed)
        surface_layer_layout.addWidget(self.surface_layer_spin)
        
        # 添加3D模式选择
        surface_layer_layout.addSpacing(20)  # 添加间距
        surface_layer_layout.addWidget(QLabel(self.lang.get_text('display_mode')))
        self.view_mode_combo = QComboBox()
        self.view_mode_combo.addItem(self.lang.get_text('top_layer'))
        self.view_mode_combo.addItem(self.lang.get_text('multi_layer'))
        self.view_mode_combo.currentIndexChanged.connect(self.view_mode_changed)
        surface_layer_layout.addWidget(self.view_mode_combo)
        
        surface_layer_layout.addStretch()
        filter_layout.addLayout(surface_layer_layout)
        
        # 创建元素选择区域（两列并排）
        elements_layout = QHBoxLayout()
        
        # 左侧：表面元素列表
        left_layout = QVBoxLayout()
        left_layout.addWidget(QLabel(self.lang.get_text('available_elements')))
        self.available_elements_list = QListWidget()
        self.available_elements_list.itemDoubleClicked.connect(self.add_filter_element)
        left_layout.addWidget(self.available_elements_list)
        elements_layout.addLayout(left_layout)
        
        # 中间：添加/移除按钮
        mid_layout = QVBoxLayout()
        mid_layout.addStretch()
        add_button = QPushButton(self.lang.get_text('add'))
        add_button.clicked.connect(self.add_filter_element)
        mid_layout.addWidget(add_button)
        remove_button = QPushButton(self.lang.get_text('remove'))
        remove_button.clicked.connect(self.remove_filter_element)
        mid_layout.addWidget(remove_button)
        clear_button = QPushButton(self.lang.get_text('clear'))
        clear_button.clicked.connect(self.clear_filter_elements)
        mid_layout.addWidget(clear_button)
        mid_layout.addStretch()
        elements_layout.addLayout(mid_layout)
        
        # 右侧：筛选元素列表
        right_layout = QVBoxLayout()
        right_layout.addWidget(QLabel(self.lang.get_text('filter_elements')))
        self.filter_elements_list = QListWidget()
        self.filter_elements_list.itemDoubleClicked.connect(self.remove_filter_element)
        right_layout.addWidget(self.filter_elements_list)
        elements_layout.addLayout(right_layout)
        
        # 添加元素选择区域到筛选布局
        filter_layout.addLayout(elements_layout)
        
        # 将筛选布局添加到筛选组
        filter_group.setLayout(filter_layout)
        control_layout.addWidget(filter_group)
        
        # 添加控制面板到主布局
        main_layout.addWidget(control_panel)
        
        # 创建图形显示区域
        self.canvas = MplCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.canvas.mpl_connect('button_press_event', self.on_click)
        
        main_layout.addWidget(self.toolbar)
        main_layout.addWidget(self.canvas)
        
        # 创建状态栏
        status_layout = QHBoxLayout()
        self.status_label = QLabel(self.lang.get_text('ready'))
        status_layout.addWidget(self.status_label)
        
        # 创建确定/取消按钮
        self.select_button = QPushButton(self.lang.get_text('select'))
        self.select_button.setEnabled(False)  # 初始时禁用
        self.select_button.clicked.connect(self.accept)
        status_layout.addWidget(self.select_button)
        
        cancel_button = QPushButton(self.lang.get_text('cancel'))
        cancel_button.clicked.connect(self.reject)
        status_layout.addWidget(cancel_button)
        
        main_layout.addLayout(status_layout)
        
        # 设置布局
        self.setLayout(main_layout)
        
        # 设置筛选控件初始状态
        self.toggle_filter_controls(Qt.Unchecked)
        
        # 填充可用元素列表
        self.populate_available_elements()
    
    def toggle_filter_controls(self, state):
        """启用或禁用筛选控件"""
        enabled = bool(state)
        self.available_elements_list.setEnabled(enabled)
        self.filter_elements_list.setEnabled(enabled)
        
        # 如果禁用筛选，则清空筛选条件并重绘
        if not enabled:
            self.clear_filter_elements()
        
    def populate_available_elements(self):
        """填充可用元素列表，基于当前选择的表面层数"""
        # 获取当前选定的表面层数
        surface_layers = self.surface_layer_spin.value()
        
        # 获取表面原子
        surface_indices = find_surface_atoms(self.substrate, surface_layers=surface_layers)
        
        # 获取表面原子的元素
        symbols = [self.substrate.get_chemical_symbols()[i] for i in surface_indices]
        unique_elements = sorted(set(symbols))
        
        # 清空现有列表
        self.available_elements_list.clear()
        
        # 添加元素到列表
        for element in unique_elements:
            self.available_elements_list.addItem(element)
    
    def add_filter_element(self):
        """添加选定的元素到筛选列表"""
        current_item = self.available_elements_list.currentItem()
        if current_item is not None:
            element = current_item.text()
            # 检查是否已在筛选列表中
            existing_items = [self.filter_elements_list.item(i).text() 
                             for i in range(self.filter_elements_list.count())]
            if element not in existing_items:
                self.filter_elements_list.addItem(element)
                self.filter_elements.append(element)
                self.apply_filter()
    
    def remove_filter_element(self):
        """从筛选列表中移除选定的元素"""
        current_item = self.filter_elements_list.currentItem()
        if current_item is not None:
            element = current_item.text()
            # 从列表和筛选条件中删除
            row = self.filter_elements_list.row(current_item)
            self.filter_elements_list.takeItem(row)
            if element in self.filter_elements:
                self.filter_elements.remove(element)
            self.apply_filter()
    
    def clear_filter_elements(self):
        """清空筛选元素列表"""
        self.filter_elements_list.clear()
        self.filter_elements = []
        self.apply_filter()
    
    def apply_filter(self):
        """应用筛选条件，更新显示的吸附位点"""
        if not self.filter_checkbox.isChecked() or not self.filter_elements:
            # 如果筛选未启用或没有筛选元素，显示所有位点
            self.filtered_sites = self.sites
            self.filtered_site_elements = self.site_elements
        else:
            # 应用筛选，只保留包含所有筛选元素的位点
            self.filtered_sites = []
            self.filtered_site_elements = []
            
            for i, site_element_list in enumerate(self.site_elements):
                # 检查此位点是否包含所有筛选元素
                if self.site_type == "top":
                    # 对于top位点，元素必须匹配
                    if site_element_list and site_element_list[0] in self.filter_elements:
                        self.filtered_sites.append(self.sites[i])
                        self.filtered_site_elements.append(site_element_list)
                else:
                    # 对于bridge和hollow位点，检查是否包含所有筛选元素
                    contains_all = all(elem in site_element_list for elem in self.filter_elements)
                    if contains_all:
                        self.filtered_sites.append(self.sites[i])
                        self.filtered_site_elements.append(site_element_list)
        
        # 更新位点显示
        self.plot_substrate_and_sites()
        
        # 更新状态标签
        self.status_label.setText(f"符合条件的位点: {len(self.filtered_sites)}/{len(self.sites)}")

    def plot_substrate_and_sites(self):
        """Plot substrate and available adsorption sites"""
        # 保存当前视角，以便重绘后恢复
        if hasattr(self.canvas.axes, 'elev') and hasattr(self.canvas.axes, 'azim'):
            current_elev = self.canvas.axes.elev
            current_azim = self.canvas.axes.azim
        else:
            current_elev = 30
            current_azim = 45
        
        # 清除当前图形
        self.canvas.axes.clear()
        
        # 获取基底原子位置和元素
        positions = self.substrate.positions
        symbols = self.substrate.get_chemical_symbols()
        
        # 使用全局元素颜色表
        self.default_color = DEFAULT_COLOR
        
        # 设置3D图形的一些属性
        self.canvas.axes.set_xlabel('X (Å)')
        self.canvas.axes.set_ylabel('Y (Å)')
        self.canvas.axes.set_zlabel('Z (Å)')
        
        # 绘制基底原子 (根据选定的表面层数)
        try:
            # 获取当前选定的表面层数
            surface_layers = self.surface_layer_spin.value()
            
            # 获取当前的显示模式
            is_3d_mode = self.view_mode_combo.currentIndex() == 1  # 1表示多层3D模式
            
            # 根据显示模式决定如何获取表面原子
            if is_3d_mode:
                # 3D模式：直接按Z坐标排序显示更多原子
                all_atoms_indices = list(range(len(self.substrate)))
                z_coords = positions[:, 2]
                
                # 计算Z坐标统计信息，用于确定显示哪些原子
                min_z = np.min(z_coords)
                max_z = np.max(z_coords)
                z_range = max_z - min_z
                
                # 通过"阈值百分比"选择表面原子（显示上面60%高度范围的原子）
                # 这样可以显示更多的原子层而不是简单地选择前N个
                z_threshold = min_z + z_range * 0.4  # 显示上面60%高度范围的原子
                
                # 选择高于阈值的所有原子
                surface_indices = np.where(z_coords > z_threshold)[0]
                
                # 如果选出的原子太少，至少选择20个
                if len(surface_indices) < 20:
                    # 按Z坐标降序排列所有原子
                    sorted_indices = [i for _, i in sorted(zip(z_coords, all_atoms_indices), key=lambda x: x[0], reverse=True)]
                    surface_indices = sorted_indices[:min(20, len(sorted_indices))]
            else:
                # 传统模式：仅显示表面层
                surface_indices = find_surface_atoms(self.substrate, surface_layers=surface_layers)
            
            # 计算Z坐标范围，用于颜色和大小渐变
            if len(surface_indices) > 0:
                z_values = positions[surface_indices, 2]
                min_z = np.min(z_values)
                max_z = np.max(z_values)
                z_range = max_z - min_z if max_z > min_z else 1.0
            else:
                min_z = max_z = 0
                z_range = 1.0
            
            # 在三维空间中绘制表面原子
            for idx in surface_indices:
                pos = positions[idx]
                symbol = symbols[idx]
                color = ELEMENT_COLORS[symbol]
                
                # 根据深度调整大小和透明度
                depth_factor = (pos[2] - min_z) / z_range if z_range > 0 else 1.0
                
                # 越深的原子越小且越透明
                size = 50 + 50 * depth_factor
                alpha = 0.5 + 0.5 * depth_factor
                
                self.canvas.axes.scatter(pos[0], pos[1], pos[2], color=color, s=size, alpha=alpha)
            
            # 添加图例 - 只添加显示的元素，而不是所有元素
            # 收集当前表面层显示的元素
            displayed_elements = set()
            for idx in surface_indices:
                displayed_elements.add(symbols[idx])
            
            # 为每种显示的元素添加图例项
            legend_elements = []
            for element in displayed_elements:
                color = ELEMENT_COLORS.get(element, DEFAULT_COLOR)
                legend_elements.append(
                    plt.Line2D([0], [0], marker='o', color='w', 
                               label=element, markerfacecolor=color, markersize=10,
                               markeredgecolor='black', markeredgewidth=1.0)
                )
            
            # 将图例放在图表右上角，但稍微调整位置避免遮挡
            if legend_elements:
                self.canvas.axes.legend(handles=legend_elements, loc='upper right', 
                                      bbox_to_anchor=(0.95, 0.95), fontsize=10, 
                                      frameon=True, framealpha=0.8, facecolor='white')
            
            # 查找并绘制吸附位点
            self.sites = []
            self.site_elements = []
            
            if self.site_type == "top":
                if is_3d_mode:
                    # 在3D模式下，使用更多的原子作为顶位
                    all_atoms_indices = list(range(len(self.substrate)))
                    z_coords = positions[:, 2]
                    
                    # 计算Z坐标统计信息
                    min_z = np.min(z_coords)
                    max_z = np.max(z_coords)
                    z_range = max_z - min_z
                    
                    # 使用与显示原子相同的阈值来选择位点
                    z_threshold = min_z + z_range * 0.4
                    top_indices = np.where(z_coords > z_threshold)[0]
                    
                    # 如果选出的位点太少，至少选择25个
                    if len(top_indices) < 25:
                        sorted_indices = [i for _, i in sorted(zip(z_coords, all_atoms_indices), key=lambda x: x[0], reverse=True)]
                        top_indices = sorted_indices[:min(25, len(sorted_indices))]
                    
                    # 生成位点
                    for idx in top_indices:
                        pos = positions[idx]
                        self.sites.append([pos[0], pos[1], pos[2]])
                        self.site_elements.append([symbols[idx]])
                else:
                    # 传统模式：Top sites are directly on surface atoms
                    for idx in surface_indices:
                        pos = positions[idx]
                        self.sites.append([pos[0], pos[1], pos[2]])
                        self.site_elements.append([symbols[idx]])
                    
            elif self.site_type == "bridge":
                # Find bridge sites
                if is_3d_mode:
                    # 在3D模式下，使用当前用户设置的层数，而不是硬编码的5
                    # 修改：使用当前设置的表面层数
                    bridge_sites, bridge_elements = find_bridge_sites(self.substrate, max_distance=4.0, surface_layers=surface_layers)
                else:
                    bridge_sites, bridge_elements = find_bridge_sites(self.substrate, surface_layers=surface_layers)
                
                for site, elements in zip(bridge_sites, bridge_elements):
                    self.sites.append(site)
                    self.site_elements.append(elements)
                    
            elif self.site_type == "hollow":
                # Find hollow sites
                if is_3d_mode:
                    # 在3D模式下，使用当前用户设置的层数，而不是硬编码的5
                    # 修改：使用当前设置的表面层数
                    hollow_sites, hollow_elements = find_hollow_sites(self.substrate, max_distance=5.5, surface_layers=surface_layers)
                else:
                    hollow_sites, hollow_elements = find_hollow_sites(self.substrate, surface_layers=surface_layers)
                
                for site, elements in zip(hollow_sites, hollow_elements):
                    self.sites.append(site)
                    self.site_elements.append(elements)
            
            # 应用筛选条件
            if not hasattr(self, 'filtered_sites') or self.filtered_sites is None or len(self.filtered_sites) == 0:
                self.apply_filter()
            
            # 绘制筛选后的位点
            if hasattr(self, 'filtered_sites') and self.filtered_sites is not None and len(self.filtered_sites) > 0:
                # 原代码：一次性绘制所有位点
                # filtered_sites_array = np.array(self.filtered_sites)
                # self.canvas.axes.scatter(
                #     filtered_sites_array[:, 0], 
                #     filtered_sites_array[:, 1], 
                #     filtered_sites_array[:, 2], 
                #     color='green', marker='x', s=50
                # )
                
                # 新代码：为每个位点单独绘制，以便根据深度调整外观
                for i, site in enumerate(self.filtered_sites):
                    x, y, z = site
                    
                    # 计算深度因子
                    depth_factor = (z - min_z) / z_range if z_range > 0 else 1.0
                    # 限制depth_factor在0-1范围内
                    depth_factor = max(0.0, min(1.0, depth_factor))
                    
                    # 调整位点大小和透明度
                    size = 30 + 40 * depth_factor  # 大小从30到70
                    marker_alpha = 0.4 + 0.6 * depth_factor  # 透明度从0.4到1.0
                    
                    # 根据位点类型设置颜色
                    if self.site_type == "top":
                        site_color = 'red'
                    elif self.site_type == "bridge":
                        site_color = 'blue'
                    else:  # hollow或其他位点
                        site_color = 'green'
                    
                    self.canvas.axes.scatter(
                        x, y, z, 
                        color=site_color, marker='x', s=size, alpha=marker_alpha
                    )
                
                # 在位点旁添加索引编号，便于用户识别
                for i, site in enumerate(self.filtered_sites):
                    x, y, z = site[0], site[1], site[2]
                    # 计算Z轴深度，使得深层位点的文字更小更透明
                    depth_factor = (z - min_z) / z_range if z_range > 0 else 1.0
                    # 限制depth_factor在0-1范围内
                    depth_factor = max(0.0, min(1.0, depth_factor))
                    text_alpha = 0.5 + 0.5 * depth_factor  # 透明度从0.5到1.0
                    text_size = 6 + 2 * depth_factor  # 字体大小从6到8
                    
                    # 设置文本颜色与位点颜色相同
                    if self.site_type == "top":
                        text_color = 'red'
                    elif self.site_type == "bridge":
                        text_color = 'blue'
                    else:  # hollow或其他位点
                        text_color = 'green'
                        
                    self.canvas.axes.text(x, y, z, f"{i+1}", 
                                       fontsize=text_size, ha='right', va='bottom',
                                       color=text_color, weight='bold', alpha=text_alpha)
            
            # 标记anchor原子位置（如果在父窗口有选择anchor原子且分子已经加载）
            if self.parent_window and hasattr(self.parent_window, 'molecule') and self.parent_window.molecule is not None and hasattr(self.parent_window, 'anchor'):
                # 获取锚点原子索引和位置
                anchor_idx = self.parent_window.anchor.value()
                if anchor_idx < len(self.parent_window.molecule):
                    # 如果有选定位点，以选定位点为基准显示anchor位置
                    if hasattr(self.parent_window, 'selected_site') and self.parent_window.selected_site is not None:
                        # 获取父窗口的当前anchor原子位置和选定位点
                        site_pos = self.parent_window.selected_site
                        
                        # 定义锚点标记颜色
                        ANCHOR_COLOR = '#FF0000'  # 鲜红色
                        
                        # 在吸附位点上方标记anchor原子位置
                        self.canvas.axes.scatter(
                            site_pos[0], site_pos[1], site_pos[2], 
                            color=ANCHOR_COLOR, marker='o', s=100, alpha=0.7,
                            edgecolors='black', linewidths=1.0,
                            label=self.lang.get_text("current_anchor_position")
                        )
            
            # 设置图形标题
            site_type_text = ""
            if self.site_type == "top":
                site_type_text = self.lang.get_text('site_types')[0]  # 顶位
            elif self.site_type == "bridge":
                site_type_text = self.lang.get_text('site_types')[1]  # 桥位
            elif self.site_type == "hollow":
                site_type_text = self.lang.get_text('site_types')[2]  # 空位
            
            self.canvas.axes.set_title(f"{self.site_type.capitalize()} {site_type_text}")
            
            # 设置适合的显示范围
            positions_array = np.array(positions)
            filtered_sites_array = np.array(self.filtered_sites) if self.filtered_sites else np.array([])
            
            if len(filtered_sites_array) > 0:
                # 计算包含基底和位点的边界框
                all_points = np.vstack((positions_array, filtered_sites_array))
                
                x_min, y_min, z_min = np.min(all_points, axis=0)
                x_max, y_max, z_max = np.max(all_points, axis=0)
                
                # 稍微放大视图
                padding = max(x_max - x_min, y_max - y_min, z_max - z_min) * 0.1
                
                self.canvas.axes.set_xlim(x_min - padding, x_max + padding)
                self.canvas.axes.set_ylim(y_min - padding, y_max + padding)
                self.canvas.axes.set_zlim(z_min - padding, z_max + padding)
            else:
                # 只显示基底
                x_min, y_min, z_min = np.min(positions_array, axis=0)
                x_max, y_max, z_max = np.max(positions_array, axis=0)
                
                padding = max(x_max - x_min, y_max - y_min, z_max - z_min) * 0.1
                
                self.canvas.axes.set_xlim(x_min - padding, x_max + padding)
                self.canvas.axes.set_ylim(y_min - padding, y_max + padding)
                self.canvas.axes.set_zlim(z_min - padding, z_max + padding)
            
            # 添加图例，但去除重复项
            handles, labels = self.canvas.axes.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            self.canvas.axes.legend(by_label.values(), by_label.keys(), 
                                  loc='upper right', fontsize='small')
            
            # 恢复之前的视角
            self.canvas.axes.view_init(elev=current_elev, azim=current_azim)
            
            # 自动调整布局
            self.canvas.figure.tight_layout()
            
        except Exception as e:
            import traceback
            traceback.print_exc()
            # 显示异常信息
            self.canvas.axes.text(0.5, 0.5, 0.5, f"Error plotting sites: {str(e)}", 
                                ha='center', va='center', transform=self.canvas.axes.transAxes)
        
        # 刷新画布
        self.canvas.draw()
    
    def on_click(self, event):
        """处理鼠标点击事件，选择最近的位点"""
        if event.inaxes != self.canvas.axes:
            return
            
        # 保存当前视角角度
        current_elev = self.canvas.axes.elev
        current_azim = self.canvas.axes.azim
        
        # 获取点击位置的2D坐标（图形坐标系）
        click_x, click_y = event.xdata, event.ydata
        
        # 添加调试输出
        print(f"{self.lang.get_text('click_position')}: x={click_x:.2f}, y={click_y:.2f}")
        print(f"{self.lang.get_text('available_sites_number')}: {len(self.filtered_sites)}")
        print(f"{self.lang.get_text('current_view')}: elev={current_elev}, azim={current_azim}")
        
        # 如果没有可选位点，直接返回
        if len(self.filtered_sites) == 0:
            return
        
        # 获取3D坐标轴
        ax = self.canvas.axes
        
        # 获取当前视图的投影变换
        proj = ax.get_proj()
        
        # 计算所有点在当前视图中的投影坐标
        min_dist = float('inf')
        nearest_idx = None
        
        # 将点击位置转换为绝对像素位置
        click_pixel = event.x, event.y
        
        for i, site in enumerate(self.filtered_sites):
            # 将3D点投影到2D平面
            x, y, z = site
            x2d, y2d, _ = proj3d.proj_transform(x, y, z, proj)
            
            # 转换为图像坐标系下的坐标
            xy_pixel = ax.transData.transform((x2d, y2d))
            
            # 计算在像素坐标系中到点击位置的距离
            dist = np.sqrt((xy_pixel[0] - click_pixel[0])**2 + 
                          (xy_pixel[1] - click_pixel[1])**2)
            
            # 调试输出
            print(f"{self.lang.get_text('site')} {i}: {self.lang.get_text('position')}={site}, {self.lang.get_text('distance')}={dist:.2f}像素")
            
            # 更新最近点
            if dist < min_dist:
                min_dist = dist
                nearest_idx = i
        
        if nearest_idx is not None:
            nearest_site = self.filtered_sites[nearest_idx]
            
            print(f"{self.lang.get_text('selected_site_index')}: {nearest_idx}, {self.lang.get_text('position')}: {nearest_site}, {self.lang.get_text('distance')}: {min_dist:.2f}像素")
            
            # 保存选中位点的信息（在绘制之前）
            self.selected_site = nearest_site
            self.select_button.setEnabled(True)
            
            # 重新绘制基底和所有位点
            self.plot_substrate_and_sites()
            
            # 恢复视角设置
            self.canvas.axes.view_init(elev=current_elev, azim=current_azim)
            
            # 突出显示选中的位点
            self.canvas.axes.scatter([nearest_site[0]], [nearest_site[1]], [nearest_site[2]],
                                     color='red', s=200, marker='o', zorder=10,
                                     edgecolors='black')
            
            # 在位点上方添加清晰的标记
            self.canvas.axes.text(nearest_site[0], nearest_site[1], nearest_site[2] + 0.5, 
                                 f"{self.lang.get_text('selected_site')} {nearest_idx+1}", 
                                 color='red', fontweight='bold', ha='center', va='bottom',
                                 zorder=11)
            
            # 刷新画布
            self.canvas.draw()
            
            # 更新状态标签显示选定位点的元素组成
            site_elements = self.filtered_site_elements[nearest_idx]
            element_str = ", ".join(site_elements)
            self.status_label.setText(f"{self.lang.get_text('selected_site')} {nearest_idx+1}/{len(self.filtered_sites)}: "
                                     f"{nearest_site} ({self.lang.get_text('element')}: {element_str}, {self.lang.get_text('distance')}: {min_dist:.2f}像素)")
    
    def accept(self):
        """Accept currently selected site"""
        super(SiteSelectionDialog, self).accept()

    def surface_layers_changed(self, value):
        """当表面层数改变时更新界面"""
        # 更新可用元素列表
        self.populate_available_elements()
        # 重新绘制位点
        self.plot_substrate_and_sites()
        # 清空筛选条件，因为表面原子可能变化
        self.clear_filter_elements()
        
    def view_mode_changed(self, index):
        """当显示模式改变时更新界面"""
        # 记录一下当前的层数值
        current_layers = self.surface_layer_spin.value()
        
        # 清除筛选条件和当前选择的位点
        if hasattr(self, 'filtered_sites'):
            self.filtered_sites = []
        self.clear_filter_elements()
        
        # 如果是从3D模式切换到传统模式，确保层数设置合理
        if index == 0 and current_layers > 3:  # 传统模式且层数大于3
            # 在传统模式下，过多的层数可能导致显示混乱，默认设为1
            self.surface_layer_spin.setValue(1)
            # 设置值会触发surface_layers_changed，它会调用plot_substrate_and_sites
        else:
            # 直接重绘位点
            self.plot_substrate_and_sites()

class AtomSelectionDialog(QDialog):
    """选择原子对话框"""
    def __init__(self, atoms, parent=None):
        super().__init__(parent)
        self.atoms = atoms
        self.selected_atom = None
        self.center_perpendicular = None
        self.c_perpendicular = None
        self.no_perpendicular = None
        self.adjustment_group = None
        # 获取语言管理器实例
        self.lang = LanguageManager()
        self.initUI()
    
    def initUI(self):
        """Initialize the dialog UI"""
        self.setWindowTitle(self.lang.get_text('select_anchor_atom'))
        self.resize(900, 600)
        
        # Main layout
        layout = QVBoxLayout()
        
        # Create top panel with controls
        top_panel = QWidget()
        top_layout = QHBoxLayout(top_panel)
        
        # Element selection
        element_label = QLabel(self.lang.get_text('element') + ":")
        self.element_combo = QComboBox()
        
        # Get unique elements
        unique_elements = sorted(set(self.atoms.get_chemical_symbols()))
        self.element_combo.addItems(unique_elements)
        
        # Connect signal to update atom indices when element changes
        self.element_combo.currentTextChanged.connect(self._update_atom_combo)
        
        # Atom selection
        atom_label = QLabel(self.lang.get_text('atom_index'))
        self.atom_combo = QComboBox()
        
        # Connect signal to update selection when atom changes
        self.atom_combo.currentTextChanged.connect(self._update_selection)
        
        # Add elements to top panel
        top_layout.addWidget(element_label)
        top_layout.addWidget(self.element_combo)
        top_layout.addWidget(atom_label)
        top_layout.addWidget(self.atom_combo)
        
        # Create info display area
        self.atom_info = QLabel(self.lang.get_text('select_an_atom_to_view_information'))
        self.atom_info.setMinimumHeight(50)
        self.atom_info.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        self.atom_info.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.atom_info.setWordWrap(True)
        
        # Create 3D visualization widget
        self.canvas = MplCanvas(self, width=7, height=5, dpi=100)
        self.canvas_widget = QWidget()
        canvas_layout = QVBoxLayout(self.canvas_widget)
        canvas_layout.addWidget(self.canvas)
        
        # Create checkbox for adjustment options
        adjustment_groupbox = QGroupBox(self.lang.get_text('direction_adjustment_options'))
        adjustment_layout = QVBoxLayout()
        
        # 创建选项按钮组，确保选项互斥
        self.adjustment_group = QButtonGroup()
        
        # 第一个选项：与分子质心连线垂直于表面
        self.center_perpendicular = QRadioButton(self.lang.get_text('anchor_perpendicular_to_molecule_center'))
        self.center_perpendicular.setChecked(True)  # 默认选中
        self.adjustment_group.addButton(self.center_perpendicular)
        adjustment_layout.addWidget(self.center_perpendicular)
        
        # 第二个选项：与C原子连线垂直于表面
        self.c_perpendicular = QRadioButton(self.lang.get_text('anchor_perpendicular_to_C_atom'))
        self.adjustment_group.addButton(self.c_perpendicular)
        adjustment_layout.addWidget(self.c_perpendicular)

        # 第三个选项：不做改变
        self.no_perpendicular = QRadioButton(self.lang.get_text('keep_current_molecule_state'))
        self.adjustment_group.addButton(self.no_perpendicular)
        adjustment_layout.addWidget(self.no_perpendicular)
        
        # 添加是否需要调整分子位置的选项
        self.adjust_checkbox = QCheckBox(self.lang.get_text('adjust_molecule_position'))
        self.adjust_checkbox.setChecked(True)  # 默认选中
        adjustment_layout.addWidget(self.adjust_checkbox)
        
        adjustment_groupbox.setLayout(adjustment_layout)
        
        # Create button layout
        button_layout = QHBoxLayout()
        self.ok_button = QPushButton(self.lang.get_text('ok'))
        self.cancel_button = QPushButton(self.lang.get_text('cancel'))
        
        # Connect buttons
        self.ok_button.clicked.connect(self.accept)
        self.cancel_button.clicked.connect(self.reject)
        
        button_layout.addStretch()
        button_layout.addWidget(self.ok_button)
        button_layout.addWidget(self.cancel_button)
        
        # Combine all components in the main layout
        layout.addWidget(top_panel)
        layout.addWidget(self.atom_info)
        layout.addWidget(self.canvas_widget)
        layout.addWidget(adjustment_groupbox)  # 添加调整选项组
        layout.addLayout(button_layout)
        
        self.setLayout(layout)
        
        # Initialize with first element
        if unique_elements:
            self._update_atom_combo(unique_elements[0])
        
        # Draw the molecule
        self.plot_molecule()
    
    def _update_atom_combo(self, element):
        """根据选择的元素更新原子索引下拉菜单"""
        self.atom_combo.clear()
        
        if not element:
            return
            
        symbols = self.atoms.get_chemical_symbols()
        
        # 查找所有匹配元素的原子索引
        for i, symbol in enumerate(symbols):
            if symbol == element:
                self.atom_combo.addItem(f"{i}: {symbol}")
                
        # 如果有原子，选择第一个
        if self.atom_combo.count() > 0:
            self.atom_combo.setCurrentIndex(0)
            self._update_selection(self.atom_combo.currentText())
    
    def _update_selection(self, atom_text):
        """根据选择的原子文本更新所选原子"""
        if not atom_text:
            return
            
        # 从文本中提取索引 (格式: "索引: 元素")
        try:
            atom_idx = int(atom_text.split(":")[0].strip())
            symbol = self.atoms.get_chemical_symbols()[atom_idx]
            
            # 更新UI
            self.atom_info.setText(f"{self.lang.get_text('selected_atom')}: {atom_idx} ({symbol})")
            self.selected_atom = atom_idx
            
            # 高亮显示选中的原子
            self.highlight_selected_atom(atom_idx)
        except (ValueError, IndexError):
            self.atom_info.setText(self.lang.get_text('error'))
            self.selected_atom = None
    
    def plot_molecule(self):
        """绘制分子结构并标注原子索引"""
        self.canvas.figure.clear()
        ax = self.canvas.figure.add_subplot(111, projection='3d')
        
        # 获取位置和元素符号
        positions = self.atoms.get_positions()
        symbols = self.atoms.get_chemical_symbols()
        
        # 使用全局定义的元素颜色表
        self.default_color = DEFAULT_COLOR  # 保存默认颜色作为实例变量以便其他方法使用
        
        # 绘制分子的键
        self.draw_molecular_bonds(self.atoms, ax)
        
        # 绘制原子
        for i, (pos, symbol) in enumerate(zip(positions, symbols)):
            x, y, z = pos
            color = ELEMENT_COLORS.get(symbol, self.default_color)
            # 添加edgecolors='black'参数为原子添加黑色边框
            ax.scatter(x, y, z, color=color, s=100, edgecolors='black', linewidths=1.0)
            # 添加原子索引标签，稍微偏移位置以避免被遮挡
            ax.text(x+0.1, y+0.1, z+0.1, f"{i}:{symbol}", size=8, zorder=5, color='black')
        
        # 添加图例
        legend_elements = []
        for element in set(symbols):
            color = ELEMENT_COLORS.get(element, self.default_color)
            legend_elements.append(
                plt.Line2D([0], [0], marker='o', color='w', 
                          label=element, markerfacecolor=color, markersize=10,
                          markeredgecolor='black', markeredgewidth=1.0)
            )
        
        # 优化图例位置，避免遮挡
        if legend_elements:
            ax.legend(handles=legend_elements, loc='upper right', 
                     bbox_to_anchor=(0.95, 0.95), fontsize=10, 
                     frameon=True, framealpha=0.8, facecolor='white')
        
        # 设置标签
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        
        # 设置等比例视图
        max_range = np.array([
            positions[:, 0].max() - positions[:, 0].min(),
            positions[:, 1].max() - positions[:, 1].min(),
            positions[:, 2].max() - positions[:, 2].min()
        ]).max() / 2.0
        
        mid_x = (positions[:, 0].max() + positions[:, 0].min()) * 0.5
        mid_y = (positions[:, 1].max() + positions[:, 1].min()) * 0.5
        mid_z = (positions[:, 2].max() + positions[:, 2].min()) * 0.5
        
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        self.canvas.draw()
    
    def highlight_selected_atom(self, atom_idx):
        """重绘图表并高亮显示选中的原子"""
        self.canvas.figure.clear()
        ax = self.canvas.figure.add_subplot(111, projection='3d')
        
        # 获取位置和元素符号
        positions = self.atoms.get_positions()
        symbols = self.atoms.get_chemical_symbols()
        
        # 使用全局定义的元素颜色表
        if not hasattr(self, 'default_color'):
            self.default_color = DEFAULT_COLOR
        
        # 定义选中原子的高亮颜色
        HIGHLIGHT_COLOR = '#FF0000'  # 鲜红色
        
        # 绘制分子的键
        self.draw_molecular_bonds(self.atoms, ax)
        
        # 绘制原子
        for i, (pos, symbol) in enumerate(zip(positions, symbols)):
            x, y, z = pos
            
            # 设置原子大小和颜色
            if i == atom_idx:
                # 高亮显示选中的原子
                color = HIGHLIGHT_COLOR
                size = 200
            else:
                color = ELEMENT_COLORS.get(symbol, self.default_color)
                size = 100
                
            # 添加edgecolors='black'参数为原子添加黑色边框
            ax.scatter(x, y, z, color=color, s=size, edgecolors='black', linewidths=1.0)
            
            # 添加原子索引标签，稍微偏移位置以避免被遮挡
            ax.text(x+0.1, y+0.1, z+0.1, f"{i}:{symbol}", size=8, zorder=5, color='black')
        
        # 添加图例
        legend_elements = []
        for element in set(symbols):
            color = ELEMENT_COLORS.get(element, self.default_color)
            legend_elements.append(
                plt.Line2D([0], [0], marker='o', color='w', 
                          label=element, markerfacecolor=color, markersize=10,
                          markeredgecolor='black', markeredgewidth=1.0)
            )
        
        # 添加选中原子的图例
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                                label=self.lang.get_text('selected_atom'), markerfacecolor=HIGHLIGHT_COLOR, markersize=10))
        
        # 优化图例位置，避免遮挡
        if legend_elements:
            ax.legend(handles=legend_elements, loc='upper right', 
                     bbox_to_anchor=(0.95, 0.95), fontsize=10, 
                     frameon=True, framealpha=0.8, facecolor='white')
        
        # 设置标签
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        
        # 设置等比例视图
        max_range = np.array([
            positions[:, 0].max() - positions[:, 0].min(),
            positions[:, 1].max() - positions[:, 1].min(),
            positions[:, 2].max() - positions[:, 2].min()
        ]).max() / 2.0
        
        mid_x = (positions[:, 0].max() + positions[:, 0].min()) * 0.5
        mid_y = (positions[:, 1].max() + positions[:, 1].min()) * 0.5
        mid_z = (positions[:, 2].max() + positions[:, 2].min()) * 0.5
        
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        self.canvas.draw()
    
    def draw_molecular_bonds(self, atoms, ax):
        """绘制分子的键连接"""
        # 获取原子位置和元素符号
        positions = atoms.positions
        symbols = atoms.get_chemical_symbols()
        
        # 定义共价半径（单位：埃）
        covalent_radii = {
            'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 
            'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 
            'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
            'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39,
            'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22,
            'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 'Br': 1.20, 'Kr': 1.16,
            'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
            'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44,
            'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
            'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01,
            'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92,
            'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75,
            'Ta': 1.70, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36,
            'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40,
            'At': 1.50, 'Rn': 1.50
        }
        
        # 定义键的绘制参数
        bond_color = '#777777'  # 灰色
        bond_linewidth = 1.5
        bond_tolerance = 1.3  # 容差因子，决定两个原子之间是否连接（两个原子的共价半径之和 * 容差因子）
        
        # 计算所有原子对之间的距离并检查是否需要绘制键
        n_atoms = len(atoms)
        for i in range(n_atoms):
            pos_i = positions[i]
            symbol_i = symbols[i]
            radius_i = covalent_radii.get(symbol_i, 0.7)  # 默认值为0.7埃
            
            for j in range(i+1, n_atoms):
                pos_j = positions[j]
                symbol_j = symbols[j]
                radius_j = covalent_radii.get(symbol_j, 0.7)  # 默认值为0.7埃
                
                # 计算两原子间距离
                distance = np.linalg.norm(pos_i - pos_j)
                
                # 检查是否应该绘制键
                if distance <= (radius_i + radius_j) * bond_tolerance:
                    # 绘制连接线
                    ax.plot([pos_i[0], pos_j[0]], [pos_i[1], pos_j[1]], [pos_i[2], pos_j[2]], 
                           color=bond_color, linewidth=bond_linewidth, alpha=0.7)
        
        return ax
    
    def get_adjustment_option(self):
        """返回是否需要调整分子位置的选项和调整类型"""
        adjust = self.adjust_checkbox.isChecked()
        
        # 确定调整类型
        if self.center_perpendicular.isChecked():
            adjust_type = "center"  # 与分子质心连线垂直
        elif self.c_perpendicular.isChecked():
            adjust_type = "carbon"  # 与相连C原子连线垂直
        else:  # self.no_perpendicular.isChecked()
            adjust_type = "none"    # 保持当前分子状态
            # 如果选择了保持当前状态，强制不调整
            adjust = False
            
        return adjust, adjust_type

class AdsorptionGUI(QMainWindow):
    """Molecular Adsorption Simulation Interface"""
    def __init__(self):
        """Initialize the AdsorptionGUI application"""
        super().__init__()
        
        # 初始化语言设置
        self.current_language = 'en'  # 默认使用英语
        # 获取语言管理器实例
        self.lang = LanguageManager()
        
        # Set up the UI
        self.initUI()
        
        # 创建菜单栏
        self.create_menu_bar()
        
        # Initialize variables
        self.molecule = None
        self.substrate = None
        self.adsorption_system = None
        self.original_positions = None
        self.selected_site = None
        
        # 添加变量保存原始分子和初始分子结构
        self.original_molecule = None  # 保存原始未经任何处理的分子
        self.initial_molecule = None   # 保存应用anchor_type后的初始分子结构
        
        # 添加跟踪垂直调整的属性
        self.vertical_adjustment_applied = False
        
        # 添加锚点调整风格记录
        self.anchor_style = None  # 'center' 或 'carbon'
        self.anchor_style_label = None  # 用于显示和导出
        
        # 添加偏移量记录
        self.x_offset = 0.0
        self.y_offset = 0.0
        
        # Connect signals to slots
        self.connect_signals()
        
        # Set initial UI state
        self.update_ui_state()
        
        # Show the main window
        self.show()
    
    def initUI(self):
        """Initialize the UI components"""
        self.setWindowTitle(self.lang.get_text('window_title'))
        self.setMinimumSize(1200, 800)
        
        # Create central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Create main layout
        main_layout = QHBoxLayout(central_widget)
        
        # Left control panel
        control_panel = QWidget()
        control_layout = QVBoxLayout(control_panel)
        control_panel.setMaximumWidth(500)
        
        # File operations section
        file_group = QGroupBox(self.lang.get_text('file_operations'))
        self.file_group = file_group  # 保存组件引用以供update_language使用
        file_layout = QVBoxLayout()
        
        # Molecule file selection
        mol_layout = QHBoxLayout()
        self.mol_label = QLabel(self.lang.get_text('molecule_file'))
        self.mol_path = QLineEdit()
        self.mol_browse = QPushButton(self.lang.get_text('browse'))
        self.mol_browse.clicked.connect(self.browse_molecule)
        mol_layout.addWidget(self.mol_label)
        mol_layout.addWidget(self.mol_path)
        mol_layout.addWidget(self.mol_browse)
        file_layout.addLayout(mol_layout)
        
        # Substrate file selection
        sub_layout = QHBoxLayout()
        self.sub_label = QLabel(self.lang.get_text('substrate_file'))
        self.sub_path = QLineEdit()
        self.sub_browse = QPushButton(self.lang.get_text('browse'))
        self.sub_browse.clicked.connect(self.browse_substrate)
        sub_layout.addWidget(self.sub_label)
        sub_layout.addWidget(self.sub_path)
        sub_layout.addWidget(self.sub_browse)
        file_layout.addLayout(sub_layout)
        
        # 添加导入/导出状态按钮
        state_layout = QHBoxLayout()
        self.import_state_btn = QPushButton(self.lang.get_text('import_state'))
        self.import_state_btn.clicked.connect(self.import_state)
        self.export_state_btn = QPushButton(self.lang.get_text('export_state'))
        self.export_state_btn.clicked.connect(self.export_state)
        state_layout.addWidget(self.import_state_btn)
        state_layout.addWidget(self.export_state_btn)
        file_layout.addLayout(state_layout)
        
        # Display style combo box
        style_layout = QHBoxLayout()
        style_layout.addWidget(QLabel(self.lang.get_text('display_style')))
        self.style_combo = QComboBox()
        self.style_combo.addItems(self.lang.get_text('display_styles'))
        self.style_combo.currentIndexChanged.connect(self.update_display)
        style_layout.addWidget(self.style_combo)
        file_layout.addLayout(style_layout)
        
        file_group.setLayout(file_layout)
        control_layout.addWidget(file_group)
        
        # Adsorption parameters section
        adsorption_group = QGroupBox(self.lang.get_text('adsorption_parameters'))
        self.adsorption_group = adsorption_group  # 保存组件引用
        adsorption_layout = QVBoxLayout()
        
        # Adsorption site selection
        site_layout = QHBoxLayout()
        self.site_label = QLabel(self.lang.get_text('adsorption_site'))
        site_layout.addWidget(self.site_label)
        self.site_combo = QComboBox()
        self.site_combo.addItems(self.lang.get_text('site_types'))
        site_layout.addWidget(self.site_combo)
        self.site_select = QPushButton(self.lang.get_text('select_site'))
        self.site_select.clicked.connect(self.open_site_selection)
        site_layout.addWidget(self.site_select)
        adsorption_layout.addLayout(site_layout)
        
        # Target element entry
        element_layout = QHBoxLayout()
        element_label = QLabel(self.lang.get_text('target_element'))
        self.element_entry = QLineEdit()
        
        # 创建一个空的widget来存放这些控件但不显示它们
        hidden_element_widget = QWidget()
        hidden_layout = QHBoxLayout(hidden_element_widget)
        hidden_layout.addWidget(element_label)
        hidden_layout.addWidget(self.element_entry)
        hidden_element_widget.setVisible(False)  # 隐藏整个widget
        adsorption_layout.addWidget(hidden_element_widget)
        
        # Height setting
        height_layout = QHBoxLayout()
        self.height_label = QLabel(self.lang.get_text('distance_from_surface'))
        height_layout.addWidget(self.height_label)
        self.height = QDoubleSpinBox()
        self.height.setRange(0.1, 10.0)
        self.height.setValue(2.0)
        self.height.setSingleStep(0.1)
        height_layout.addWidget(self.height)
        adsorption_layout.addLayout(height_layout)
        
        # Vacuum layer height
        vacuum_layout = QHBoxLayout()
        self.vacuum_label = QLabel(self.lang.get_text('minimum_vacuum'))
        vacuum_layout.addWidget(self.vacuum_label)
        self.vacuum = QDoubleSpinBox()
        self.vacuum.setRange(0.0, 30.0)
        self.vacuum.setValue(10.0)
        self.vacuum.setSingleStep(1.0)
        vacuum_layout.addWidget(self.vacuum)
        adsorption_layout.addLayout(vacuum_layout)
        
        # Anchor Atom
        self.anchor_label = QLabel(self.lang.get_text('anchor_atom'))
        self.anchor = QSpinBox()
        self.anchor.setMinimum(0)
        self.anchor.setMaximum(999)
        self.anchor_element_label = QLabel("")
        self.select_anchor_btn = QPushButton(self.lang.get_text('select_by_structure'))
        self.select_anchor_btn.clicked.connect(self.open_anchor_selection)
        
        anchor_layout = QHBoxLayout()
        anchor_layout.addWidget(self.anchor_label)
        anchor_layout.addWidget(self.anchor)
        anchor_layout.addWidget(self.anchor_element_label)
        anchor_layout.addWidget(self.select_anchor_btn)
        adsorption_layout.addLayout(anchor_layout)
        
        adsorption_group.setLayout(adsorption_layout)
        control_layout.addWidget(adsorption_group)
        
        # Add molecule adjustment section
        adjust_group = QGroupBox(self.lang.get_text('molecule_adjustment'))
        self.adjust_group = adjust_group  # 保存组件引用
        adjust_layout = QVBoxLayout()
        
        # Rotation adjustment
        rotate_group = QGroupBox(self.lang.get_text('rotation_controls'))
        self.rotate_group = rotate_group  # 保存组件引用
        rotate_layout = QVBoxLayout()
        
        # X-axis rotation
        x_rotate_layout = QHBoxLayout()
        self.x_axis_label = QLabel(self.lang.get_text('x_axis'))
        x_rotate_layout.addWidget(self.x_axis_label)
        self.x_slider = QSlider(Qt.Horizontal)
        self.x_slider.setRange(0, 360)
        self.x_slider.setValue(0)
        self.x_label = QLabel("0°")
        x_rotate_layout.addWidget(self.x_slider)
        x_rotate_layout.addWidget(self.x_label)
        rotate_layout.addLayout(x_rotate_layout)
        
        # Y-axis rotation
        y_rotate_layout = QHBoxLayout()
        self.y_axis_label = QLabel(self.lang.get_text('y_axis'))
        y_rotate_layout.addWidget(self.y_axis_label)
        self.y_slider = QSlider(Qt.Horizontal)
        self.y_slider.setRange(0, 360)
        self.y_slider.setValue(0)
        self.y_label = QLabel("0°")
        y_rotate_layout.addWidget(self.y_slider)
        y_rotate_layout.addWidget(self.y_label)
        rotate_layout.addLayout(y_rotate_layout)
        
        # Z-axis rotation
        z_rotate_layout = QHBoxLayout()
        self.z_axis_label = QLabel(self.lang.get_text('z_axis'))
        z_rotate_layout.addWidget(self.z_axis_label)
        self.z_slider = QSlider(Qt.Horizontal)
        self.z_slider.setRange(0, 360)
        self.z_slider.setValue(0)
        self.z_label = QLabel("0°")
        z_rotate_layout.addWidget(self.z_slider)
        z_rotate_layout.addWidget(self.z_label)
        rotate_layout.addLayout(z_rotate_layout)
        
        rotate_group.setLayout(rotate_layout)
        adjust_layout.addWidget(rotate_group)
        
        # Translation adjustment
        translate_group = QGroupBox(self.lang.get_text('translation_adjustment'))
        self.translate_group = translate_group  # 保存组件引用
        translate_layout = QVBoxLayout()
        
        # 添加保持垂直调整的选项
        self.keep_vertical = QCheckBox(self.lang.get_text('keep_vertical'))
        self.keep_vertical.setChecked(False)
        translate_layout.addWidget(self.keep_vertical)
        
        # XY plane check and auto adjustment
        self.xy_check = QCheckBox(self.lang.get_text('check_xy_plane'))
        self.xy_check.setChecked(True)
        translate_layout.addWidget(self.xy_check)
        
        # Auto adjust position button
        self.adjust_btn = QPushButton(self.lang.get_text('auto_adjust'))
        translate_layout.addWidget(self.adjust_btn)
        
        translate_group.setLayout(translate_layout)
        adjust_layout.addWidget(translate_group)
        
        adjust_group.setLayout(adjust_layout)
        control_layout.addWidget(adjust_group)
        
        # Create action buttons
        action_layout = QHBoxLayout()
        # 修改按钮文本
        self.create_btn = QPushButton(self.lang.get_text('reset_system'))
        self.validate_btn = QPushButton(self.lang.get_text('validate_structure'))
        self.export_btn = QPushButton(self.lang.get_text('export_structure'))
        
        action_layout.addWidget(self.create_btn)
        action_layout.addWidget(self.validate_btn)
        action_layout.addWidget(self.export_btn)
        
        control_layout.addLayout(action_layout)
        
        # Add status information
        self.status_label = QLabel(self.lang.get_text('status_ready'))
        control_layout.addWidget(self.status_label)
        
        # 添加批处理区域
        batch_group = QGroupBox(self.lang.get_text('batch_processing'))
        self.batch_group = batch_group  # 保存组件引用
        batch_layout = QVBoxLayout()
        
        # 添加三个批处理按钮
        self.site_batch_btn = QPushButton(self.lang.get_text('site_batch'))
        self.anchor_batch_btn = QPushButton(self.lang.get_text('anchor_batch'))
        self.path_batch_btn = QPushButton(self.lang.get_text('path_batch'))
        
        batch_layout.addWidget(self.site_batch_btn)
        batch_layout.addWidget(self.anchor_batch_btn)
        batch_layout.addWidget(self.path_batch_btn)
        
        batch_group.setLayout(batch_layout)
        control_layout.addWidget(batch_group)
        
        # Spring: fill control panel bottom space
        control_layout.addStretch(1)
        
        # Create display area
        display_widget = QWidget()
        display_layout = QHBoxLayout(display_widget)  # 改为水平布局
        
        # 创建3D视图容器
        view3d_container = QWidget()
        view3d_layout = QVBoxLayout()
        
        # 3D显示画布
        self.canvas = MplCanvas(self, width=8, height=6, dpi=100)
        self.canvas.setMinimumSize(400, 300)
        
        # 添加保存图片按钮
        self.save_fig_btn = QPushButton(self.lang.get_text('save_current_view'))
        self.save_fig_btn.clicked.connect(self.save_current_view)
        
        # 将3D视图和保存按钮添加到布局中
        view_layout = QVBoxLayout()
        view_layout.addWidget(self.canvas)
        view_layout.addWidget(self.save_fig_btn)
        view_layout.setStretch(0, 1)  # 3D视图占据更多空间
        view_layout.setStretch(1, 0)  # 按钮占据较少空间
        
        view3d_layout.addLayout(view_layout)
        
        # 导航工具栏
        self.toolbar = NavigationToolbar(self.canvas, self)
        view3d_layout.addWidget(self.toolbar)
        
        # 添加视图控制按钮区域
        view_control_group = QGroupBox(self.lang.get_text('view_controls'))
        self.view_control_group = view_control_group  # 保存组件引用
        view_control_layout = QVBoxLayout()
        
        # 缩放控制
        zoom_layout = QHBoxLayout()
        self.zoom_in_btn = QPushButton(self.lang.get_text('zoom_in'))
        self.zoom_out_btn = QPushButton(self.lang.get_text('zoom_out'))
        self.reset_view_btn = QPushButton(self.lang.get_text('reset_view'))
        zoom_layout.addWidget(self.zoom_in_btn)
        zoom_layout.addWidget(self.zoom_out_btn)
        zoom_layout.addWidget(self.reset_view_btn)
        view_control_layout.addLayout(zoom_layout)
        
        # 视角控制
        view_layout = QHBoxLayout()
        self.view_direction_label = QLabel(self.lang.get_text('view_direction'))
        view_layout.addWidget(self.view_direction_label)
        self.view_direction = QComboBox()
        self.view_direction.addItems(self.lang.get_text('view_directions'))
        view_layout.addWidget(self.view_direction)
        view_control_layout.addLayout(view_layout)
        
        view_control_group.setLayout(view_control_layout)
        view3d_layout.addWidget(view_control_group)
        
        view3d_container.setLayout(view3d_layout)
        
        # 设置3D视图的大小策略
        view3d_container.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        # 投影视图容器
        projection_container = QWidget()
        projection_layout = QVBoxLayout()
        projection_layout.setContentsMargins(0, 0, 0, 0)  # 移除内边距
        
        # 创建三个独立的投影画布
        self.xy_canvas = ProjectionCanvas(self, "XY", width=10, height=6, with_legend=True)  # 顶部画布较大，带图例
        self.xz_canvas = ProjectionCanvas(self, "XZ", width=10, height=4, with_legend=False)  # 中间画布
        self.yz_canvas = ProjectionCanvas(self, "YZ", width=10, height=4, with_legend=False)  # 底部画布
        
        projection_layout.addWidget(self.xy_canvas)
        projection_layout.addWidget(self.xz_canvas)
        projection_layout.addWidget(self.yz_canvas)
        
        projection_container.setLayout(projection_layout)
        
        # 设置投影视图容器的大小策略
        projection_container.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        projection_container.setFixedWidth(400)  # 固定投影视图宽度
        
        # 创建一个包含3D视图和投影视图的水平布局
        content_layout = QHBoxLayout()
        content_layout.setSpacing(20)  # 设置3D视图和投影视图之间的间距
        content_layout.addWidget(view3d_container)
        content_layout.addWidget(projection_container)
        
        # 创建一个容器来包含这个布局
        content_widget = QWidget()
        content_widget.setLayout(content_layout)
        content_widget.setMinimumWidth(1200)  # 增加最小总宽度，防止过度压缩
        
        # 将左侧控制面板和内容区域添加到主布局
        main_layout = QHBoxLayout()
        main_layout.addWidget(control_panel)
        main_layout.addWidget(content_widget, 1)
        
        # 创建中央窗口部件
        central_widget = QWidget()
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)
        
        # 设置窗口初始大小和最小大小
        self.setGeometry(100, 100, 1600, 900)  # 增加初始窗口大小
        self.setMinimumSize(1400, 800)  # 设置窗口最小大小
        
        # 添加仅旋转分子按钮区域
        rotation_buttons_layout = QHBoxLayout()
        self.reset_rotation_btn = QPushButton(self.lang.get_text('reset_pose'))
        rotation_buttons_layout.addWidget(self.reset_rotation_btn)
        rotate_layout.addLayout(rotation_buttons_layout)
        
        rotate_group.setLayout(rotate_layout)
        adjust_layout.addWidget(rotate_group)
        
        # 添加XY偏移量控制组
        offset_group = QGroupBox(self.lang.get_text('xy_offset'))
        self.offset_group = offset_group  # 保存组件引用
        offset_layout = QVBoxLayout()
        
        # X偏移
        x_offset_layout = QHBoxLayout()
        self.x_offset_label = QLabel(self.lang.get_text('x_offset'))
        x_offset_layout.addWidget(self.x_offset_label)
        self.x_offset_spin = QDoubleSpinBox()
        self.x_offset_spin.setRange(-10.0, 10.0)
        self.x_offset_spin.setValue(0.0)
        self.x_offset_spin.setSingleStep(0.1)
        self.x_offset_spin.setDecimals(2)
        x_offset_layout.addWidget(self.x_offset_spin)
        offset_layout.addLayout(x_offset_layout)
        
        # Y偏移
        y_offset_layout = QHBoxLayout()
        self.y_offset_label = QLabel(self.lang.get_text('y_offset'))
        y_offset_layout.addWidget(self.y_offset_label)
        self.y_offset_spin = QDoubleSpinBox()
        self.y_offset_spin.setRange(-10.0, 10.0)
        self.y_offset_spin.setValue(0.0)
        self.y_offset_spin.setSingleStep(0.1)
        self.y_offset_spin.setDecimals(2)
        y_offset_layout.addWidget(self.y_offset_spin)
        offset_layout.addLayout(y_offset_layout)
        
        offset_group.setLayout(offset_layout)
        adjust_layout.addWidget(offset_group)
        
    def connect_signals(self):
        """连接所有信号"""
        # 保留其他信号连接
        self.style_combo.currentIndexChanged.connect(self.update_display)
        
        # Adsorption parameters
        self.site_select.clicked.connect(self.open_site_selection)
        
        # 连接批处理按钮信号
        self.site_batch_btn.clicked.connect(self.batch_site_process)
        self.anchor_batch_btn.clicked.connect(self.batch_anchor_process)
        self.path_batch_btn.clicked.connect(self.batch_path_process)
        
        # 连接Distance from Surface滑块的信号
        self.height.valueChanged.connect(self.update_distance_from_surface)
        
        # 连接Minimum Vacuum Height滑块的信号
        self.vacuum.valueChanged.connect(self.update_minimum_vacuum)
        
        # 连接锚点原子变化的信号
        self.anchor.valueChanged.connect(self.update_anchor_atom)
        
        # Molecule adjustment
        self.x_slider.valueChanged.connect(self.update_rotation_label)
        self.y_slider.valueChanged.connect(self.update_rotation_label)
        self.z_slider.valueChanged.connect(self.update_rotation_label)
        
        # Apply rotation when any slider changes
        self.x_slider.valueChanged.connect(self.apply_rotation)
        self.y_slider.valueChanged.connect(self.apply_rotation)
        self.z_slider.valueChanged.connect(self.apply_rotation)
        
        # 连接重置旋转按钮
        self.reset_rotation_btn.clicked.connect(self.reset_rotation)
        
        # 连接保持垂直调整复选框
        self.keep_vertical.stateChanged.connect(self.update_vertical_adjustment_state)
        
        # 连接XY偏移控件
        self.x_offset_spin.valueChanged.connect(self.apply_offset)
        self.y_offset_spin.valueChanged.connect(self.apply_offset)
        
        # Buttons
        self.adjust_btn.clicked.connect(self.auto_adjust_position)
        self.create_btn.clicked.connect(self.create_adsorption_system)
        self.validate_btn.clicked.connect(self.validate_adsorption_system)
        self.export_btn.clicked.connect(self.export_structure)
        
        # View controls
        self.zoom_in_btn.clicked.connect(self.zoom_in)
        self.zoom_out_btn.clicked.connect(self.zoom_out)
        self.reset_view_btn.clicked.connect(self.reset_view)
        self.view_direction.currentIndexChanged.connect(self.change_view_direction)
        
        # Check box for XY plane checking
        self.xy_check.stateChanged.connect(self.update_projection_view)
    
    def update_rotation_label(self):
        """Update rotation angle labels"""
        self.x_label.setText(self.lang.get_text('rotation_angle').format(self.x_slider.value()))
        self.y_label.setText(self.lang.get_text('rotation_angle').format(self.y_slider.value()))
        self.z_label.setText(self.lang.get_text('rotation_angle').format(self.z_slider.value()))
    
    def update_ui_state(self):
        """Update UI state based on loaded data"""
        molecule_loaded = self.molecule is not None
        substrate_loaded = self.substrate is not None
        system_created = hasattr(self, 'adsorption_system') and self.adsorption_system is not None
        
        # Update file paths
        if molecule_loaded and hasattr(self, 'molecule_file'):
            self.mol_path.setText(self.molecule_file)
        else:
            self.mol_path.setText("")
        
        if substrate_loaded and hasattr(self, 'substrate_file'):
            self.sub_path.setText(self.substrate_file)
        else:
            self.sub_path.setText("")
        
        # Enable/disable buttons based on state
        self.site_select.setEnabled(substrate_loaded)
        self.element_entry.setEnabled(substrate_loaded)
        self.create_btn.setEnabled(molecule_loaded and substrate_loaded)
        
        self.validate_btn.setEnabled(system_created)
        self.export_btn.setEnabled(system_created)
        self.adjust_btn.setEnabled(system_created)
        
        # Enable/disable rotation controls
        rotation_enabled = system_created
        self.x_slider.setEnabled(rotation_enabled)
        self.y_slider.setEnabled(rotation_enabled)
        self.z_slider.setEnabled(rotation_enabled)
        
        # Update anchor atom if molecule is loaded
        if molecule_loaded:
            max_atom = len(self.molecule) - 1
            self.anchor.setMaximum(max_atom)
            # Update element label if anchor is within valid range
            if self.anchor.value() <= max_atom:
                element = self.molecule.get_chemical_symbols()[self.anchor.value()]
                self.anchor_element_label.setText(f"({element})")
            else:
                self.anchor_element_label.setText("")
        else:
            self.anchor_element_label.setText("")
        
        # Update display
        if system_created:
            self.plot_structure(self.adsorption_system)
        elif substrate_loaded and molecule_loaded:
            # Show both separately
            combined = self.substrate.copy()
            combined += self.molecule
            self.plot_structure(combined)
        elif substrate_loaded:
            self.plot_structure(self.substrate)
        elif molecule_loaded:
            self.plot_structure(self.molecule)
        else:
            # Clear canvas
            self.canvas.figure.clear()
            self.canvas.draw()
        
        # Update projection view
        self.update_projection_view()
    
    def browse_molecule(self):
        """Browse and load molecule file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, self.lang.get_text('file_dialog_title').format('Molecule'), "", 
            self.lang.get_text('file_dialog_filter')
        )
        
        if not file_path:
            return
        
        try:
            self.molecule = read(file_path)
            self.molecule_file = file_path
            self.mol_path.setText(file_path)
            self.status_label.setText(self.lang.get_text('status_loaded_molecule').format(len(self.molecule)))
            
            # 保存原始分子结构
            self.original_molecule = self.molecule.copy()
            self.initial_molecule = None  # 重置初始分子结构
            
            # Update UI
            self.update_ui_state()
            
            # Update anchor range
            max_atom = len(self.molecule) - 1
            self.anchor.setMaximum(max_atom)
            if self.anchor.value() <= max_atom:
                element = self.molecule.get_chemical_symbols()[self.anchor.value()]
                self.anchor_element_label.setText(self.lang.get_text('anchor_element').format(element))
            
            # 如果基底已加载，自动创建吸附系统
            if self.substrate is not None:
                self.create_adsorption_system()
            
        except Exception as e:
            QMessageBox.warning(self, self.lang.get_text('loading_error'), 
                              self.lang.get_text('error_loading_molecule').format(str(e)))
    
    def browse_substrate(self):
        """Browse and load substrate file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, self.lang.get_text('file_dialog_title').format('Substrate'), "", 
            self.lang.get_text('file_dialog_filter')
        )
        
        if not file_path:
            return
        
        try:
            self.substrate = read(file_path)
            self.substrate_file = file_path
            self.sub_path.setText(file_path)
            self.status_label.setText(self.lang.get_text('status_loaded_substrate').format(len(self.substrate)))
            
            # Update UI
            self.update_ui_state()
            
            # 如果分子已加载，自动创建吸附系统
            if self.molecule is not None:
                self.create_adsorption_system()
            
        except Exception as e:
            QMessageBox.warning(self, self.lang.get_text('loading_error'), 
                              self.lang.get_text('error_loading_substrate').format(str(e)))
    
    def open_site_selection(self):
        """打开吸附位点选择对话框"""
        if self.substrate is None:
            QMessageBox.warning(self, "警告", "请先加载基底结构")
            return
        
        # 断开旧的连接以防重复
        try:
            self.site_select.clicked.disconnect(self.open_site_selection)
        except:
            pass
            
        # 获取位点类型
        site_type = self.site_combo.currentText().lower()
        
        # 获取目标元素（如果指定）
        element = self.element_entry.text().strip()
        if not element:
            element = None
        
        # 创建并打开对话框
        dialog = SiteSelectionDialog(self.substrate, site_type, element, self)
        
        # 模态显示对话框
        result = dialog.exec_()
        
        # 处理结果
        if result == QDialog.Accepted and dialog.selected_site is not None:
            self.selected_site = dialog.selected_site
            self.status_label.setText(f"状态: 已选择{site_type}位点 {self.selected_site}")
            
            # 如果同时有分子和基底，可以显示预览
            if self.molecule is not None:
                # 创建临时组合系统用于可视化
                combined = self.substrate.copy()
                
                # 在选定位点放置分子
                height = self.height.value()
                anchor = self.anchor.value()
                
                try:
                    # 获取当前的旋转角度以保持旋转状态
                    x_angle = self.x_slider.value()
                    y_angle = self.y_slider.value()
                    z_angle = self.z_slider.value()
                    
                    # 从原始分子或initial_molecule开始，这取决于是否应用了垂直调整
                    if hasattr(self, 'initial_molecule') and self.initial_molecule is not None and self.vertical_adjustment_applied:
                        # 使用initial_molecule（已应用了垂直调整的基础分子）
                        temp_mol = self.initial_molecule.copy()
                    else:
                        # 使用原始分子
                        temp_mol = self.molecule.copy()
                    
                    # 如果之前应用了旋转，应用相同的旋转
                    if x_angle != 0 or y_angle != 0 or z_angle != 0:
                        # 使用锚点为中心进行旋转
                        anchor_pos = temp_mol[anchor].position.copy()
                        
                        # 将锚点移动到原点
                        temp_mol.positions -= anchor_pos
                        
                        # 创建旋转对象并应用旋转
                        rot = Rotation.from_euler('xyz', [x_angle, y_angle, z_angle], degrees=True)
                        temp_mol.positions = rot.apply(temp_mol.positions)
                        
                        # 将锚点移回原位置
                        temp_mol.positions += anchor_pos
                    
                    # 放置在指定位点
                    positioned_mol = position_molecule_on_site(
                        temp_mol, 
                        self.substrate, 
                        self.selected_site, 
                        height=height, 
                        anchor_atom=anchor
                    )
                    
                    # 不再重新应用垂直调整，因为temp_mol已经包含了所有姿态信息
                    
                    # 将分子添加到预览系统
                    combined += positioned_mol
                    
                    # 更新当前的预览系统，使其成为内部存储的当前系统
                    self.adsorption_system = combined.copy()
                    self.original_positions = self.adsorption_system.get_positions().copy()
                    
                    # 检查并调整真空高度
                    self.check_and_adjust_vacuum_height()
                    
                    # 更新3D显示
                    self.plot_structure(self.adsorption_system)
                    
                    # 清除当前投影视图并强制重新创建
                    self.xy_canvas.clear_plot()
                    self.xz_canvas.clear_plot()
                    self.yz_canvas.clear_plot()
                    
                    # 更新投影视图，确保2D视图被刷新
                    self.update_projection_view()
                    
                except Exception as e:
                    QMessageBox.warning(self, "定位错误", 
                                       f"定位分子时出错: {str(e)}")
                    import traceback
                    traceback.print_exc()
            
        # 重新连接信号
        self.site_select.clicked.connect(self.open_site_selection)
    
    def create_adsorption_system(self):
        """Create adsorption system with the current parameters"""
        if not self.molecule or not self.substrate:
            QMessageBox.warning(self, self.lang.get_text('warning'), self.lang.get_text('please_load_substrate_and_molecule_first'))
            return
        
        try:
            # 如果没有选择位点，提示用户
            if not hasattr(self, 'selected_site') or self.selected_site is None:
                QMessageBox.warning(self, self.lang.get_text('warning'), self.lang.get_text('please_select_adsorption_site_first'))
                return

            # 使用原始分子作为基础
            if not hasattr(self, 'original_molecule') or self.original_molecule is None:
                self.original_molecule = self.molecule.copy()
            
            # 获取参数值
            height = self.height.value()
            vacuum = self.vacuum.value()
            anchor = self.anchor.value()
            
            # 获取当前滑块的旋转角度
            x_angle = self.x_slider.value()
            y_angle = self.y_slider.value()
            z_angle = self.z_slider.value()
                
            # 获取当前偏移量
            x_offset = self.x_offset if hasattr(self, 'x_offset') else 0.0
            y_offset = self.y_offset if hasattr(self, 'y_offset') else 0.0
            
            # 获取垂直调整风格
            anchor_style = "none"
            if hasattr(self, 'anchor_style'):
                anchor_style = self.anchor_style
            
            self.status_label.setText(self.lang.get_text('status_processing'))
            
            # 调用自定义函数创建吸附系统
            system, is_valid, issues = create_adsorption_system_custom(
                self.original_molecule, 
                self.substrate,
                self.selected_site,
                height=height,
                x_angle=x_angle,
                y_angle=y_angle,
                z_angle=z_angle,
                x_offset=x_offset,
                y_offset=y_offset,
                anchor_atom=anchor,
                anchor_style=anchor_style,
                min_vacuum_above=vacuum,
                check_position=self.xy_check.isChecked(),
                validate=True,
                adjust_vacuum=True
            )
            
            # 保存系统并更新显示
            self.adsorption_system = system
            self.original_positions = system.positions.copy()
                
            # 显示问题（如果有）
            if not is_valid:
                issue_text = "\n".join(issues)
                QMessageBox.warning(
                    self, 
                    self.lang.get_text('structure_issues'), 
                    self.lang.get_text('system_issues').format(issue_text)
                )
                
                self.status_label.setText(self.lang.get_text('status_error') + ": " + ", ".join(issues))
            else:
                self.status_label.setText(self.lang.get_text('status_created_adsorption_system_msg'))
                
            # 更新UI
            self.update_ui_state()
                
            # 更新显示
            self.plot_structure(self.adsorption_system)
                
            # 更新投影视图
            self.update_projection_view()
                
        except Exception as e:
            QMessageBox.warning(
                self, 
                self.lang.get_text('error'), 
                self.lang.get_text('error_occurred_when_creating_adsorption_system').format(str(e))
            )
            self.status_label.setText(self.lang.get_text('status_error') + ": " + str(e))
            import traceback
            traceback.print_exc()
    
    def validate_adsorption_system(self):
        """Validate current adsorption system's reasonability"""
        if not self.adsorption_system:
            QMessageBox.warning(self, self.lang.get_text('error'), self.lang.get_text('create_system_first'))
            return
        
        try:
            # Validate system
            is_valid, issues = validate_system(
                self.adsorption_system,
                n_substrate_atoms=len(self.substrate),
                verbose=True
            )
            
            # Display validation results
            if is_valid:
                QMessageBox.information(
                    self, 
                    self.lang.get_text('validation_result'), 
                    self.lang.get_text('validation_passed')
                )
                self.status_label.setText(self.lang.get_text('status_completed'))
            else:
                issues_str = '\n'.join(issues)
                QMessageBox.warning(
                    self, 
                    self.lang.get_text('validation_result'), 
                    self.lang.get_text('validation_failed').format(issues_str)
                )
                self.status_label.setText(self.lang.get_text('status_error'))
            
        except Exception as e:
            QMessageBox.warning(
                self, 
                self.lang.get_text('validation_error'), 
                self.lang.get_text('error_validating_system').format(str(e))
            )
    
    def export_structure(self):
        """导出当前结构"""
        if not hasattr(self, 'adsorption_system') or self.adsorption_system is None:
            QMessageBox.warning(self, self.lang.get_text('warning'), self.lang.get_text('no_structure_to_save'))
            return
            
        # 创建导出选项对话框
        dialog = ExportOptionsDialog(self)
        if dialog.exec_() != QDialog.Accepted:
            return
            
        # 获取导出选项
        options = dialog.get_export_options()
        
        # 创建输出目录
        base_dir = "Adsorption_Results"
        os.makedirs(base_dir, exist_ok=True)
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        output_dir = os.path.join(base_dir, f"export_results_{timestamp}")
        os.makedirs(output_dir, exist_ok=True)
        
        # 获取分子和基底名称
        mole_name = os.path.splitext(os.path.basename(self.molecule_file))[0]
        slab_name = os.path.splitext(os.path.basename(self.substrate_file))[0]
        
        try:
            # 根据导出选项保存结构
            if options['export_adsorption']:
                # 保存吸附系统结构
                tail = options['format'].split('_')[0]
                output_file = os.path.join(output_dir, f"combine-{mole_name}_on_{slab_name}.{tail}")
                if options['format'] in ['vasp_cartesian', 'vasp_fractional']:
                    self.adsorption_system.write(output_file, format='vasp', direct=(options['format'] == 'vasp_fractional'))
                else:
                    self.adsorption_system.write(output_file, format=options['format'])
                print(self.lang.get_text('saved_system_to').format(output_file))
            
            if options['export_structures']:
                # 保存基底结构
                tail = options['format'].split('_')[0]
                slab_with_cell = self.adsorption_system[:len(self.substrate)]
                slab_file = os.path.join(output_dir, f"slab-{mole_name}_on_{slab_name}.{tail}")
                if options['format'] in ['vasp_cartesian', 'vasp_fractional']:
                    write(slab_file, slab_with_cell, format='vasp', direct=(options['format'] == 'vasp_fractional'))
                else:
                    write(slab_file, slab_with_cell, format=options['format'])
                print(self.lang.get_text('saved_substrate_to').format(slab_file))
                
                # 保存分子结构
                mole_with_cell = self.adsorption_system[len(self.substrate):]
                mole_file = os.path.join(output_dir, f"molecule-{mole_name}_on_{slab_name}.{tail}")
                if options['format'] in ['vasp_cartesian', 'vasp_fractional']:
                    write(mole_file, mole_with_cell, format='vasp', direct=(options['format'] == 'vasp_fractional'))
                else:
                    write(mole_file, mole_with_cell, format=options['format'])
                print(self.lang.get_text('saved_molecule_to').format(mole_file))
            
            # 保存图片
            if options['export_figure']:
                # 保存3D视图
                self.plot_structure(self.adsorption_system)
                self.canvas.figure.savefig(os.path.join(output_dir, f"{mole_name}_on_{slab_name}_3d.png"), dpi=300, bbox_inches='tight')
                
                # 保存投影图
                for proj_type in ['xy', 'xz', 'yz']:
                    canvas = ProjectionCanvas(projection_type=proj_type.upper(), with_legend=True)
                    canvas.plot_projection(self.substrate, self.adsorption_system[len(self.substrate):])
                    canvas.figure.savefig(os.path.join(output_dir, f"{mole_name}_on_{slab_name}_{proj_type}_proj.png"), dpi=300, bbox_inches='tight')
                    canvas.clear_plot()
                    del canvas
                
                print(self.lang.get_text('images_saved_to_dir').format(output_dir))
            
            # 保存JSON状态
            if options['export_json']:
                json_file = os.path.join(output_dir, f"{mole_name}_on_{slab_name}_state.json")
                self.mute_export_json(json_file)
                print(self.lang.get_text('saved_state_to').format(json_file))
            
            # 显示成功消息
            QMessageBox.information(
                self, 
                self.lang.get_text('export_success'), 
                self.lang.get_text('export_success_msg').format(os.path.abspath(output_dir))
            )
            
        except Exception as e:
            QMessageBox.critical(
                self, 
                self.lang.get_text('export_error'), 
                self.lang.get_text('export_error_msg').format(str(e))
            )
            import traceback
            print(self.lang.get_text('error') + ":")
            print(traceback.format_exc())
    
    def export_state(self):
        """导出当前吸附配置状态为JSON文件"""
        if not hasattr(self, 'substrate') or self.substrate is None:
            QMessageBox.warning(self, self.lang.get_text('error'), self.lang.get_text('warning_load_substrate_first'))
            return
            
        if not hasattr(self, 'selected_site') or self.selected_site is None:
            QMessageBox.warning(self, self.lang.get_text('error'), self.lang.get_text('please_select_adsorption_site_first'))
            return
            
        # 生成默认文件名（基于substrate文件名和时间戳）
        import os
        import datetime
        
        substrate_filename = os.path.basename(self.substrate_file) if hasattr(self, 'substrate_file') else "substrate"
        # 移除可能的文件扩展名（如.vasp）
        substrate_name = os.path.splitext(substrate_filename)[0]
        # 添加当前时间戳
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        default_filename = f"{substrate_name}_{timestamp}.json"
        
        # 打开文件保存对话框
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getSaveFileName(
            self, self.lang.get_text('export_adsorption_state'), default_filename, self.lang.get_text('json_file_filter')
        )
        
        if not file_path:
            return
            
        try:
            import json
            
            # 计算分数坐标
            fractional_coord = None
            if hasattr(self, 'selected_site') and self.selected_site is not None and hasattr(self, 'substrate'):
                cell = self.substrate.get_cell()
                # 转换x,y,z坐标到分数坐标
                cart_pos = np.array([self.selected_site[0], self.selected_site[1], self.selected_site[2]])
                # 创建3x3的坐标变换矩阵（考虑x,y,z三个方向）
                transform_matrix = np.linalg.inv(cell)[:3, :3]
                # 使用完整的xyz坐标
                xyz_pos = xyz_pos = np.array([self.selected_site[0], self.selected_site[1], self.selected_site[2]])
                frac_xyz = np.dot(transform_matrix, xyz_pos)
                fractional_coord = frac_xyz.tolist()
            
            # 使用存储的anchor_style属性
            anchor_style_str = None
            if hasattr(self, 'anchor_style') and self.anchor_style:
                if hasattr(self, 'anchor_style_label') and self.anchor_style_label:
                    anchor_style_str = self.anchor_style_label
                else:
                    # 如果没有标签，则使用直接的值
                    if self.anchor_style == "center":
                        anchor_style_str = "mass center"
                    elif self.anchor_style == "carbon":
                        anchor_style_str = "nearest carbon"
                    elif self.anchor_style == "none":
                        anchor_style_str = "keep current"
            
            # 构建状态数据
            state_data = {
                "substrate_file_path": self.substrate_file if hasattr(self, 'substrate_file') else None,
                "molecule_file_path": self.molecule_file if hasattr(self, 'molecule_file') else None,
                "adsorption_site_type": self.site_combo.currentText().lower(),
                "adsorption_coord_cartesian": list(self.selected_site) if hasattr(self, 'selected_site') else None,  # 保存完整的x,y,z坐标
                "adsorption_coord_fractional": fractional_coord,
                "adsorption_distance": self.height.value(),
                "vacuum_thickness": self.vacuum.value(),
                "anchor_style": anchor_style_str,
                "keep_vertical_enabled": self.keep_vertical.isChecked() if hasattr(self, 'keep_vertical') else False,
                
                # 保存旋转信息
                "rotation": {
                    # 角度值，便于人类阅读
                    "x": self.x_slider.value(),
                    "y": self.y_slider.value(),
                    "z": self.z_slider.value(),
                },
                
                # 保存偏移量信息
                "offset": {
                    "x": self.x_offset,
                    "y": self.y_offset
                }
            }
            
            # 保存到JSON文件
            with open(file_path, 'w', encoding='utf-8') as f:
                json.dump(state_data, f, indent=4, ensure_ascii=False)
                
            self.status_label.setText(self.lang.get_text('status_state_exported').format(file_path))
            QMessageBox.information(
                self, 
                self.lang.get_text('export_success'), 
                self.lang.get_text('state_exported_successfully').format(file_path)
            )
            
        except Exception as e:
            QMessageBox.warning(
                self, 
                self.lang.get_text('export_error'), 
                self.lang.get_text('error_occurred_when_exporting_adsorption_state').format(str(e))
            )
            import traceback
            traceback.print_exc()
    
    def update_display(self):
        """Update the structure display based on the selected style"""
        if hasattr(self, 'adsorption_system') and self.adsorption_system is not None:
            self.plot_structure(self.adsorption_system)
        elif hasattr(self, 'molecule') and self.molecule is not None:
            self.plot_structure(self.molecule)
        elif hasattr(self, 'substrate') and self.substrate is not None:
            self.plot_structure(self.substrate)
    
    def plot_structure(self, atoms):
        """Draw atomic structure"""
        if not atoms:
            print(f"{self.lang.get_text('warning')}: {self.lang.get_text('no_atoms_to_plot')}")
            return
        
        # 保存当前视角（如果存在）
        current_view = None
        if hasattr(self, 'canvas') and hasattr(self.canvas, 'axes'):
            try:
                current_elev, current_azim = self.canvas.axes.elev, self.canvas.axes.azim
                current_view = (current_elev, current_azim)
                print(self.lang.get_text('saving_current_view_msg').format(current_elev, current_azim))
            except:
                pass
        
        # 重新创建3D图形以避免渲染问题
        self.canvas.figure.clear()
        ax = self.canvas.figure.add_subplot(111, projection='3d')
        ax.set_facecolor('white')  # 设置背景颜色
        
        # 确保坐标轴有标签
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        
        # Get atom positions and element symbols
        positions = atoms.positions
        symbols = atoms.get_chemical_symbols()
        
        print(self.lang.get_text('plotting_structure_msg').format(len(atoms)))
        print(self.lang.get_text('first_few_positions_msg').format(positions[:3]))
        
        # Get selected display style
        display_style = self.style_combo.currentText()
        
        # 使用全局定义的元素颜色表
        # 不再需要创建color mapping
        
        # Determine substrate and molecule parts if applicable
        is_adsorption_system = hasattr(self, 'substrate') and self.substrate is not None and len(atoms) > len(self.substrate)
        n_substrate = len(self.substrate) if is_adsorption_system else 0
        
        # 绘制晶胞
        if hasattr(atoms, 'get_cell') and not np.allclose(atoms.get_cell(), np.zeros((3, 3))):
            cell = atoms.get_cell()
            # 创建晶胞顶点
            vertices = np.array([
                [0, 0, 0],
                [cell[0, 0], cell[0, 1], cell[0, 2]],
                [cell[0, 0] + cell[1, 0], cell[0, 1] + cell[1, 1], cell[0, 2] + cell[1, 2]],
                [cell[1, 0], cell[1, 1], cell[1, 2]],
                [0, 0, 0],  # 回到起点以闭合底面
                [cell[2, 0], cell[2, 1], cell[2, 2]],  # 到上表面的一个顶点
                [cell[0, 0] + cell[2, 0], cell[0, 1] + cell[2, 1], cell[0, 2] + cell[2, 2]],  # 下一个顶点
                [cell[0, 0] + cell[1, 0] + cell[2, 0], cell[0, 1] + cell[1, 1] + cell[2, 1], cell[0, 2] + cell[1, 2] + cell[2, 2]],  # 下一个顶点
                [cell[1, 0] + cell[2, 0], cell[1, 1] + cell[2, 1], cell[1, 2] + cell[2, 2]],  # 下一个顶点
                [cell[2, 0], cell[2, 1], cell[2, 2]],  # 回到上表面起点
            ])
            # 绘制底面轮廓
            ax.plot(vertices[:5, 0], vertices[:5, 1], vertices[:5, 2], 'k--', alpha=0.5, linewidth=1)
            # 绘制连接底面和顶面的线
            for i in range(4):
                ax.plot([vertices[i, 0], vertices[i+5, 0]], 
                        [vertices[i, 1], vertices[i+5, 1]], 
                        [vertices[i, 2], vertices[i+5, 2]], 'k--', alpha=0.5, linewidth=1)
            # 绘制顶面轮廓
            ax.plot(vertices[5:10, 0], vertices[5:10, 1], vertices[5:10, 2], 'k--', alpha=0.5, linewidth=1)
        
        # Dictionary to store legend elements for each element
        legend_elements = {}
        
        # Draw atoms based on selected style
        if display_style == "Ball and Stick":
            # Draw atoms as spheres
            for i, (pos, symbol) in enumerate(zip(positions, symbols)):
                # 获取标准元素颜色或使用默认颜色
                color = ELEMENT_COLORS.get(symbol, DEFAULT_COLOR)
                size = 50 if is_adsorption_system and i >= n_substrate else 30
                
                # Different marker for substrate and molecule if it's an adsorption system
                if is_adsorption_system:
                    if i < n_substrate:
                        # Substrate atom
                        marker = 'o'
                        edgecolor = 'black'
                        alpha = 0.7
                    else:
                        # Molecule atom
                        marker = 'o'
                        edgecolor = 'black'
                        alpha = 1.0
                else:
                    marker = 'o'
                    edgecolor = 'black'
                    alpha = 0.9
                
                ax.scatter(
                    pos[0], pos[1], pos[2], 
                    color=color, s=size, marker=marker, 
                    edgecolors=edgecolor, alpha=alpha
                )
            
            # Draw bonds
            self.draw_bonds(atoms, ax)
            
            # 为Ball and Stick模式创建固定大小的图例
            for symbol in set(symbols):
                color = ELEMENT_COLORS.get(symbol, DEFAULT_COLOR)
                legend_elements[symbol] = plt.Line2D(
                    [0], [0], 
                    marker='o', 
                    color='w', 
                    markerfacecolor=color,
                    markeredgecolor='black',
                    markersize=8,
                    linestyle=''
                )
            
        elif display_style == "Space Filling":
            # Draw atoms as larger spheres without bonds
            from ase.data import covalent_radii
            
            for i, (pos, symbol) in enumerate(zip(positions, symbols)):
                # Get atom number and covalent radius
                atom_number = atoms[i].number
                radius = covalent_radii[atom_number]
                
                # 获取标准元素颜色或使用默认颜色
                color = ELEMENT_COLORS.get(symbol, DEFAULT_COLOR)
                # 修改缩放因子，使其更合理
                size = (radius * 25)**2  # 减小缩放因子
                
                # Draw atom
                ax.scatter(
                    pos[0], pos[1], pos[2], 
                    color=color, s=size, marker='o', 
                    edgecolors='black', alpha=0.8
                )
            
            # 为Space Filling模式创建固定大小的图例元素
            for symbol in set(symbols):
                color = ELEMENT_COLORS.get(symbol, DEFAULT_COLOR)
                legend_elements[symbol] = plt.Line2D(
                    [0], [0], 
                    marker='o', 
                    color='w', 
                    markerfacecolor=color,
                    markeredgecolor='black',
                    markersize=8,
                    linestyle=''
                )
                    
        elif display_style == "Wireframe":
            # Draw atoms as small points with bonds
            for i, (pos, symbol) in enumerate(zip(positions, symbols)):
                # 获取标准元素颜色或使用默认颜色
                color = ELEMENT_COLORS.get(symbol, DEFAULT_COLOR)
                size = 20  # Smaller size for wireframe
                
                # Draw atom
                ax.scatter(
                    pos[0], pos[1], pos[2], 
                    color=color, s=size, marker='o',
                    alpha=0.8
                )
            
            # Draw all bonds (thinner for wireframe)
            self.draw_bonds(atoms, ax, wireframe=True)
            
            # 为Wireframe模式创建固定大小的图例
            for symbol in set(symbols):
                color = ELEMENT_COLORS.get(symbol, DEFAULT_COLOR)
                legend_elements[symbol] = plt.Line2D(
                    [0], [0], 
                    marker='o', 
                    color='w', 
                    markerfacecolor=color,
                    markersize=8,
                    linestyle=''
                )
        
        # Add legend - 移到不会遮挡的位置
        if legend_elements:
            ax.legend(
                legend_elements.values(), 
                legend_elements.keys(), 
                loc='upper right', 
                fontsize='xx-small',
                framealpha=0.7,
                ncol=2
            )
        
        # Calculate suitable viewpoint
        max_range = max([
            atoms.positions[:, 0].ptp(),
            atoms.positions[:, 1].ptp(),
            atoms.positions[:, 2].ptp()
        ])
        mid_x = (atoms.positions[:, 0].min() + atoms.positions[:, 0].max()) / 2
        mid_y = (atoms.positions[:, 1].min() + atoms.positions[:, 1].max()) / 2
        mid_z = (atoms.positions[:, 2].min() + atoms.positions[:, 2].max()) / 2
        
        # 设置视图范围
        ax.set_xlim(mid_x - max_range/1.8, mid_x + max_range/1.8)
        ax.set_ylim(mid_y - max_range/1.8, mid_y + max_range/1.8)
        ax.set_zlim(mid_z - max_range/1.8, mid_z + max_range/1.8)
        
        # 如果有保存的视角，恢复它；否则使用默认视角
        if current_view:
            ax.view_init(elev=current_view[0], azim=current_view[1])
            print(f"恢复之前的视角: 仰角={current_view[0]}, 方位角={current_view[1]}")
        else:
            # 强制设置视角，以确保3D图形可见 (仅首次绘制时)
            ax.view_init(elev=30, azim=30)
            print("使用默认视角: 仰角=30, 方位角=30")
        
        # 设置等比例 - 修复可能的广播错误
        try:
            ax.set_box_aspect([1.0, 1.0, 1.0])
        except Exception as e:
            print(f"Warning: Could not set box aspect: {str(e)}")
            try:
                ax.set_aspect('equal')
            except:
                pass
        
        # Set title
        if is_adsorption_system:
            # Set title for adsorption system
            ax.set_title(f"Adsorption System ({len(atoms)} atoms)")
        else:
            # Set title for molecule or substrate
            structure_type = "Molecule" if self.molecule is not None and atoms.get_chemical_formula() == self.molecule.get_chemical_formula() else "Substrate"
            ax.set_title(f"{structure_type} Structure ({len(atoms)} atoms)")
        
        # 禁用自动缩放，避免渲染问题
        ax.autoscale(enable=False)
        
        # 更新canvas变量以指向新的axes
        self.canvas.axes = ax
        
        # 确保刷新图形
        self.canvas.figure.tight_layout()
        self.canvas.draw_idle()  # 使用draw_idle而非draw
        print("Plot structure completed")

    def draw_bonds(self, atoms, ax, wireframe=False):
        """Draw bonds between atoms and show bond lengths"""
        from ase.neighborlist import NeighborList
        from ase.data import covalent_radii
        
        # 创建基于共价半径的邻居列表，略微增加半径以捕获键
        cutoffs = [covalent_radii[a.number] * 1.5 for a in atoms]
        nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl.update(atoms)
        
        # Define a smaller font for bond labels
        small_font = {'fontsize': 8}
        
        # Get positions
        positions = atoms.get_positions()
        
        # Line style settings based on display mode
        if wireframe:
            linestyle = '-'
            linewidth = 1
            alpha = 0.5
        else:
            linestyle = '-'
            linewidth = 2
            alpha = 0.7
        
        # Draw bonds between atoms
        drawn_bonds = set()  # Keep track of bonds we've already drawn
        
        if hasattr(self, 'substrate') and self.substrate is not None and len(atoms) > len(self.substrate):
            n_substrate = len(self.substrate)
            
            # 绘制基底内部的键连接
            for i in range(n_substrate):
                indices, offsets = nl.get_neighbors(i)
                
                for j, offset in zip(indices, offsets):
                    # 只绘制基底内部的键
                    if j < n_substrate:
                        # 创建唯一的键标识符 (始终使用较小的索引在前)
                        bond_id = tuple(sorted([i, j]))
                        
                        if bond_id not in drawn_bonds:
                            pos_i = positions[i]
                            pos_j = positions[j] + np.dot(offset, atoms.get_cell())
                            
                            # 计算距离
                            distance = np.linalg.norm(pos_i - pos_j)
                            
                            # 只绘制合理距离的键 (避免长距离)
                            # 使用更精确的共价半径和距离判断
                            atom_i_radius = covalent_radii[atoms[i].number]
                            atom_j_radius = covalent_radii[atoms[j].number]
                            max_bond_length = (atom_i_radius + atom_j_radius) * 1.5
                            
                            if distance < max_bond_length:
                                # 绘制基底键 - 使用灰色
                                ax.plot([pos_i[0], pos_j[0]], 
                                                    [pos_i[1], pos_j[1]], 
                                                    [pos_i[2], pos_j[2]], 
                                                    color='grey', linestyle=linestyle, 
                                                    linewidth=linewidth, alpha=alpha*0.7)
                                
                                drawn_bonds.add(bond_id)
            
            # 绘制分子-基底之间的键
            for i in range(n_substrate, len(atoms)):
                indices, offsets = nl.get_neighbors(i)
                
                for j, offset in zip(indices, offsets):
                    # 只绘制分子和基底之间的键
                    if j < n_substrate:
                        # 创建唯一的键标识符
                        bond_id = tuple(sorted([i, j]))
                        
                        if bond_id not in drawn_bonds:
                            pos_i = positions[i]
                            pos_j = positions[j] + np.dot(offset, atoms.get_cell())
                            
                            # 计算距离
                            distance = np.linalg.norm(pos_i - pos_j)
                            
                            # 使用更精确的共价半径和距离判断
                            atom_i_radius = covalent_radii[atoms[i].number]
                            atom_j_radius = covalent_radii[atoms[j].number]
                            max_bond_length = (atom_i_radius + atom_j_radius) * 1.5
                            
                            if distance < max_bond_length:
                                # 绘制键线
                                ax.plot([pos_i[0], pos_j[0]], 
                                                    [pos_i[1], pos_j[1]], 
                                                    [pos_i[2], pos_j[2]], 
                                                    color='k', linestyle=linestyle, 
                                                    linewidth=linewidth, alpha=alpha)
                                
                                # 添加距离标签（在中点位置） - 只标记分子与表面的成键
                                midpoint = (pos_i + pos_j) / 2
                                ax.text(midpoint[0], midpoint[1], midpoint[2], 
                                                    f"{distance:.2f} Å", **small_font)
                                
                                drawn_bonds.add(bond_id)
            
            # 绘制分子内部的键连接
            for i in range(n_substrate, len(atoms)):
                indices, offsets = nl.get_neighbors(i)
                
                for j, offset in zip(indices, offsets):
                    # 只绘制分子内部的键
                    if j >= n_substrate:
                        # 创建唯一的键标识符
                        bond_id = tuple(sorted([i, j]))
                        
                        if bond_id not in drawn_bonds:
                            pos_i = positions[i]
                            pos_j = positions[j] + np.dot(offset, atoms.get_cell())
                            
                            # 计算距离
                            distance = np.linalg.norm(pos_i - pos_j)
                            
                            # 使用更精确的共价半径和距离判断
                            atom_i_radius = covalent_radii[atoms[i].number]
                            atom_j_radius = covalent_radii[atoms[j].number]
                            max_bond_length = (atom_i_radius + atom_j_radius) * 1.5
                            
                            if distance < max_bond_length:
                                # 绘制分子内部键 - 蓝色
                                ax.plot([pos_i[0], pos_j[0]], 
                                                    [pos_i[1], pos_j[1]], 
                                                    [pos_i[2], pos_j[2]], 
                                                    color='b', linestyle=linestyle,
                                                    linewidth=linewidth, alpha=alpha)
                                
                                # 不添加距离标签，只标记分子与表面的成键距离
                                
                                drawn_bonds.add(bond_id)
        else:
            # 绘制单个结构时的所有键连接
            for i in range(len(atoms)):
                indices, offsets = nl.get_neighbors(i)
                
                for j, offset in zip(indices, offsets):
                    # 创建唯一的键标识符
                    bond_id = tuple(sorted([i, j]))
                    
                    if bond_id not in drawn_bonds:
                        pos_i = positions[i]
                        pos_j = positions[j] + np.dot(offset, atoms.get_cell())
                        
                        # 计算距离
                        distance = np.linalg.norm(pos_i - pos_j)
                        
                        # 使用更精确的共价半径和距离判断
                        atom_i_radius = covalent_radii[atoms[i].number]
                        atom_j_radius = covalent_radii[atoms[j].number]
                        max_bond_length = (atom_i_radius + atom_j_radius) * 1.5
                        
                        if distance < max_bond_length:
                            # 绘制键线
                            ax.plot([pos_i[0], pos_j[0]], 
                                            [pos_i[1], pos_j[1]], 
                                            [pos_i[2], pos_j[2]], 
                                            color='grey', linestyle=linestyle,
                                            linewidth=linewidth, alpha=alpha)
                            
                            drawn_bonds.add(bond_id)

    # 添加新的视图控制方法
    def zoom_in(self):
        """放大视图"""
        xlim = self.canvas.axes.get_xlim()
        ylim = self.canvas.axes.get_ylim()
        zlim = self.canvas.axes.get_zlim()
        
        # 计算中心点
        x_center = (xlim[0] + xlim[1]) / 2
        y_center = (ylim[0] + ylim[1]) / 2
        z_center = (zlim[0] + zlim[1]) / 2
        
        # 缩小范围 (放大视图)
        x_range = (xlim[1] - xlim[0]) * 0.8
        y_range = (ylim[1] - ylim[0]) * 0.8
        z_range = (zlim[1] - zlim[0]) * 0.8
        
        # 设置新的范围
        self.canvas.axes.set_xlim(x_center - x_range/2, x_center + x_range/2)
        self.canvas.axes.set_ylim(y_center - y_range/2, y_center + y_range/2)
        self.canvas.axes.set_zlim(z_center - z_range/2, z_center + z_range/2)
        
        self.canvas.draw()
    
    def zoom_out(self):
        """缩小视图"""
        xlim = self.canvas.axes.get_xlim()
        ylim = self.canvas.axes.get_ylim()
        zlim = self.canvas.axes.get_zlim()
        
        # 计算中心点
        x_center = (xlim[0] + xlim[1]) / 2
        y_center = (ylim[0] + ylim[1]) / 2
        z_center = (zlim[0] + zlim[1]) / 2
        
        # 扩大范围 (缩小视图)
        x_range = (xlim[1] - xlim[0]) * 1.25
        y_range = (ylim[1] - ylim[0]) * 1.25
        z_range = (zlim[1] - zlim[0]) * 1.25
        
        # 设置新的范围
        self.canvas.axes.set_xlim(x_center - x_range/2, x_center + x_range/2)
        self.canvas.axes.set_ylim(y_center - y_range/2, y_center + y_range/2)
        self.canvas.axes.set_zlim(z_center - z_range/2, z_center + z_range/2)
        
        self.canvas.draw()
    
    def reset_view(self):
        """重置视图到适合当前结构的范围"""
        if not hasattr(self, 'adsorption_system') and not hasattr(self, 'molecule') and not hasattr(self, 'substrate'):
            return
        
        # 确定当前显示的结构
        atoms = None
        if hasattr(self, 'adsorption_system') and self.adsorption_system is not None:
            atoms = self.adsorption_system
        elif hasattr(self, 'molecule') and self.molecule is not None:
            atoms = self.molecule
        elif hasattr(self, 'substrate') and self.substrate is not None:
            atoms = self.substrate
        
        if atoms is None:
            return
        
        # 计算合适的视图范围
        positions = atoms.positions
        max_range = max([
            positions[:, 0].ptp(),
            positions[:, 1].ptp(),
            positions[:, 2].ptp()
        ])
        mid_x = (positions[:, 0].min() + positions[:, 0].max()) / 2
        mid_y = (positions[:, 1].min() + positions[:, 1].max()) / 2
        mid_z = (positions[:, 2].min() + positions[:, 2].max()) / 2
        
        # 设置范围，略微扩大以便看清
        self.canvas.axes.set_xlim(mid_x - max_range*0.6, mid_x + max_range*0.6)
        self.canvas.axes.set_ylim(mid_y - max_range*0.6, mid_y + max_range*0.6)
        self.canvas.axes.set_zlim(mid_z - max_range*0.6, mid_z + max_range*0.6)
        
        # 重置视角到默认方向
        self.canvas.axes.view_init(elev=30, azim=30)
        
        self.canvas.draw()
    
    def change_view_direction(self):
        """根据选择更改视图方向"""
        view_option = self.view_direction.currentText()
        view_directions = self.lang.get_text('view_directions')
        
        if view_option == view_directions[0]:  # Default
            self.canvas.axes.view_init(elev=30, azim=30)
        elif view_option == view_directions[1]:  # A Axis
            self.canvas.axes.view_init(elev=0, azim=0)
        elif view_option == view_directions[2]:  # B Axis
            self.canvas.axes.view_init(elev=0, azim=90)
        elif view_option == view_directions[3]:  # C Axis
            self.canvas.axes.view_init(elev=90, azim=0)
        elif view_option == view_directions[4]:  # A-B Plane
            self.canvas.axes.view_init(elev=0, azim=45)
        elif view_option == view_directions[5]:  # B-C Plane
            self.canvas.axes.view_init(elev=45, azim=90)
        elif view_option == view_directions[6]:  # A-C Plane
            self.canvas.axes.view_init(elev=45, azim=0)
        
        self.canvas.draw()
        self.status_label.setText(self.lang.get_text('status_view_direction_changed').format(view_option))

    def open_anchor_selection(self):
        """打开锚点原子选择对话框"""
        if not self.molecule:
            QMessageBox.warning(self, self.lang.get_text('error'), self.lang.get_text('please_load_molecule_first'))
            return
        
        try:
            # 创建对话框
            dialog = AtomSelectionDialog(self.molecule)
            
            # 如果对话框被接受
            if dialog.exec_() == QDialog.Accepted:
                selected_atom = dialog.selected_atom
                adjust, adjust_type = dialog.get_adjustment_option()
                
                if selected_atom is not None:
                    # 更新UI
                    self.anchor.setValue(selected_atom)
            
                    # 获取选中原子的元素符号
                    element = self.molecule.get_chemical_symbols()[selected_atom]
                    self.anchor_element_label.setText(self.lang.get_text('anchor_element').format(element))
            
                    # 如果需要调整分子方向
                    if adjust:
                        # 更新锚点样式
                        self.update_anchor_style(adjust_type, True)
                        
                        # 调整分子方向
                        self.adjust_molecule_orientation(selected_atom, adjust_type)
                    else:
                        # 清除锚点样式
                        self.update_anchor_style("none", True)
            
        except Exception as e:
            QMessageBox.warning(
                self, 
                self.lang.get_text('selection_error'), 
                self.lang.get_text('error_selecting_anchor').format(str(e))
            )
            import traceback
            traceback.print_exc()

    def adjust_molecule_orientation(self, anchor_atom_idx, adjustment_type):
        """
        调整分子方向，使得质心/碳原子与锚点的连线垂直于表面
        
        参数:
            anchor_atom_idx: 锚点原子索引
            adjustment_type: 调整类型，"center"表示质心，"carbon"表示最近碳原子
        """
        try:
            # 如果没有原始分子结构，保存当前结构
            if not hasattr(self, 'original_molecule') or self.original_molecule is None:
                self.original_molecule = self.molecule.copy()
                
            # 确保初始分子结构存在
            if not hasattr(self, 'initial_molecule') or self.initial_molecule is None:
                self.create_initial_molecule_structure()
            
            # 标记已经应用了垂直调整
            self.vertical_adjustment_applied = True
            
            # 重新创建吸附系统
            if hasattr(self, 'adsorption_system') and self.adsorption_system is not None and hasattr(self, 'selected_site'):
                self.create_adsorption_system()
            
        except Exception as e:
            QMessageBox.warning(
                self, 
                self.lang.get_text('orientation_error'), 
                self.lang.get_text('error_adjusting_molecule_orientation').format(str(e))
            )
            import traceback
            traceback.print_exc()

    def update_projection_view(self):
        """Update all projection views"""
        if not self.molecule and not self.substrate and not self.adsorption_system:
            # 如果没有任何结构，清除所有投影视图
            self.xy_canvas.clear_plot()
            self.xz_canvas.clear_plot()
            self.yz_canvas.clear_plot()
            return
        
        try:
            # 如果存在吸附系统，优先显示它
            if self.adsorption_system is not None:
                # 分离基底和分子
                n_substrate = len(self.substrate) if self.substrate is not None else 0
                substrate_part = self.adsorption_system[:n_substrate].copy() if n_substrate > 0 else None
                molecule_part = self.adsorption_system[n_substrate:].copy() if len(self.adsorption_system) > n_substrate else None
                
                # 更新所有投影视图
                self.xy_canvas.plot_projection(substrate_part, molecule_part, check_position=self.xy_check.isChecked())
                self.xz_canvas.plot_projection(substrate_part, molecule_part)
                self.yz_canvas.plot_projection(substrate_part, molecule_part)
            
            # 如果没有吸附系统但有基底和分子
            elif self.substrate is not None and self.molecule is not None:
                self.xy_canvas.plot_projection(self.substrate, self.molecule, check_position=self.xy_check.isChecked())
                self.xz_canvas.plot_projection(self.substrate, self.molecule)
                self.yz_canvas.plot_projection(self.substrate, self.molecule)
            
            # 如果只有基底
            elif self.substrate is not None:
                self.xy_canvas.plot_projection(self.substrate, None)
                self.xz_canvas.plot_projection(self.substrate, None)
                self.yz_canvas.plot_projection(self.substrate, None)
            
            # 如果只有分子
            elif self.molecule is not None:
                self.xy_canvas.plot_projection(None, self.molecule)
                self.xz_canvas.plot_projection(None, self.molecule)
                self.yz_canvas.plot_projection(None, self.molecule)
        
        except Exception as e:
            print(f"{self.lang.get_text('error_updating_projection_views')}: {str(e)}")
            import traceback
            traceback.print_exc()

    def apply_rotation(self):
        """Apply rotation to the molecule part of the adsorption system"""
        if not self.adsorption_system:
            return
        
        # Get current rotation values
        x_rotation = self.x_slider.value()
        y_rotation = self.y_slider.value()
        z_rotation = self.z_slider.value()
        
        # Get anchor atom index
        anchor = self.anchor.value()
        
        try:
            # 确保已创建初始分子结构
            if not hasattr(self, 'initial_molecule') or self.initial_molecule is None:
                self.create_initial_molecule_structure()
                
            # 使用initial_molecule的副本进行旋转，不修改initial_molecule本身
            rotated_molecule = self.initial_molecule.copy()
            
            # 获取锚点原子位置
            anchor_atom_pos = rotated_molecule[anchor].position.copy()
            
            # 打印正常旋转前的分子信息
            #print("正常旋转 - 旋转前分子信息（锚点将移至原点）：")
            #print(f"锚点原子索引: {anchor}, 元素: {rotated_molecule.get_chemical_symbols()[anchor]}")
            #print("元素  |  X坐标  |  Y坐标  |  Z坐标")
            # 计算将锚点移至原点的坐标
            temp_positions = rotated_molecule.positions - anchor_atom_pos
            #for i, (symbol, pos) in enumerate(zip(rotated_molecule.get_chemical_symbols(), temp_positions)):
            #    print(f"{i}: {symbol} | {pos[0]:.6f} | {pos[1]:.6f} | {pos[2]:.6f}")
            
            # 创建旋转对象并应用旋转
            #print('normal rotation', x_rotation, y_rotation, z_rotation)
            rotation = Rotation.from_euler('xyz', [x_rotation, y_rotation, z_rotation], degrees=True)
            rotated_positions = rotation.apply(rotated_molecule.positions - anchor_atom_pos) + anchor_atom_pos
            rotated_molecule.positions = rotated_positions
            
            # 将旋转后的分子正确定位到吸附位点上
            if hasattr(self, 'selected_site') and self.selected_site is not None:
                # 准备目标位点位置
                site_pos = self.selected_site
                target_pos = site_pos.copy()
                adsorption_distance = self.height.value()
                target_pos[2] += adsorption_distance  # 在z方向上增加吸附距离
                
                # 添加x和y方向的偏移量
                x_offset = self.x_offset.value() if hasattr(self.x_offset, 'value') else self.x_offset
                y_offset = self.y_offset.value() if hasattr(self.y_offset, 'value') else self.y_offset
                target_pos[0] += x_offset  # x方向偏移
                target_pos[1] += y_offset  # y方向偏移
                
                # 计算需要移动的向量
                move_vector = target_pos - rotated_molecule[anchor].position
                
                # 移动整个分子到正确位置
                rotated_molecule.positions += move_vector
            
            # 构建新的吸附系统
            n_substrate = len(self.substrate)
            combined = self.substrate.copy()
            combined += rotated_molecule
            
            # 更新吸附系统
            self.adsorption_system = combined.copy()
            
            # 保存原始位置以供后续参考
            self.original_positions = self.adsorption_system.get_positions().copy()
            
            # 检查并调整真空高度
            self.check_and_adjust_vacuum_height()
            
            # 更新3D显示
            self.plot_structure(self.adsorption_system)
            
            # 更新投影视图
            self.update_projection_view()
            
            if self.keep_vertical.isChecked() and self.vertical_adjustment_applied:
                self.status_label.setText(self.lang.get_text('status_rotation_applied_with_vertical_adjustment'))
            else:
                self.status_label.setText(self.lang.get_text('status_rotation_applied'))
                
        except Exception as e:
            QMessageBox.warning(
                self, 
                self.lang.get_text('rotation_error'), 
                self.lang.get_text('error_applying_rotation').format(str(e))
            )

    def auto_adjust_position(self):
        """调整分子高度，确保与表面保持最小距离"""
        if not self.adsorption_system or not self.molecule:
            return
            
        try:
            # 获取基底部分
            substrate_atoms = self.substrate.copy()
            
            # 找到表面原子（z轴最高的原子）
            substrate_positions = substrate_atoms.positions
            max_z_substrate = np.max(substrate_positions[:, 2])
            max_z_indices = np.where(substrate_positions[:, 2] >= max_z_substrate - 0.1)[0]
            
            print(f"{self.lang.get_text('surface_atoms_highest_z')}: {max_z_substrate}")
            
            # 获取当前分子部分
            n_substrate = len(self.substrate)
            molecule_indices = list(range(n_substrate, len(self.adsorption_system)))
            molecule_part = self.adsorption_system[molecule_indices].copy()
            
            # 找到分子中z轴最低的原子
            molecule_positions = molecule_part.positions
            min_z_molecule = np.min(molecule_positions[:, 2])
            min_z_index = np.argmin(molecule_positions[:, 2])
            
            print(f"{self.lang.get_text('molecule_lowest_z_atom_index')}: {min_z_index}, {self.lang.get_text('position')}: {molecule_positions[min_z_index]}")
            
            # 计算当前分子最低点到表面的距离
            current_distance = min_z_molecule - max_z_substrate
            print(f"{self.lang.get_text('current_molecule_distance_to_surface')}: {current_distance:.4f} Å")
            
            # 如果距离小于1.5A，需要调整
            if current_distance < 1.5:
                # 计算需要的位移
                z_shift = 1.5 - current_distance
                print(f"{self.lang.get_text('need_upward_shift')}: {z_shift:.4f} Å")
                
                # 应用位移到所有分子原子
                for i, idx in enumerate(molecule_indices):
                    self.adsorption_system[idx].position[2] += z_shift
                
                # 更新原始位置
                self.original_positions = self.adsorption_system.get_positions().copy()
                
                # 更新吸附距离
                current_height = self.height.value()
                new_height = current_height + z_shift
                self.height.blockSignals(True)
                self.height.setValue(new_height)
                self.height.blockSignals(False)
                
                # 检查并调整真空高度
                self.check_and_adjust_vacuum_height()
                
                # 更新3D显示
                self.plot_structure(self.adsorption_system)
                
                # 更新投影视图
                self.update_projection_view()
                
                print(f"{self.lang.get_text('auto_adjusted_molecule_height')}: {z_shift:.4f} Å")
                self.status_label.setText(self.lang.get_text('status_molecule_adjusted_to_safe_distance'))
            else:
                print(self.lang.get_text('molecule_surface_distance_reasonable'))
                
        except Exception as e:
            print(f"{self.lang.get_text('error_auto_adjusting_position')}: {str(e)}")
            import traceback
            traceback.print_exc()

    def update_distance_from_surface(self):
        """实时更新分子到表面的距离"""
        if not self.molecule or not self.substrate or not hasattr(self, 'selected_site') or self.selected_site is None:
            return
        
        # 获取当前的Distance from Surface值
        height = self.height.value()
        
        # 获取锚点原子索引
        anchor = self.anchor.value()
        
        try:
            # 使用原始分子作为基础
            if not hasattr(self, 'original_molecule') or self.original_molecule is None:
                self.original_molecule = self.molecule.copy()
            
            # 获取当前的旋转角度
            x_angle = self.x_slider.value()
            y_angle = self.y_slider.value()
            z_angle = self.z_slider.value()
            
            # 获取当前偏移量
            x_offset = self.x_offset if hasattr(self, 'x_offset') else 0.0
            y_offset = self.y_offset if hasattr(self, 'y_offset') else 0.0
            
            # 获取垂直调整风格
            anchor_style = "none"
            if hasattr(self, 'anchor_style'):
                anchor_style = self.anchor_style
            
            # 调用自定义函数创建吸附系统
            system, is_valid, issues = create_adsorption_system_custom(
                self.original_molecule, 
                self.substrate, 
                self.selected_site, 
                height=height,
                x_angle=x_angle,
                y_angle=y_angle,
                z_angle=z_angle,
                x_offset=x_offset,
                y_offset=y_offset,
                anchor_atom=anchor,
                anchor_style=anchor_style,
                min_vacuum_above=self.vacuum.value(),
                check_position=self.xy_check.isChecked(),
                validate=False,  # 高度调整时不需要验证
                adjust_vacuum=True
            )
            
            # 保存系统并更新显示
            self.adsorption_system = system
            self.original_positions = system.positions.copy()
            
            # 更新3D显示
            self.plot_structure(self.adsorption_system)
            
            # 更新投影视图
            self.update_projection_view()
            
            # 更新状态栏
            self.status_label.setText(self.lang.get_text('status_distance_adjusted').format(height))
            
        except Exception as e:
            print(f"{self.lang.get_text('error_updating_distance')}: {str(e)}")
            import traceback
            traceback.print_exc()

    def check_and_adjust_vacuum_height(self):
        """检查并调整真空高度，确保满足最小真空层要求"""
        if not hasattr(self, 'adsorption_system') or self.adsorption_system is None:
            return
        
        try:
            # 获取新的最小真空高度
            vacuum_height = self.vacuum.value()
            
            # 获取当前晶胞
            cell = self.adsorption_system.get_cell()
            
            # 获取原子位置
            positions = self.adsorption_system.get_positions()
            
            # 找到z方向的最高原子
            max_z = np.max(positions[:, 2])
            
            # 计算当前真空高度
            current_vacuum = cell[2, 2] - max_z
            
            # 如果当前真空高度小于设定的最小值，则调整晶胞
            if current_vacuum != vacuum_height:
                # 计算需要增加的高度
                additional_height = vacuum_height - current_vacuum
                
                # 调整晶胞
                new_cell = cell.copy()
                new_cell[2, 2] += additional_height
                self.adsorption_system.set_cell(new_cell)
                
                # 输出调整信息
                print(f"晶胞高度已调整: +{additional_height:.2f}Å (新高度: {new_cell[2, 2]:.2f}Å)")
                
                return True
            return False
        
        except Exception as e:
            print(f"Error adjusting vacuum height: {str(e)}")
            import traceback
            traceback.print_exc()
            return False

    def update_minimum_vacuum(self):
        """更新最小真空高度设置"""
        if not hasattr(self, 'adsorption_system') or self.adsorption_system is None:
            return
        
        # 获取新的真空高度值
        vacuum_height = self.vacuum.value()
        
        # 调用独立的真空高度调整函数
        adjusted_system = adjust_vacuum_height(self.adsorption_system, vacuum_height)
        
        # 检查是否有调整
        original_cell = self.adsorption_system.get_cell()
        new_cell = adjusted_system.get_cell()
        
        if not np.array_equal(original_cell, new_cell):
            # 更新系统
            self.adsorption_system = adjusted_system
            
            # 更新3D显示
            self.plot_structure(self.adsorption_system)
            
            # 更新投影视图
            self.update_projection_view()
            
            # 更新状态栏
            self.status_label.setText(f"状态: 已调整晶胞高度以满足最小真空高度 {vacuum_height} Å")
        else:
            # 获取原子位置
            positions = self.adsorption_system.get_positions()
            
            # 找到z方向的最高原子
            max_z = np.max(positions[:, 2])
            
            # 计算当前真空高度
            current_vacuum = original_cell[2, 2] - max_z
            
            self.status_label.setText(f"状态: 当前真空高度 {current_vacuum:.2f} Å 已满足最小要求 {vacuum_height} Å")

    # 添加一个新方法来更新垂直调整状态
    def update_vertical_adjustment_state(self, state):
        """Update the vertical adjustment state based on checkbox"""
        if state == Qt.Checked:
            self.vertical_adjustment_applied = True
            # 如果用户勾选了复选框，并且已加载吸附系统，则立即应用垂直调整
            if hasattr(self, 'adsorption_system') and self.adsorption_system is not None:
                self.auto_adjust_position()
        else:
            # 取消垂直调整标志
            self.vertical_adjustment_applied = False
        
    def update_rotation_label(self):
        """Update rotation angle labels"""
        self.x_label.setText(self.lang.get_text('rotation_angle').format(self.x_slider.value()))
        self.y_label.setText(self.lang.get_text('rotation_angle').format(self.y_slider.value()))
        self.z_label.setText(self.lang.get_text('rotation_angle').format(self.z_slider.value()))

    # 添加重置旋转方法
    def reset_rotation(self):
        """重置所有旋转角度和偏移量到零，恢复分子原始姿态"""
        if not self.adsorption_system:
            return
            
        # 重置旋转
        # 暂时阻断信号连接，以避免多次触发apply_rotation
        old_x_connection = self.x_slider.blockSignals(True)
        old_y_connection = self.y_slider.blockSignals(True)
        old_z_connection = self.z_slider.blockSignals(True)
        
        # 重置所有旋转滑块到零
        self.x_slider.setValue(0)
        self.y_slider.setValue(0)
        self.z_slider.setValue(0)
        
        # 更新标签
        self.x_label.setText("0°")
        self.y_label.setText("0°")
        self.z_label.setText("0°")
        
        # 恢复信号连接
        self.x_slider.blockSignals(old_x_connection)
        self.y_slider.blockSignals(old_y_connection)
        self.z_slider.blockSignals(old_z_connection)
        
        # 重置偏移量
        # 暂时阻断信号连接
        old_x_offset_connection = self.x_offset_spin.blockSignals(True)
        old_y_offset_connection = self.y_offset_spin.blockSignals(True)
        
        # 重置偏移量控件到零
        self.x_offset_spin.setValue(0.0)
        self.y_offset_spin.setValue(0.0)
        
        # 恢复信号连接
        self.x_offset_spin.blockSignals(old_x_offset_connection)
        self.y_offset_spin.blockSignals(old_y_offset_connection)
        
        # 重置内部偏移量
        self.x_offset = 0.0
        self.y_offset = 0.0
        
        # 使用新的旋转逻辑，从初始分子重新创建吸附系统
        try:
            # 确保已创建初始分子结构
            if self.initial_molecule is None:
                self.create_initial_molecule_structure()
                
            # 重新创建吸附系统
            self.create_adsorption_system()
            
            self.status_label.setText(self.lang.get_text('status_reset_molecule_pose_msg'))
        except Exception as e:
            QMessageBox.warning(
                self, 
                self.lang.get_text('reset_error'), 
                self.lang.get_text('reset_molecule_pose_error_msg').format(str(e))
            )
            import traceback
            traceback.print_exc()

    def create_initial_molecule_structure(self):
        """
        创建初始分子结构（应用垂直调整）
        """
        try:
            # 使用原始分子作为基础
            if not hasattr(self, 'original_molecule') or self.original_molecule is None:
                self.original_molecule = self.molecule.copy()
            
            # 应用垂直调整逻辑
            anchor_atom = self.anchor.value()
            if not hasattr(self, 'anchor_style'):
                self.anchor_style = "none"
            
            # 使用新函数进行垂直调整
            if self.anchor_style != "none":
                try:
                    reference_type = "center" if self.anchor_style == "center" else "carbon"
                    self.initial_molecule = align_molecule_with_anchor_at_bottom(
                        self.original_molecule, 
                        anchor_atom, 
                        reference_type
                    )
                    self.vertical_adjustment_applied = True
                    
                    # 更新状态
                    if self.anchor_style == "center":
                        style_text = self.lang.get_text('anchor_perpendicular_to_molecule_center')
                    elif self.anchor_style == "carbon":
                        style_text = self.lang.get_text('anchor_perpendicular_to_C_atom')
                    else:
                        style_text = self.lang.get_text('keep_current_molecule_state')
                    
                    self.status_label.setText(self.lang.get_text('status_applied_vertical_adjustment_msg').format(style_text))
                except ValueError as e:
                    print(self.lang.get_text('vertical_adjustment_error_msg').format(str(e)))
                    self.initial_molecule = self.original_molecule.copy()
                    self.vertical_adjustment_applied = False
                    self.status_label.setText(self.lang.get_text('status_vertical_adjustment_failed_msg').format(str(e)))
            else:
                # 如果不需要垂直调整，直接使用原始分子
                self.initial_molecule = self.original_molecule.copy()
                self.vertical_adjustment_applied = False
                self.status_label.setText(self.lang.get_text('status_no_vertical_adjustment_applied_msg'))
                
        except Exception as e:
            QMessageBox.warning(
                self, 
                self.lang.get_text('error'), 
                self.lang.get_text('status_error_occurred_when_creating_initial_molecule_structure').format(str(e))
            )
            self.status_label.setText(self.lang.get_text('status_error') + ": " + str(e))
            import traceback
            traceback.print_exc()

    def apply_offset(self):
        """Apply offset to the molecule"""
        if not self.adsorption_system or not self.molecule:
            return
        
        # 获取新的偏移量
        new_x_offset = self.x_offset_spin.value()
        new_y_offset = self.y_offset_spin.value()
        
        # 计算增量偏移
        delta_x = new_x_offset - self.x_offset
        delta_y = new_y_offset - self.y_offset
        
        # 保存新的偏移量
        self.x_offset = new_x_offset
        self.y_offset = new_y_offset
        
        # 如果没有变化，直接返回
        if delta_x == 0 and delta_y == 0:
            return
            
        # 创建偏移向量
        offset_vector = np.array([delta_x, delta_y, 0.0])
        
        # 只对分子部分应用偏移
        if hasattr(self, 'adsorption_system') and self.adsorption_system is not None:
            # 获取基底原子数量
            n_substrate = len(self.substrate) if self.substrate is not None else 0
            
            # 仅对分子部分应用偏移
            for i in range(n_substrate, len(self.adsorption_system)):
                self.adsorption_system[i].position += offset_vector
        
        # 更新3D显示
        self.plot_structure(self.adsorption_system)
        self.update_projection_view()
        
        # 更新状态标签
        status_text = f"{self.lang.get_text('status_ready')}: {self.lang.get_text('xy_offset')} (X: {self.x_offset:.2f} Å, Y: {self.y_offset:.2f} Å)"
        self.status_label.setText(status_text)

    def import_state(self):
        """Import adsorption state from JSON file"""
        if not hasattr(self, 'substrate') or self.substrate is None:
            QMessageBox.warning(self, self.lang.get_text('error'), self.lang.get_text('warning_load_substrate_first'))
            return
            
        if not hasattr(self, 'selected_site') or self.selected_site is None:
            QMessageBox.warning(self, self.lang.get_text('error'), self.lang.get_text('please_select_adsorption_site_first'))
            return
            
        # 打开文件选择对话框
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(
            self, self.lang.get_text('import_adsorption_state'), "", self.lang.get_text('json_file_filter')
        )
        
        if not file_path:
            return
            
        try:
            import json
            
            # 读取JSON文件
            with open(file_path, 'r', encoding='utf-8') as f:
                state_data = json.load(f)
            
            # 解析数据
            substrate_file_path = state_data.get("substrate_file_path")
            molecule_file_path = state_data.get("molecule_file_path")
            adsorption_site_type = state_data.get("adsorption_site_type")
            adsorption_coord_cartesian = state_data.get("adsorption_coord_cartesian")
            adsorption_coord_fractional = state_data.get("adsorption_coord_fractional")
            adsorption_distance = state_data.get("adsorption_distance")
            vacuum_thickness = state_data.get("vacuum_thickness")
            anchor_style = state_data.get("anchor_style")
            keep_vertical_enabled = state_data.get("keep_vertical_enabled")
            
            # 更新UI
            self.update_ui_state()
            
            # 设置吸附位点类型
            self.site_combo.setCurrentText(adsorption_site_type)
            
            # 设置吸附距离
            self.height.setValue(adsorption_distance)
            
            # 设置真空层高度
            self.vacuum.setValue(vacuum_thickness)
            
            # 设置垂直调整风格
            if anchor_style:
                if anchor_style == "mass center":
                    self.anchor_style = "center"
                    self.anchor_style_label = "mass center"
                elif anchor_style == "nearest carbon":
                    self.anchor_style = "carbon"
                    self.anchor_style_label = "nearest carbon"
                else:
                    self.anchor_style = "none"
                    self.anchor_style_label = "keep current"
            
            # 设置垂直调整状态
            if keep_vertical_enabled is not None:
                self.keep_vertical.setChecked(keep_vertical_enabled)
            
            # 设置旋转角度
            rotation = state_data.get("rotation")
            if rotation:
                self.x_slider.setValue(rotation.get("x", 0))
                self.y_slider.setValue(rotation.get("y", 0))
                self.z_slider.setValue(rotation.get("z", 0))
            
            # 设置偏移量
            offset = state_data.get("offset")
            if offset:
                self.x_offset_spin.setValue(offset.get("x", 0.0))
                self.y_offset_spin.setValue(offset.get("y", 0.0))
            
            # 设置吸附位点
            if adsorption_coord_cartesian:
                self.selected_site = adsorption_coord_cartesian
                self.status_label.setText(self.lang.get_text('status_state_imported').format(file_path))
            elif adsorption_coord_fractional:
                # 转换分数坐标到笛卡尔坐标
                cell = self.substrate.get_cell()
                transform_matrix = np.linalg.inv(cell)[:3, :3]
                adsorption_coord_cartesian = np.dot(transform_matrix, adsorption_coord_fractional)
                adsorption_coord_cartesian = adsorption_coord_cartesian.tolist()
                
                self.selected_site = adsorption_coord_cartesian
                self.status_label.setText(self.lang.get_text('status_state_imported').format(file_path))
            else:
                QMessageBox.warning(self, self.lang.get_text('import_error'), self.lang.get_text('no_valid_adsorption_coordinates_found'))
            
            # 重新创建吸附系统
            self.create_adsorption_system()
            
        except Exception as e:
            QMessageBox.warning(self, self.lang.get_text('import_error'), self.lang.get_text('error_occurred_when_importing_adsorption_state').format(str(e)))
            import traceback
            traceback.print_exc()

    def update_anchor_atom(self):
        """Update the anchor atom when it changes"""
        if hasattr(self, 'adsorption_system') and self.adsorption_system is not None:
            # 重新创建吸附系统
            self.create_adsorption_system()
            
    def update_anchor_style(self, style, checked):
        """更新锚点样式并应用到系统"""
        if not checked:
            return
            
        if style == "center":
            self.anchor_style = "center"
            self.anchor_style_label = "mass center"
        elif style == "carbon":
            self.anchor_style = "carbon"
            self.anchor_style_label = "nearest carbon"
        elif style == "none":
            self.anchor_style = "none"
            self.anchor_style_label = "keep current"
        
        if hasattr(self, 'adsorption_system') and self.adsorption_system is not None:
            # 重新创建吸附系统
            self.create_adsorption_system()

    def batch_site_process(self):
        """处理批量位点选择"""
        try:
            from batch_dialogs import BatchSiteSelectionDialog
        except ImportError:
            from .batch_dialogs import BatchSiteSelectionDialog
            
        if self.substrate is None:
            QMessageBox.warning(self, self.lang.get_text('warning'), self.lang.get_text('warning_load_substrate_first'))
            return
            
        # 创建对话框
        dialog = BatchSiteSelectionDialog(self.substrate, self.site_combo.currentText().lower(), self)
        
        if dialog.exec_() == QDialog.Accepted:
            # 获取选中的位点和导出选项
            selected_sites = dialog.selected_indices
            export_options = dialog.get_export_options()
            
            # 创建输出目录（使用时间戳）
            base_dir = "Adsorption_Results"
            os.makedirs(base_dir, exist_ok=True)
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            output_dir = os.path.join(base_dir, f"batch_sites_{timestamp}")
            os.makedirs(output_dir, exist_ok=True)
            
            # 创建state目录
            state_dir = os.path.join(output_dir, "state")
            os.makedirs(state_dir, exist_ok=True)
            
            if export_options['export_images']:
                os.makedirs(os.path.join(output_dir, "3D"))
                os.makedirs(os.path.join(output_dir, "proj"))
                
            # 获取分子和基底名称
            mole_name = os.path.splitext(os.path.basename(self.molecule_file))[0]
            slab_name = os.path.splitext(os.path.basename(self.substrate_file))[0]
            
            # 创建进度对话框
            progress = QProgressDialog(self.lang.get_text('status_processing'), self.lang.get_text('cancel'), 0, len(selected_sites), self)
            progress.setWindowModality(Qt.WindowModal)
            progress.setMinimumDuration(0)  # 立即显示进度条
            progress.setAutoClose(True)  # 完成后自动关闭
            progress.setAutoReset(True)  # 完成后自动重置
            
            # 处理每个选中的位点
            for i, site in enumerate(selected_sites):
                if progress.wasCanceled():
                    break
                    
                try:
                    print(f"\n处理位点 {site}:")
                    # 获取原子位置作为吸附位点
                    site_pos = np.array(self.substrate.positions[site])
                    print(f"吸附位点坐标: {site_pos}")
                    
                    # 创建吸附系统
                    result = create_adsorption_system_custom(
                        self.molecule.copy(),
                        self.substrate.copy(),
                        site_pos,
                        height=self.height.value(),
                        x_angle=self.x_slider.value(),
                        y_angle=self.y_slider.value(),
                        z_angle=self.z_slider.value(),
                        x_offset=self.x_offset,
                        y_offset=self.y_offset,
                        anchor_atom=self.anchor.value(),
                        anchor_style=self.anchor_style if hasattr(self, 'anchor_style') else "none",
                        min_vacuum_above=self.vacuum.value(),
                        check_position=True,
                        validate=True,
                        adjust_vacuum=True
                    )
                    
                    # 从返回值中获取吸附系统
                    adsorption_system = result[0]
                    
                    # 获取位点元素
                    site_element = self.substrate.get_chemical_symbols()[site]
                    
                    # 检查分子是否合理
                    molecule_part = adsorption_system[len(self.substrate):]
                    
                    # 检查各个投影
                    status_xy, _, _ = is_molecule_pbc_valid(molecule_part, self.substrate, 'xy')
                    status_xz, _, _ = is_molecule_pbc_valid(molecule_part, self.substrate, 'xz')
                    status_yz, _, _ = is_molecule_pbc_valid(molecule_part, self.substrate, 'yz')
                    
                    # 如果任何一个投影显示overlap，添加标记
                    is_valid = all(status != 'overlap' for status in [status_xy, status_xz, status_yz])
                    validity_mark = "" if is_valid else "_invalid"
                    
                    # 构建输出文件名的基础部分
                    base_filename = f"{mole_name}_on_{slab_name}_site_{site_element}_{i+1}{validity_mark}"
                    
                    # 导出吸附结构
                    if export_options['export_adsorption']:
                        tail = export_options['export_format'].split('_')[0]
                        output_file = os.path.join(output_dir, f"combine-{base_filename}.{tail}")
                        print(self.lang.get_text('saved_system_to').format(output_file))
                        if export_options['export_format'] in ['vasp_cartesian', 'vasp_fractional']:
                            write(output_file, adsorption_system, format='vasp', direct=(export_options['export_format'] == 'vasp_fractional'))
                        else:
                            write(output_file, adsorption_system, format=export_options['export_format'])
                    
                    # 导出基底和分子结构
                    if export_options['export_structures']:
                        slab_with_cell = adsorption_system[:len(self.substrate)]
                        mole_with_cell = adsorption_system[len(self.substrate):]

                        # 导出基底结构
                        tail = export_options['export_format'].split('_')[0]
                        slab_file = os.path.join(output_dir, f"slab-{base_filename}.{tail}")
                        if export_options['export_format'] in ['vasp_cartesian', 'vasp_fractional']:
                            write(slab_file, slab_with_cell, format='vasp', direct=(export_options['export_format'] == 'vasp_fractional'))
                        else:
                            write(slab_file, slab_with_cell, format=export_options['export_format'])
                        
                        # 导出分子结构
                        tail = export_options['export_format'].split('_')[0]
                        mole_file = os.path.join(output_dir, f"molecule-{base_filename}.{tail}")
                        if export_options['export_format'] in ['vasp_cartesian', 'vasp_fractional']:
                            write(mole_file, mole_with_cell, format='vasp', direct=(export_options['export_format'] == 'vasp_fractional'))
                        else:
                            write(mole_file, mole_with_cell, format=export_options['export_format'])
                    
                    # 导出图像
                    if export_options['export_images']:
                        # 保存3D视图
                        self.plot_structure(adsorption_system)
                        self.canvas.figure.savefig(os.path.join(output_dir, "3D", f"{mole_name}_on_{slab_name}_site_{i+1}{validity_mark}_3d.png"), dpi=300, bbox_inches='tight')
                        
                        # 保存投影图
                        for proj_type in ['xy', 'xz', 'yz']:
                            canvas = ProjectionCanvas(projection_type=proj_type.upper())
                            canvas.plot_projection(self.substrate, molecule_part)
                            canvas.figure.savefig(os.path.join(output_dir, "proj", f"{mole_name}_on_{slab_name}_site_{i+1}{validity_mark}_{proj_type}_proj.png"), dpi=300, bbox_inches='tight')
                            canvas.figure.clear()
                            plt.close(canvas.figure)
                    
                    # 保存JSON状态
                    if export_options['export_json']:
                        json_file = os.path.join(state_dir, f"{base_filename}_state.json")
                        self.mute_export_json(json_file)
                        print(self.lang.get_text('saved_state_to').format(json_file))
                    
                    # 更新进度
                    progress.setValue(i + 1)
                    QApplication.processEvents()  # 确保UI更新
                        
                except Exception as e:
                    print(f"处理位点 {i + 1} 时出错: {str(e)}")
                    import traceback
                    print(self.lang.get_text('error'))
                    print(traceback.format_exc())
                    continue
            
            if not progress.wasCanceled():
                QMessageBox.information(self, self.lang.get_text('completed'), self.lang.get_text('batch_site_processing_completed_msg').format(output_dir))

    def batch_anchor_process(self):
        """处理批量锚点选择"""
        try:
            from batch_dialogs import BatchAnchorSelectionDialog
        except ImportError:
            from .batch_dialogs import BatchAnchorSelectionDialog
            
        if self.molecule is None:
            QMessageBox.warning(self, self.lang.get_text('warning'), self.lang.get_text('please_load_molecule_first'))
            return
            
        if self.selected_site is None:
            QMessageBox.warning(self, self.lang.get_text('warning'), self.lang.get_text('please_select_adsorption_site_first'))
            return
            
        # 创建对话框
        dialog = BatchAnchorSelectionDialog(self.molecule, self)
        
        if dialog.exec_() == QDialog.Accepted:
            # 获取选中的锚点和导出选项
            selected_anchors = dialog.selected_indices
            export_options = dialog.get_export_options()
            
            # 创建输出目录（使用时间戳）
            base_dir = "Adsorption_Results"
            os.makedirs(base_dir, exist_ok=True)
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            output_dir = os.path.join(base_dir, f"batch_anchors_{timestamp}")
            os.makedirs(output_dir, exist_ok=True)
            
            # 创建state目录
            state_dir = os.path.join(output_dir, "state")
            os.makedirs(state_dir, exist_ok=True)
            
            if export_options['export_images']:
                os.makedirs(os.path.join(output_dir, "3D"))
                os.makedirs(os.path.join(output_dir, "proj"))
                
            # 获取分子和基底名称
            mole_name = os.path.splitext(os.path.basename(self.molecule_file))[0]
            slab_name = os.path.splitext(os.path.basename(self.substrate_file))[0]

            # 创建进度对话框
            progress = QProgressDialog(self.lang.get_text('status_processing'), self.lang.get_text('cancel'), 0, len(selected_anchors), self)
            progress.setWindowModality(Qt.WindowModal)
            progress.setMinimumDuration(0)  # 立即显示进度条
            progress.setAutoClose(True)  # 完成后自动关闭
            progress.setAutoReset(True)  # 完成后自动重置
            
            # 处理每个选中的锚点
            for i, anchor in enumerate(selected_anchors):
                if progress.wasCanceled():
                    break
                    
                try:
                    progress.setLabelText(f"正在处理锚点 {i + 1}/{len(selected_anchors)}...")
                    print(f"\n处理锚点 {i + 1}...")
                    
                    # 使用当前锚点创建吸附系统
                    result = create_adsorption_system_custom(
                        self.molecule, self.substrate,
                        site_position=self.selected_site,
                        height=self.height.value(),
                        x_angle=self.x_slider.value(),
                        y_angle=self.y_slider.value(),
                        z_angle=self.z_slider.value(),
                        x_offset=self.x_offset,
                        y_offset=self.y_offset,
                        anchor_atom=anchor,
                        anchor_style=self.anchor_style if hasattr(self, 'anchor_style') else "none",
                        min_vacuum_above=self.vacuum.value(),
                        check_position=True,
                        validate=True,
                        adjust_vacuum=True
                    )
                    
                    # 从返回值中获取吸附系统
                    adsorption_system = result[0]
                    
                    # 获取锚点元素
                    anchor_element = self.molecule.get_chemical_symbols()[anchor]
                    
                    # 检查分子是否合理
                    molecule_part = adsorption_system[len(self.substrate):]
                    
                    # 检查各个投影
                    status_xy, _, _ = is_molecule_pbc_valid(molecule_part, self.substrate, 'xy')
                    status_xz, _, _ = is_molecule_pbc_valid(molecule_part, self.substrate, 'xz')
                    status_yz, _, _ = is_molecule_pbc_valid(molecule_part, self.substrate, 'yz')
                    
                    # 如果任何一个投影显示overlap，添加标记
                    is_valid = all(status != 'overlap' for status in [status_xy, status_xz, status_yz])
                    validity_mark = "" if is_valid else "_invalid"
                    
                    # 构建输出文件名的基础部分
                    base_filename = f"{mole_name}_on_{slab_name}_anchor_{anchor_element}_{i+1}{validity_mark}"
                    
                    # 导出吸附结构
                    if export_options['export_adsorption']:
                        tail = export_options['export_format'].split('_')[0]
                        output_file = os.path.join(output_dir, f"combine-{base_filename}.{tail}")
                        print(f"保存吸附结构到: {output_file}")
                        if export_options['export_format'] in ['vasp_cartesian', 'vasp_fractional']:
                            write(output_file, adsorption_system, format='vasp', direct=(export_options['export_format'] == 'vasp_fractional'))
                        else:
                            write(output_file, adsorption_system, format=export_options['export_format'])
                    
                    # 导出基底和分子结构
                    if export_options['export_structures']:
                        slab_with_cell = adsorption_system[:len(self.substrate)]
                        mole_with_cell = adsorption_system[len(self.substrate):]

                        # 导出基底结构
                        tail = export_options['export_format'].split('_')[0]
                        slab_file = os.path.join(output_dir, f"slab-{base_filename}.{tail}")
                        if export_options['export_format'] in ['vasp_cartesian', 'vasp_fractional']:
                            write(slab_file, slab_with_cell, format='vasp', direct=(export_options['export_format'] == 'vasp_fractional'))
                        else:
                            write(slab_file, slab_with_cell, format=export_options['export_format'])
                        
                        # 导出分子结构
                        tail = export_options['export_format'].split('_')[0]
                        mole_file = os.path.join(output_dir, f"molecule-{base_filename}.{tail}")
                        if export_options['export_format'] in ['vasp_cartesian', 'vasp_fractional']:
                            write(mole_file, mole_with_cell, format='vasp', direct=(export_options['export_format'] == 'vasp_fractional'))
                        else:
                            write(mole_file, mole_with_cell, format=export_options['export_format'])
                    
                    # 导出图像
                    if export_options['export_images']:
                        # 保存3D视图
                        self.plot_structure(adsorption_system)
                        self.canvas.figure.savefig(os.path.join(output_dir, "3D", f"{mole_name}_on_{slab_name}_anchor_{i+1}{validity_mark}_3d.png"), dpi=300, bbox_inches='tight')
                        
                        # 保存投影图
                        for proj_type in ['xy', 'xz', 'yz']:
                            canvas = ProjectionCanvas(projection_type=proj_type.upper())
                            canvas.plot_projection(self.substrate, molecule_part)
                            canvas.figure.savefig(os.path.join(output_dir, "proj", f"{mole_name}_on_{slab_name}_anchor_{i+1}{validity_mark}_{proj_type}_proj.png"), dpi=300, bbox_inches='tight')
                            canvas.figure.clear()
                            plt.close(canvas.figure)
                    
                    # 保存JSON状态
                    if export_options['export_json']:
                        json_file = os.path.join(state_dir, f"{base_filename}_state.json")
                        self.mute_export_json(json_file)
                        print(f"保存状态到: {json_file}")
                    
                    # 更新进度
                    progress.setValue(i + 1)
                    QApplication.processEvents()  # 确保UI更新
                        
                except Exception as e:
                    print(f"处理锚点 {i + 1} 时出错: {str(e)}")
                    import traceback
                    print("错误堆栈:")
                    print(traceback.format_exc())
                    continue
            
            if not progress.wasCanceled():
                QMessageBox.information(self, "完成", f"批量锚点处理完成，结果保存在 {output_dir} 目录下")

    def batch_path_process(self):
        """处理批量路径插值"""
        try:
            from batch_dialogs import BatchPathDialog
        except ImportError:
            from .batch_dialogs import BatchPathDialog
            
        if self.substrate is None or self.molecule is None:
            QMessageBox.warning(self, "警告", "请先加载基底和分子")
            return
            
        # 创建对话框
        dialog = BatchPathDialog(self.molecule, self.anchor.value(), self)
        
        if dialog.exec_() == QDialog.Accepted:
            # 获取路径数据和导出选项
            path_data = dialog.get_path_data()
            export_options = dialog.get_export_options()
            
            # 创建输出目录（使用时间戳）
            base_dir = "Adsorption_Results"
            os.makedirs(base_dir, exist_ok=True)
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            output_dir = os.path.join(base_dir, f"batch_paths_{timestamp}")
            os.makedirs(output_dir, exist_ok=True)
            
            # 创建state目录
            state_dir = os.path.join(output_dir, "state")
            os.makedirs(state_dir, exist_ok=True)
            
            if export_options['export_images']:
                os.makedirs(os.path.join(output_dir, "3D"))
                os.makedirs(os.path.join(output_dir, "proj"))
                
            # 获取分子和基底名称
            mole_name = os.path.splitext(os.path.basename(self.molecule_file))[0]
            slab_name = os.path.splitext(os.path.basename(self.substrate_file))[0]

            # 处理每对相邻的状态
            for i in range(len(path_data) - 1):
                
                json1, anchor1, _ = path_data[i]
                json2, anchor2, steps = path_data[i + 1]
                # 为每条路径创建进度对话框
                progress = QProgressDialog(f"正在处理路径 {i + 1}/{len(path_data) - 1}...", "取消", 0, steps + 1, self)
                progress.setWindowModality(Qt.WindowModal)
                progress.setMinimumDuration(0)  # 立即显示进度条
                progress.setAutoClose(True)  # 完成后自动关闭
                progress.setAutoReset(True)  # 完成后自动重置    
                    
                try:
                    print(f"\n处理路径 {i + 1}/{len(path_data) - 1}: {anchor1} -> {anchor2}")
                    
                    # 加载两个状态的参数
                    with open(json1, 'r') as f:
                        state1 = json.load(f)
                    with open(json2, 'r') as f:
                        state2 = json.load(f)
                    
                    # 对每个参数进行插值
                    for step in range(steps + 1):
                        if progress.wasCanceled():
                            break
                            
                        progress.setLabelText(f"正在处理路径 {i + 1}/{len(path_data) - 1} 的插值点 {step + 1}/{steps + 1}...")
                        
                        t = step / steps if steps > 0 else 0
                        
                        # 线性插值所有参数
                        height = state1['adsorption_distance'] * (1 - t) + state2['adsorption_distance'] * t
                        x_angle = state1['rotation']['x'] * (1 - t) + state2['rotation']['x'] * t
                        y_angle = state1['rotation']['y'] * (1 - t) + state2['rotation']['y'] * t
                        z_angle = state1['rotation']['z'] * (1 - t) + state2['rotation']['z'] * t
                        x_offset = state1['offset']['x'] * (1 - t) + state2['offset']['x'] * t
                        y_offset = state1['offset']['y'] * (1 - t) + state2['offset']['y'] * t

                        # 插值吸附位点坐标
                        site1 = np.array(state1['adsorption_coord_cartesian'])
                        site2 = np.array(state2['adsorption_coord_cartesian'])
                        interpolated_site = site1 * (1 - t) + site2 * t
                        
                        # 使用插值后的参数创建吸附系统
                        result = create_adsorption_system_custom(
                            self.molecule, self.substrate,
                            site_position=interpolated_site,  # 使用插值后的位点
                            height=height, 
                            x_angle=x_angle, y_angle=y_angle, z_angle=z_angle,
                            x_offset=x_offset, y_offset=y_offset,
                            anchor_atom=anchor1 if t < 0.5 else anchor2,  # 使用较近的锚点
                            anchor_style=self.anchor_style if hasattr(self, 'anchor_style') else "none",
                            min_vacuum_above=self.vacuum.value(),
                            check_position=True,
                            validate=True,
                            adjust_vacuum=True
                        )
                        
                        # 从返回值中获取吸附系统
                        adsorption_system = result[0]  # 第一个元素是吸附系统
                        
                        # 检查分子是否合理
                        molecule_part = adsorption_system[len(self.substrate):]
                        
                        # 检查各个投影
                        status_xy, _, _ = is_molecule_pbc_valid(molecule_part, self.substrate, 'xy')
                        status_xz, _, _ = is_molecule_pbc_valid(molecule_part, self.substrate, 'xz')
                        status_yz, _, _ = is_molecule_pbc_valid(molecule_part, self.substrate, 'yz')
                        
                        # 如果任何一个投影显示overlap，添加标记
                        is_valid = all(status != 'overlap' for status in [status_xy, status_xz, status_yz])
                        validity_mark = "" if is_valid else "_inproper"
                        
                        # 构建输出文件名
                        step_index = i * steps + step
                        base_filename = f"{mole_name}_on_{slab_name}_interpola_{step_index}{validity_mark}"
                        
                        # 导出吸附结构
                        if export_options['export_adsorption']:
                            tail = export_options['export_format'].split('_')[0]
                            output_file = os.path.join(output_dir, f"combine-{base_filename}.{tail}")
                            print(f"保存吸附结构到: {output_file}")
                            if export_options['export_format'] in ['vasp_cartesian', 'vasp_fractional']:
                                write(output_file, adsorption_system, format='vasp', direct=(export_options['export_format'] == 'vasp_fractional'))
                            else:
                                write(output_file, adsorption_system, format=export_options['export_format'])
                        
                        # 导出基底和分子结构
                        if export_options['export_structures']:
                            slab_with_cell = adsorption_system[:len(self.substrate)]
                            #slab_with_cell.set_cell(adsorption_system.get_cell())
                            #slab_with_cell.set_pbc(adsorption_system.get_pbc())
                            
                            mole_with_cell = adsorption_system[len(self.substrate):]
                            #mole_with_cell.set_cell(adsorption_system.get_cell())
                            #mole_with_cell.set_pbc(adsorption_system.get_pbc())
                            
                            # 导出基底结构
                            tail = export_options['export_format'].split('_')[0]
                            slab_file = os.path.join(output_dir, f"slab-{base_filename}.{tail}")
                            if export_options['export_format'] in ['vasp_cartesian', 'vasp_fractional']:
                                write(slab_file, slab_with_cell, format='vasp', direct=(export_options['export_format'] == 'vasp_fractional'))
                            else:
                                write(slab_file, slab_with_cell, format=export_options['export_format'])
                            
                            # 导出分子结构
                            tail = export_options['export_format'].split('_')[0]
                            mole_file = os.path.join(output_dir, f"molecule-{base_filename}.{tail}")
                            if export_options['export_format'] in ['vasp_cartesian', 'vasp_fractional']:
                                write(mole_file, mole_with_cell, format='vasp', direct=(export_options['export_format'] == 'vasp_fractional'))
                            else:
                                write(mole_file, mole_with_cell, format=export_options['export_format'])
                        
                        # 导出图像
                        if export_options['export_images']:
                            # 保存3D视图
                            self.plot_structure(adsorption_system)
                            self.canvas.figure.savefig(os.path.join(output_dir, "3D", f"{mole_name}_on_{slab_name}_interpola_{step_index}{validity_mark}_3d.png"), dpi=300, bbox_inches='tight')
                            
                            # 保存投影图
                            for proj_type in ['xy', 'xz', 'yz']:
                                canvas = ProjectionCanvas(projection_type=proj_type.upper())
                                canvas.plot_projection(self.substrate, molecule_part)
                                canvas.figure.savefig(os.path.join(output_dir, "proj", f"{mole_name}_on_{slab_name}_interpola_{step_index}{validity_mark}_{proj_type}_proj.png"), dpi=300, bbox_inches='tight')
                                canvas.figure.clear()
                                plt.close(canvas.figure)
                        
                        # 保存JSON状态
                        if export_options['export_json']:
                            json_file = os.path.join(state_dir, f"{base_filename}_state.json")
                            self.mute_export_json(json_file)
                            print(f"保存状态到: {json_file}")

                        # 更新进度
                        progress.setValue(step + 1)
                        QApplication.processEvents()  # 确保UI更新
                    
                    if progress.wasCanceled():
                        break
                        
                except Exception as e:
                    print(f"处理路径 {i + 1} 时出错: {str(e)}")
                    import traceback
                    print("错误堆栈:")
                    print(traceback.format_exc())
                    continue
            
            if not progress.wasCanceled():
                QMessageBox.information(self, "完成", f"批量路径处理完成，结果保存在 {output_dir} 目录下")

    def save_current_view(self):
        """保存当前3D视图和投影图"""
        if not hasattr(self, 'adsorption_system') or self.adsorption_system is None:
            QMessageBox.warning(self, self.lang.get_text('warning'), self.lang.get_text('no_structure_to_save'))
            return
            
        # 获取文件名前缀
        prefix, ok = QInputDialog.getText(self, self.lang.get_text('save_current_view'), self.lang.get_text('base_name'))
        if not ok or not prefix:
            return
            
        # 创建输出目录
        base_dir = "Adsorption_Results"
        os.makedirs(base_dir, exist_ok=True)
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        output_dir = os.path.join(base_dir, f"uniq_figures_{timestamp}")
        os.makedirs(output_dir, exist_ok=True)
            
        # 保存3D视图
        self.plot_structure(self.adsorption_system)
        self.canvas.figure.savefig(os.path.join(output_dir, f"{prefix}_3d.png"), dpi=300, bbox_inches='tight')
        
        # 保存投影图
        molecule_part = self.adsorption_system[len(self.substrate):]
        for proj_type in ['xy', 'xz', 'yz']:
            canvas = ProjectionCanvas(projection_type=proj_type.upper())
            canvas.plot_projection(self.substrate, molecule_part)
            canvas.figure.savefig(os.path.join(output_dir, f"{prefix}_{proj_type}_proj.png"), dpi=300, bbox_inches='tight')
            canvas.figure.clear()
            plt.close(canvas.figure)
            
        QMessageBox.information(self, self.lang.get_text('completed'), self.lang.get_text('images_saved_to_dir').format(output_dir))

    def mute_export_json(self, json_file):
        """导出当前吸附配置状态为JSON文件，不进行目录管理"""
        if not hasattr(self, 'substrate') or self.substrate is None:
            return
        try:
            import json
            
            # 计算分数坐标
            fractional_coord = None
            if hasattr(self, 'selected_site') and self.selected_site is not None and hasattr(self, 'substrate'):
                cell = self.substrate.get_cell()
                # 转换x,y,z坐标到分数坐标
                cart_pos = np.array([self.selected_site[0], self.selected_site[1], self.selected_site[2]])
                # 创建3x3的坐标变换矩阵（考虑x,y,z三个方向）
                transform_matrix = np.linalg.inv(cell)[:3, :3]
                # 使用完整的xyz坐标
                xyz_pos = xyz_pos = np.array([self.selected_site[0], self.selected_site[1], self.selected_site[2]])
                frac_xyz = np.dot(transform_matrix, xyz_pos)
                fractional_coord = frac_xyz.tolist()
            
            # 使用存储的anchor_style属性
            anchor_style_str = None
            if hasattr(self, 'anchor_style') and self.anchor_style:
                if hasattr(self, 'anchor_style_label') and self.anchor_style_label:
                    anchor_style_str = self.anchor_style_label
                else:
                    # 如果没有标签，则使用直接的值
                    if self.anchor_style == "center":
                        anchor_style_str = "mass center"
                    elif self.anchor_style == "carbon":
                        anchor_style_str = "nearest carbon"
                    elif self.anchor_style == "none":
                        anchor_style_str = "keep current"
            
            # 构建状态数据
            state_data = {
                "substrate_file_path": self.substrate_file if hasattr(self, 'substrate_file') else None,
                "molecule_file_path": self.molecule_file if hasattr(self, 'molecule_file') else None,
                "adsorption_site_type": self.site_combo.currentText().lower(),
                "adsorption_coord_cartesian": list(self.selected_site) if hasattr(self, 'selected_site') else None,  # 保存完整的x,y,z坐标
                "adsorption_coord_fractional": fractional_coord,
                "adsorption_distance": self.height.value(),
                "vacuum_thickness": self.vacuum.value(),
                "anchor_style": anchor_style_str,
                "keep_vertical_enabled": self.keep_vertical.isChecked() if hasattr(self, 'keep_vertical') else False,
                
                # 保存旋转信息
                "rotation": {
                    # 角度值，便于人类阅读
                    "x": self.x_slider.value(),
                    "y": self.y_slider.value(),
                    "z": self.z_slider.value(),
                },
                
                # 保存偏移量信息
                "offset": {
                    "x": self.x_offset,
                    "y": self.y_offset
                }
            }

                   # 保存到JSON文件
            with open(json_file, 'w', encoding='utf-8') as f:
                json.dump(state_data, f, indent=4, ensure_ascii=False)
                
            self.status_label.setText(self.lang.get_text('status_state_exported').format(json_file))
            #QMessageBox.information(self, "导出成功", f"吸附状态已成功导出到:\n{json_file}")
            
        except Exception as e:
            #QMessageBox.warning(self, "导出错误", f"导出状态时发生错误: {str(e)}")
            import traceback
            print(self.lang.get_text('error'))
            print(traceback.format_exc())

    def create_menu_bar(self):
        """创建菜单栏"""
        menubar = self.menuBar()
        
        # 文件菜单
        file_menu = menubar.addMenu(TRANSLATIONS[self.current_language]['file'])
        
        # 导入分子
        import_molecule_action = QAction(TRANSLATIONS[self.current_language]['import_molecule'], self)
        import_molecule_action.triggered.connect(self.browse_molecule)
        file_menu.addAction(import_molecule_action)
        
        # 导入基底
        import_substrate_action = QAction(TRANSLATIONS[self.current_language]['import_substrate'], self)
        import_substrate_action.triggered.connect(self.browse_substrate)
        file_menu.addAction(import_substrate_action)
        
        file_menu.addSeparator()
        
        # 导入状态
        import_state_action = QAction(TRANSLATIONS[self.current_language]['import_state'], self)
        import_state_action.triggered.connect(self.import_state)
        file_menu.addAction(import_state_action)
        
        # 导出状态
        export_state_action = QAction(TRANSLATIONS[self.current_language]['export_state'], self)
        export_state_action.triggered.connect(self.export_state)
        file_menu.addAction(export_state_action)
        
        file_menu.addSeparator()
        
        # 导出结构
        export_structure_action = QAction(TRANSLATIONS[self.current_language]['export_structure'], self)
        export_structure_action.triggered.connect(self.export_structure)
        file_menu.addAction(export_structure_action)
        
        # 视图菜单
        view_menu = menubar.addMenu(TRANSLATIONS[self.current_language]['view'])
        
        # 保存当前视图
        save_view_action = QAction(TRANSLATIONS[self.current_language]['save_current_view'], self)
        save_view_action.triggered.connect(self.save_current_view)
        view_menu.addAction(save_view_action)
        
        # 重置视图
        reset_view_action = QAction(TRANSLATIONS[self.current_language]['reset_view'], self)
        reset_view_action.triggered.connect(self.reset_view)
        view_menu.addAction(reset_view_action)
        
        # 语言菜单
        language_menu = menubar.addMenu(TRANSLATIONS[self.current_language]['language'])
        
        # 英语
        english_action = QAction('English', self)
        english_action.triggered.connect(lambda: self.change_language('en'))
        language_menu.addAction(english_action)
        
        # 中文
        chinese_action = QAction('中文', self)
        chinese_action.triggered.connect(lambda: self.change_language('zh'))
        language_menu.addAction(chinese_action)

    def change_language(self, language):
        """
        切换界面语言
        
        参数:
            language (str): 目标语言代码 ('en' 或 'zh')
        """
        if language in ['en', 'zh']:
            self.current_language = language
            self.update_language()

    def update_language(self):
        """更新所有UI元素的语言"""
        # 检查必要的属性是否存在
        required_attrs = [
            'file_group', 'mol_label', 'mol_browse', 'sub_label', 'sub_browse',
            'import_state_btn', 'export_state_btn', 'adsorption_group', 'site_label',
            'site_select', 'height_label', 'vacuum_label', 'anchor_label',
            'select_anchor_btn', 'adjust_group', 'rotate_group', 'x_axis_label',
            'y_axis_label', 'z_axis_label', 'translate_group', 'keep_vertical',
            'xy_check', 'adjust_btn', 'create_btn', 'validate_btn', 'export_btn',
            'status_label', 'site_batch_btn', 'anchor_batch_btn', 'path_batch_btn',
            'view_control_group', 'zoom_in_btn', 'zoom_out_btn', 'reset_view_btn',
            'offset_group', 'x_offset_label', 'y_offset_label', 'view_direction_label'
        ]
        
        # 检查缺失的属性
        missing_attrs = [attr for attr in required_attrs if not hasattr(self, attr)]
        if missing_attrs:
            print(f"Warning: Missing UI attributes: {missing_attrs}")
            return
            
        # 更新窗口标题
        self.setWindowTitle(TRANSLATIONS[self.current_language]['window_title'])
        
        # 添加试图方向标签的更新
        self.view_direction_label.setText(TRANSLATIONS[self.current_language]['view_direction'])
        
        # 更新文件操作组
        self.file_group.setTitle(TRANSLATIONS[self.current_language]['file_operations'])
        self.mol_label.setText(TRANSLATIONS[self.current_language]['molecule_file'])
        self.mol_browse.setText(TRANSLATIONS[self.current_language]['browse'])
        self.sub_label.setText(TRANSLATIONS[self.current_language]['substrate_file'])
        self.sub_browse.setText(TRANSLATIONS[self.current_language]['browse'])
        self.import_state_btn.setText(TRANSLATIONS[self.current_language]['import_state'])
        self.export_state_btn.setText(TRANSLATIONS[self.current_language]['export_state'])
        
        # 更新吸附参数组
        self.adsorption_group.setTitle(TRANSLATIONS[self.current_language]['adsorption_parameters'])
        self.site_label.setText(TRANSLATIONS[self.current_language]['adsorption_site'])
        self.site_select.setText(TRANSLATIONS[self.current_language]['select_site'])
        self.height_label.setText(TRANSLATIONS[self.current_language]['distance_from_surface'])
        self.vacuum_label.setText(TRANSLATIONS[self.current_language]['minimum_vacuum'])
        self.anchor_label.setText(TRANSLATIONS[self.current_language]['anchor_atom'])
        self.select_anchor_btn.setText(TRANSLATIONS[self.current_language]['select_by_structure'])
        
        # 更新分子调整组
        self.adjust_group.setTitle(TRANSLATIONS[self.current_language]['molecule_adjustment'])
        self.rotate_group.setTitle(TRANSLATIONS[self.current_language]['rotation_controls'])
        self.x_axis_label.setText(TRANSLATIONS[self.current_language]['x_axis'])
        self.y_axis_label.setText(TRANSLATIONS[self.current_language]['y_axis'])
        self.z_axis_label.setText(TRANSLATIONS[self.current_language]['z_axis'])
        
        # 更新平移调整组
        self.translate_group.setTitle(TRANSLATIONS[self.current_language]['translation_adjustment'])
        self.keep_vertical.setText(TRANSLATIONS[self.current_language]['keep_vertical'])
        self.xy_check.setText(TRANSLATIONS[self.current_language]['check_xy_plane'])
        self.adjust_btn.setText(TRANSLATIONS[self.current_language]['auto_adjust'])
        
        # 更新操作按钮
        self.create_btn.setText(TRANSLATIONS[self.current_language]['reset_system'])
        self.validate_btn.setText(TRANSLATIONS[self.current_language]['validate_structure'])
        self.export_btn.setText(TRANSLATIONS[self.current_language]['export_structure'])
        
        # 更新状态标签
        self.status_label.setText(TRANSLATIONS[self.current_language]['status_ready'])
        
        # 更新批处理按钮
        self.site_batch_btn.setText(TRANSLATIONS[self.current_language]['site_batch'])
        self.anchor_batch_btn.setText(TRANSLATIONS[self.current_language]['anchor_batch'])
        self.path_batch_btn.setText(TRANSLATIONS[self.current_language]['path_batch'])
        
        # 更新视图控制组
        self.view_control_group.setTitle(TRANSLATIONS[self.current_language]['view_controls'])
        self.zoom_in_btn.setText(TRANSLATIONS[self.current_language]['zoom_in'])
        self.zoom_out_btn.setText(TRANSLATIONS[self.current_language]['zoom_out'])
        self.reset_view_btn.setText(TRANSLATIONS[self.current_language]['reset_view'])
        
        # 更新偏移组
        self.offset_group.setTitle(TRANSLATIONS[self.current_language]['xy_offset'])
        self.x_offset_label.setText(TRANSLATIONS[self.current_language]['x_offset'])
        self.y_offset_label.setText(TRANSLATIONS[self.current_language]['y_offset'])
        
        # 更新重置姿态按钮和保存视图按钮
        self.reset_rotation_btn.setText(TRANSLATIONS[self.current_language]['reset_pose'])
        self.save_fig_btn.setText(TRANSLATIONS[self.current_language]['save_current_view'])
        
        # 更新显示样式和吸附位点下拉框
        self.style_combo.clear()
        self.style_combo.addItems(TRANSLATIONS[self.current_language]['display_styles'])
        self.site_combo.clear()
        self.site_combo.addItems(TRANSLATIONS[self.current_language]['site_types'])
        
        # 更新菜单栏文本
        menubar = self.menuBar()
        if menubar and len(menubar.actions()) >= 3:
            # 获取主菜单
            file_menu = menubar.actions()[0].menu()
            view_menu = menubar.actions()[1].menu()
            lang_menu = menubar.actions()[2].menu()
            
            # 更新菜单标题
            menubar.actions()[0].setText(TRANSLATIONS[self.current_language]['file'])
            menubar.actions()[1].setText(TRANSLATIONS[self.current_language]['view'])
            menubar.actions()[2].setText(TRANSLATIONS[self.current_language]['language'])
            
            # 更新文件菜单项
            if file_menu and len(file_menu.actions()) >= 7:
                file_menu.actions()[0].setText(TRANSLATIONS[self.current_language]['import_molecule'])
                file_menu.actions()[1].setText(TRANSLATIONS[self.current_language]['import_substrate'])
                file_menu.actions()[3].setText(TRANSLATIONS[self.current_language]['import_state'])
                file_menu.actions()[4].setText(TRANSLATIONS[self.current_language]['export_state'])
                file_menu.actions()[6].setText(TRANSLATIONS[self.current_language]['export_structure'])
            
            # 更新视图菜单项
            if view_menu and len(view_menu.actions()) >= 2:
                view_menu.actions()[0].setText(TRANSLATIONS[self.current_language]['save_current_view'])
                view_menu.actions()[1].setText(TRANSLATIONS[self.current_language]['reset_view'])
        else:
            # 如果菜单栏不存在或结构不匹配，则创建新的
            print("无法更新菜单栏，重新创建")
            menubar.clear()
            self.create_menu_bar()


# 添加处理中文字体的全局设置
def setup_chinese_font():
    """设置支持中文的字体"""
    try:
        # 设置matplotlib使用黑体或者微软雅黑
        plt.rcParams['font.family'] = ['SimHei', 'Microsoft YaHei', 'sans-serif']
        # 指定后备字体
        plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
        # 修复负号显示问题
        plt.rcParams['axes.unicode_minus'] = False
        from molecular_adsorption.languages import LanguageManager
        lang = LanguageManager()
        print(lang.get_text('chinese_font_setup_completed'))
    except Exception as e:
        print(lang.get_text('error_setting_chinese_font').format(str(e)))
        # 回退到英文标签
        pass

def angle_between(v1, v2):
    """计算两个向量之间的角度"""
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def align_molecule_with_anchor_at_bottom(molecule, anchor_atom_index, reference_point_type="center"):
    """
    调整分子方向，使锚点原子在下，分子主体在上，并与Z轴对齐
    
    参数:
        molecule (ase.Atoms): 分子结构
        anchor_atom_index (int): 锚点原子索引
        reference_point_type (str): 参考点类型，"center"表示质心，"carbon"表示最近碳原子
        
    返回:
        ase.Atoms: 调整后的分子结构
    """
    adjusted_molecule = molecule.copy()
    
    # 获取锚点原子位置
    anchor_pos = adjusted_molecule[anchor_atom_index].position
    
    # 确定参考点位置
    if reference_point_type == "center":
        # 使用质心
        reference_point = adjusted_molecule.get_center_of_mass()
    elif reference_point_type == "carbon":
        # 寻找最近的碳原子
        carbon_indices = [i for i, symbol in enumerate(adjusted_molecule.get_chemical_symbols()) if symbol == 'C']
        if carbon_indices:
            # 计算所有C原子到锚点的距离
            carbon_positions = adjusted_molecule.positions[carbon_indices]
            distances = np.linalg.norm(carbon_positions - anchor_pos, axis=1)
            nearest_carbon_idx = carbon_indices[np.argmin(distances)]
            reference_point = adjusted_molecule[nearest_carbon_idx].position
        else:
            # 没有C原子，回退到质心
            reference_point = adjusted_molecule.get_center_of_mass()
    else:
        raise ValueError(f"不支持的参考点类型: {reference_point_type}")
    
    # 计算方向向量（从锚点到参考点）
    direction_vector = reference_point - anchor_pos
    direction_length = np.linalg.norm(direction_vector)
    
    if direction_length < 1e-6:
        raise ValueError("锚点与参考点过于接近，无法确定方向")
    
    # 单位化
    direction_unit = direction_vector / direction_length
    
    # 目标向量 (Z轴正方向)
    z_axis = np.array([0, 0, 1])
    
    # 计算当前向量与Z轴的夹角
    cos_angle = np.clip(np.dot(direction_unit, z_axis), -1.0, 1.0)
    angle = np.arccos(cos_angle)
    
    # 若已对齐，无需旋转
    if np.isclose(cos_angle, 1.0, atol=1e-6):
        print("分子已经与Z轴正向对齐")
        return adjusted_molecule
    
    # 若反向对齐，旋转180度（绕X轴或Y轴）
    if np.isclose(cos_angle, -1.0, atol=1e-6):
        print("分子与Z轴反向对齐，旋转180度")
        rotation = Rotation.from_euler('x', 180, degrees=True)
    else:
        # 计算旋转轴（垂直于方向向量和Z轴）
        rotation_axis = np.cross(direction_unit, z_axis)
        rotation_axis /= np.linalg.norm(rotation_axis)
                    
        # 创建旋转对象
        rotation = Rotation.from_rotvec(rotation_axis * angle)
                    
                    # 将锚点移动到原点
    adjusted_molecule.positions -= anchor_pos
    
    # 应用旋转
    adjusted_molecule.positions = rotation.apply(adjusted_molecule.positions)

                    # 将锚点移回原位置
    adjusted_molecule.positions += anchor_pos
    
    # 验证旋转结果
    if reference_point_type == "center":
        new_reference = adjusted_molecule.get_center_of_mass()
    else:
        # 重新找最近碳原子
        carbon_indices = [i for i, symbol in enumerate(adjusted_molecule.get_chemical_symbols()) if symbol == 'C']
        if carbon_indices:
            carbon_positions = adjusted_molecule.positions[carbon_indices]
            distances = np.linalg.norm(carbon_positions - adjusted_molecule[anchor_atom_index].position, axis=1)
            nearest_carbon_idx = carbon_indices[np.argmin(distances)]
            new_reference = adjusted_molecule[nearest_carbon_idx].position
        else:
            new_reference = adjusted_molecule.get_center_of_mass()
    
    # 计算新的方向向量
    new_direction = new_reference - adjusted_molecule[anchor_atom_index].position
    new_direction_unit = new_direction / np.linalg.norm(new_direction)
    
    # 验证Z分量为正
    if new_direction_unit[2] < 0:
        print("警告：旋转后方向向量Z分量为负，锚点可能在上方")
    
    # 验证与Z轴夹角
    new_cos_angle = np.dot(new_direction_unit, z_axis)
    new_angle_deg = np.degrees(np.arccos(np.clip(new_cos_angle, -1.0, 1.0)))
    
    if new_angle_deg > 5.0:
        print(f"警告：旋转后方向向量与Z轴夹角为 {new_angle_deg:.2f}°，未完全对齐")
    else:
        print(f"旋转成功：方向向量与Z轴夹角为 {new_angle_deg:.2f}°")
        
    # 验证锚点是否在分子底部
    min_z = np.min(adjusted_molecule.positions[:, 2])
    anchor_z = adjusted_molecule[anchor_atom_index].position[2]
    
    if anchor_z > min_z + 0.1:  # 允许0.1埃的误差
        print(f"警告：锚点不是最低点，锚点Z={anchor_z:.2f}，最低点Z={min_z:.2f}")
    
    return adjusted_molecule

def adjust_vacuum_height(system, min_vacuum_above=10.0):
    """
    调整晶胞高度，确保分子顶部与晶胞顶部有足够的空间
    
    参数:
        system (ase.Atoms): 原子系统
        min_vacuum_above (float): 分子顶部与晶胞顶部的最小距离(Å)
        
    返回:
        ase.Atoms: 调整后的系统
    """
    system_copy = system.copy()
    
    # 获取原子位置
    positions = system_copy.get_positions()
    
    # 获取晶胞信息
    cell = system_copy.get_cell()
    
    # 找到z方向的最高原子
    max_z = np.max(positions[:, 2])
    
    # 计算当前真空高度
    current_vacuum = cell[2, 2] - max_z
    
    # 如果当前真空高度小于设定值，则调整晶胞
    if current_vacuum != min_vacuum_above:
        additional_height = min_vacuum_above - current_vacuum
        
        # 调整晶胞
        new_cell = cell.copy()
        new_cell[2, 2] += additional_height
        system_copy.set_cell(new_cell)
        
        print(f"晶胞高度已调整: +{additional_height:.2f}Å (新高度: {new_cell[2, 2]:.2f}Å)")
    
    return system_copy

def create_adsorption_system_custom(original_molecule, substrate, site_position, 
                            height=2.0, x_angle=0, y_angle=0, z_angle=0, 
                            x_offset=0.0, y_offset=0.0, anchor_atom=0,
                            anchor_style="none", min_vacuum_above=10.0, 
                            check_position=True, validate=True, adjust_vacuum=True,
                            lang=None):
    """
    创建分子吸附系统(完整版)
    
    参数:
        original_molecule (ase.Atoms): 原始分子结构
        substrate (ase.Atoms): 基底结构
        site_position (numpy.ndarray): 吸附位置的3D坐标
        height (float): 分子底部距离表面的高度(Å)
        x_angle, y_angle, z_angle (float): XYZ三个方向的旋转角度(度)
        x_offset, y_offset (float): XY平面上的偏移量(Å)
        anchor_atom (int): 分子中的锚定原子索引
        anchor_style (str): 垂直调整风格 - "center"(质心垂直), "carbon"(最近碳垂直), "none"(不调整)
        min_vacuum_above (float): 分子顶部与晶胞顶部的最小距离(Å)
        check_position (bool): 是否检查并调整分子位置，使其在基底内
        validate (bool): 是否验证吸附结构的合理性
        adjust_vacuum (bool): 是否调整真空层高度
        lang (LanguageManager): 语言管理器实例，如果为None则创建一个新实例
        
    返回:
        tuple: (吸附系统, 是否合理, 问题列表)
    """
    # 如果没有提供语言管理器，创建一个新实例
    if lang is None:
        lang = LanguageManager()
    
    # 创建基底的副本
    substrate_copy = substrate.copy()
    
    # 第1步：创建初始分子结构（应用垂直调整）
    initial_molecule = original_molecule.copy()
    
    # 只有在需要垂直调整时才执行以下逻辑
    if anchor_style != "none":
        try:
            # 使用新函数进行垂直调整
            reference_type = "center" if anchor_style == "center" else "carbon"
            initial_molecule = align_molecule_with_anchor_at_bottom(
                initial_molecule, 
                anchor_atom, 
                reference_type
            )
        except ValueError as e:
            print(f"{lang.get_text('vertical_adjustment_error_msg').format(str(e))}")
            # 如果调整失败，使用原始分子
            initial_molecule = original_molecule.copy()
    
    # 第2步：应用旋转角度
    molecule_to_position = initial_molecule.copy()
    
    # 获取锚点原子位置
    anchor_pos = molecule_to_position[anchor_atom].position.copy()
    
    # 将锚点移动到原点
    molecule_to_position.positions -= anchor_pos
    
    # 应用旋转角度
    if x_angle != 0 or y_angle != 0 or z_angle != 0:
        # 创建旋转对象并应用旋转
        rot = Rotation.from_euler('xyz', [x_angle, y_angle, z_angle], degrees=True)
        molecule_to_position.positions = rot.apply(molecule_to_position.positions)
    
    # 将锚点移回原位置
    molecule_to_position.positions += anchor_pos
    
    # 第3步：将分子定位到指定位点（包含XY偏移）
    # 首先确保site_position是一个numpy数组
    if not isinstance(site_position, np.ndarray):
        site_position = np.array(site_position)
    
    # 打印定位前的关键信息
    print(f"{lang.get_text('positioning_info')}:")
    print(f"- {lang.get_text('anchor_atom_index')}: {anchor_atom}")
    print(f"- {lang.get_text('selected_site_coordinates')}: {site_position}")
    print(f"- {lang.get_text('current_anchor_position')}: {molecule_to_position[anchor_atom].position}")
    print(f"- {lang.get_text('height')}: {height} Å")
    print(f"- {lang.get_text('offset')}: X={x_offset}, Y={y_offset} Å")
    
    # 计算新的锚定原子位置（包含XY偏移和高度）
    new_anchor_pos = np.array([
        site_position[0] + x_offset, 
        site_position[1] + y_offset, 
        site_position[2] + height  # 在Z方向上增加吸附高度
    ])
    
    # 计算需要移动的向量
    move_vector = new_anchor_pos - molecule_to_position[anchor_atom].position
    
    # 平移所有原子
    molecule_to_position.positions += move_vector
    
    # 验证定位结果
    actual_anchor_pos = molecule_to_position[anchor_atom].position
    expected_pos = np.array([site_position[0] + x_offset, site_position[1] + y_offset, site_position[2] + height])
    position_error = np.linalg.norm(actual_anchor_pos - expected_pos)
    if position_error > 1e-10:
        print(f"{lang.get_text('warning')}: {lang.get_text('positioning_error_msg')} {position_error:.10f} Å")
        print(f"- {lang.get_text('expected_position')}: {expected_pos}")
        print(f"- {lang.get_text('actual_position')}: {actual_anchor_pos}")
    else:
        print(f"{lang.get_text('positioning_success_msg')}")
    
    # 检查是否需要调整分子在XY平面内的位置
    if check_position:
        in_bounds = is_molecule_inside_substrate_xy(molecule_to_position, substrate_copy)
        if not in_bounds:
            print(f"{lang.get_text('molecule_outside_substrate_xy_msg')}")
            
            # 获取基底的晶胞参数
            cell = substrate_copy.get_cell()
            
            # 保存原始锚点位置，确保调整后锚点仍在原位置
            original_anchor_pos = molecule_to_position[anchor_atom].position.copy()
            
            # 获取分子XY坐标
            mol_xy = molecule_to_position.positions[:, :2]
            
            # 获取XY周期性范围
            cell_x_max = cell[0, 0]
            cell_y_max = cell[1, 1]
            
            # 对每个原子进行周期性调整，确保在XY范围内
            for i, pos in enumerate(molecule_to_position.positions):
                # 周期性调整X
                while pos[0] < 0:
                    pos[0] += cell_x_max
                while pos[0] >= cell_x_max:
                    pos[0] -= cell_x_max
                    
                # 周期性调整Y
                while pos[1] < 0:
                    pos[1] += cell_y_max
                while pos[1] >= cell_y_max:
                    pos[1] -= cell_y_max
                    
                # 更新位置
                molecule_to_position.positions[i] = pos
            
            # 重置锚点位置到原始位置
            anchor_shift = original_anchor_pos - molecule_to_position[anchor_atom].position
            molecule_to_position.positions += anchor_shift
            
            # 验证锚点是否保持不变
            anchor_deviation = np.linalg.norm(molecule_to_position[anchor_atom].position - original_anchor_pos)
            if anchor_deviation > 1e-10:
                print(f"{lang.get_text('warning')}: {lang.get_text('anchor_deviation_after_adjustment_msg')} {anchor_deviation:.10f} Å")
            else:
                print(f"{lang.get_text('adjustment_completed_anchor_unchanged_msg')}")
    
    # 第4步：创建最终吸附系统
    # 合并分子和基底
    combined = substrate_copy.copy()
    combined += molecule_to_position
    
    # 检查并调整晶胞高度
    if adjust_vacuum and min_vacuum_above > 0:
        combined = adjust_vacuum_height(combined, min_vacuum_above)
    
    # 第5步：验证系统
    is_valid = True
    issues = []
    if validate:
        is_valid, issues = validate_system(
            combined, 
            n_substrate_atoms=len(substrate), 
            scale_factor=0.7,  # 使用默认值
            min_mol_surf_dist=1.5,  # 使用默认值
            verbose=False
        )
    
    return combined, is_valid, issues

def main():
    """Main function to start the application"""
    app = QApplication(sys.argv)
    
    # 设置中文字体
    setup_chinese_font()
    
    # Create and show main window
    main_window = AdsorptionGUI()
    main_window.show()
    
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()

