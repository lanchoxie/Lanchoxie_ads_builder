from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QWidget,
                             QComboBox, QPushButton, QGroupBox, QListWidget,
                             QDialogButtonBox, QMessageBox, QCheckBox, QTableWidget,
                             QTableWidgetItem, QHeaderView, QSpinBox, QFileDialog,
                             QLabel)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
import numpy as np
import os
import sys
import importlib.util

# 判断脚本运行环境来决定导入方式
if __name__ == "__main__":
    # 当直接运行时
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if parent_dir not in sys.path:
        sys.path.insert(0, parent_dir)
    
    # 尝试导入英文版组件
    try:
        from molecular_adsorption.adsorption_gui_eng import MplCanvas, ELEMENT_COLORS
        DEFAULT_COLOR = '#CCCCCC'  # 灰色作为默认色
    except ImportError:
        # 尝试导入中文版组件
        try:
            from molecular_adsorption.adsorption_gui_chn import MplCanvas, ELEMENT_COLORS
            DEFAULT_COLOR = '#CCCCCC'  # 灰色作为默认色
        except ImportError:
            # 尝试直接导入模块
            try:
                import molecular_adsorption.adsorption_gui_eng as gui_mod
                MplCanvas = gui_mod.MplCanvas
                ELEMENT_COLORS = gui_mod.ELEMENT_COLORS
                DEFAULT_COLOR = '#CCCCCC'  # 灰色作为默认色
            except ImportError:
                try:
                    import molecular_adsorption.adsorption_gui_chn as gui_mod
                    MplCanvas = gui_mod.MplCanvas
                    ELEMENT_COLORS = gui_mod.ELEMENT_COLORS
                    DEFAULT_COLOR = '#CCCCCC'  # 灰色作为默认色
                except ImportError:
                    # 最后尝试本地导入
                    try:
                        from adsorption_gui_eng import MplCanvas, ELEMENT_COLORS
                        DEFAULT_COLOR = '#CCCCCC'  # 灰色作为默认色
                    except ImportError:
                        from adsorption_gui_chn import MplCanvas, ELEMENT_COLORS
                        DEFAULT_COLOR = '#CCCCCC'  # 灰色作为默认色
else:
    # 当作为模块导入时
    try:
        # 尝试相对导入英文版
        from .adsorption_gui_eng import MplCanvas, ELEMENT_COLORS
        DEFAULT_COLOR = '#CCCCCC'  # 灰色作为默认色
    except ImportError:
        # 尝试相对导入中文版
        from .adsorption_gui_chn import MplCanvas, ELEMENT_COLORS
        DEFAULT_COLOR = '#CCCCCC'  # 灰色作为默认色

from ase.data import covalent_radii

# 定义共价半径字典
COVALENT_RADII = {symbol: radius for symbol, radius in zip(['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                                                          'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
                                                          'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                                                          'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
                                                          'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
                                                          'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
                                                          'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
                                                          'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
                                                          'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
                                                          'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
                                                          'Md', 'No', 'Lr'], covalent_radii)}

class BatchSiteSelectionDialog(QDialog):
    """批量选择吸附位点的对话框"""
    
    def __init__(self, substrate, site_type, element=None, parent=None):
        super().__init__(parent)
        self.substrate = substrate
        self.site_type = site_type
        self.element = element
        self.parent = parent
        self.selected_sites = []  # 存储选中的位点
        self.selected_indices = set()  # 存储选中的原子索引
        self.surface_only = True  # 默认只显示表面原子
        self.initUI()
        
    def initUI(self):
        """初始化UI"""
        self.setWindowTitle("批量选择吸附位点")
        self.setMinimumSize(1000, 800)
        
        # 创建主布局
        layout = QHBoxLayout()
        
        # 左侧控制面板
        control_panel = QWidget()
        control_layout = QVBoxLayout()
        control_panel.setMaximumWidth(300)
        
        # 表面选择复选框
        self.surface_check = QCheckBox("仅表面原子")
        self.surface_check.setChecked(True)
        self.surface_check.stateChanged.connect(self.update_view)
        control_layout.addWidget(self.surface_check)
        
        # 元素筛选组
        filter_group = QGroupBox("元素筛选")
        filter_layout = QVBoxLayout()
        
        # 可用元素下拉框
        self.element_combo = QComboBox()
        self.populate_available_elements()
        filter_layout.addWidget(self.element_combo)
        
        filter_group.setLayout(filter_layout)
        control_layout.addWidget(filter_group)
        
        # 选中原子列表
        selected_group = QGroupBox("已选原子")
        selected_layout = QVBoxLayout()
        
        # 添加列表显示
        self.selected_list = QListWidget()
        selected_layout.addWidget(self.selected_list)
        
        # 添加操作按钮
        btn_layout = QHBoxLayout()
        self.add_element_btn = QPushButton("添加当前元素原子")
        self.remove_element_btn = QPushButton("删除当前元素原子")
        btn_layout.addWidget(self.add_element_btn)
        btn_layout.addWidget(self.remove_element_btn)
        selected_layout.addLayout(btn_layout)
        
        btn_layout2 = QHBoxLayout()
        self.add_all_btn = QPushButton("添加所有原子")
        self.clear_all_btn = QPushButton("清空所有原子")
        btn_layout2.addWidget(self.add_all_btn)
        btn_layout2.addWidget(self.clear_all_btn)
        selected_layout.addLayout(btn_layout2)
        
        # 添加删除框框内原子按钮
        self.remove_selected_btn = QPushButton("删除框框内原子")
        selected_layout.addWidget(self.remove_selected_btn)
        
        selected_group.setLayout(selected_layout)
        control_layout.addWidget(selected_group)
        
        # 导出选项组
        export_group = QGroupBox("导出选项")
        export_layout = QVBoxLayout()
        
        # 导出吸附结构复选框
        self.export_adsorption = QCheckBox("导出吸附结构")
        self.export_adsorption.setChecked(True)
        export_layout.addWidget(self.export_adsorption)
        
        # 导出基底和分子结构复选框
        self.export_structures = QCheckBox("导出基底结构和分子结构")
        self.export_structures.setChecked(True)
        export_layout.addWidget(self.export_structures)
        
        # 导出格式选择
        format_layout = QHBoxLayout()
        format_layout.addWidget(QLabel("导出格式:"))
        self.format_combo = QComboBox()
        self.format_combo.addItems(["vasp_cartesian", "vasp_fractional", "xsf", "xyz", "cif"])
        format_layout.addWidget(self.format_combo)
        export_layout.addLayout(format_layout)
        
        # 导出图像复选框
        self.export_images = QCheckBox("导出图像")
        self.export_images.setChecked(True)
        export_layout.addWidget(self.export_images)
        
        # 导出JSON复选框
        self.export_json = QCheckBox("导出当前状态(JSON)")
        self.export_json.setChecked(True)
        export_layout.addWidget(self.export_json)
        
        export_group.setLayout(export_layout)
        control_layout.addWidget(export_group)
        
        # 确定和取消按钮
        button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        )
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        control_layout.addWidget(button_box)
        
        control_panel.setLayout(control_layout)
        
        # 右侧显示区域
        display_widget = QWidget()
        display_layout = QVBoxLayout()
        
        # 3D显示画布
        self.canvas = MplCanvas(self)
        display_layout.addWidget(self.canvas)
        
        # 导航工具栏
        toolbar = NavigationToolbar2QT(self.canvas, self)
        display_layout.addWidget(toolbar)
        
        display_widget.setLayout(display_layout)
        
        # 添加到主布局
        layout.addWidget(control_panel)
        layout.addWidget(display_widget, stretch=1)
        
        self.setLayout(layout)
        
        # 连接信号
        self.add_element_btn.clicked.connect(self.add_current_element)
        self.remove_element_btn.clicked.connect(self.remove_current_element)
        self.add_all_btn.clicked.connect(self.add_all_atoms)
        self.clear_all_btn.clicked.connect(self.clear_all_atoms)
        self.remove_selected_btn.clicked.connect(self.remove_selected_atoms)
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.selected_list.itemClicked.connect(self.on_list_item_clicked)
        
        # 绘制初始结构
        self.plot_substrate()
        
    def get_visible_indices(self):
        """获取当前可见的原子索引"""
        if self.surface_only:
            return set(self.get_surface_atoms())
        return set(range(len(self.substrate)))
        
    def update_view(self, state):
        """更新视图显示（表面/3D）"""
        self.surface_only = state == Qt.Checked
        # 更新选中列表，移除不可见的原子
        visible_indices = self.get_visible_indices()
        self.selected_indices &= visible_indices
        self.update_selected_list()
        self.plot_substrate()
        
    def get_surface_atoms(self):
        """获取表面原子的索引"""
        positions = self.substrate.positions
        cell = self.substrate.cell
        max_z = np.max(positions[:, 2])
        surface_indices = [i for i, pos in enumerate(positions) 
                         if abs(pos[2] - max_z) < 0.1]  # 0.1埃的容差
        return surface_indices
        
    def on_list_item_clicked(self, item):
        """处理列表项点击事件"""
        # 从列表项文本中提取原子索引
        text = item.text()
        idx = int(''.join(filter(str.isdigit, text)))
        # 从选中集合中移除该原子
        if idx in self.selected_indices:
            self.selected_indices.remove(idx)
            # 更新列表和显示
            self.update_selected_list()
            self.plot_substrate()
            
    def remove_selected_atoms(self):
        """删除框框内选中的原子"""
        # 获取当前选中的列表项
        selected_items = self.selected_list.selectedItems()
        if not selected_items:
            return
            
        # 从选中集合中移除这些原子
        for item in selected_items:
            text = item.text()
            idx = int(''.join(filter(str.isdigit, text)))
            self.selected_indices.discard(idx)
            
        # 更新列表和显示
        self.update_selected_list()
        self.plot_substrate()
        
    def on_click(self, event):
        """处理鼠标点击事件"""
        if event.inaxes != self.canvas.axes:
            return
            
        if not event.button == 1:  # 只处理左键点击
            return
            
        # 获取点击位置
        x, y = event.xdata, event.ydata
        
        # 找到最近的原子
        positions = self.substrate.positions
        visible_indices = self.get_visible_indices()
        visible_positions = positions[list(visible_indices)]
        
        distances = np.sqrt((visible_positions[:, 0] - x)**2 + (visible_positions[:, 1] - y)**2)
        nearest_idx = list(visible_indices)[np.argmin(distances)]
        
        # 切换选中状态
        if nearest_idx in self.selected_indices:
            self.selected_indices.remove(nearest_idx)
        else:
            self.selected_indices.add(nearest_idx)
            
        self.update_selected_list()
        self.plot_substrate()
        
    def add_current_element(self):
        """添加当前选中元素的所有可见原子"""
        element = self.element_combo.currentText()
        visible_indices = self.get_visible_indices()
        symbols = self.substrate.get_chemical_symbols()
        
        for idx in visible_indices:
            if symbols[idx] == element:
                self.selected_indices.add(idx)
                
        self.update_selected_list()
        self.plot_substrate()
        
    def remove_current_element(self):
        """删除当前选中元素的所有可见原子"""
        element = self.element_combo.currentText()
        visible_indices = self.get_visible_indices()
        symbols = self.substrate.get_chemical_symbols()
        
        for idx in visible_indices:
            if symbols[idx] == element:
                self.selected_indices.discard(idx)
                
        self.update_selected_list()
        self.plot_substrate()
        
    def add_all_atoms(self):
        """添加所有可见原子"""
        visible_indices = self.get_visible_indices()
        self.selected_indices.update(visible_indices)
        self.update_selected_list()
        self.plot_substrate()
        
    def clear_all_atoms(self):
        """清空所有选中的原子"""
        self.selected_indices.clear()
        self.update_selected_list()
        self.plot_substrate()
        
    def plot_substrate(self):
        """绘制基底结构，标记选中的原子"""
        self.canvas.figure.clear()
        ax = self.canvas.figure.add_subplot(111, projection='3d')
        
        # 获取要显示的原子索引
        visible_indices = self.get_visible_indices()
        
        # 绘制所有可见的原子
        positions = self.substrate.positions
        symbols = self.substrate.get_chemical_symbols()
        
        # 绘制未选中的原子
        for i in visible_indices:
            if i not in self.selected_indices:
                pos = positions[i]
                symbol = symbols[i]
                color = ELEMENT_COLORS.get(symbol, DEFAULT_COLOR)
                ax.scatter(pos[0], pos[1], pos[2], c=color, s=100)
        
        # 绘制选中的原子（红色）
        for i in self.selected_indices & visible_indices:
            pos = positions[i]
            ax.scatter(pos[0], pos[1], pos[2], c='red', s=100)
        
        # 绘制晶格
        if hasattr(self.substrate, 'cell'):
            cell = self.substrate.cell
            for i in range(2):
                for j in range(2):
                    ax.plot([0, cell[0][0]], [j*cell[1][1], j*cell[1][1]], [i*cell[2][2], i*cell[2][2]], 'k--', alpha=0.5)
                    ax.plot([j*cell[0][0], j*cell[0][0]], [0, cell[1][1]], [i*cell[2][2], i*cell[2][2]], 'k--', alpha=0.5)
                    ax.plot([j*cell[0][0], j*cell[0][0]], [i*cell[1][1], i*cell[1][1]], [0, cell[2][2]], 'k--', alpha=0.5)
        
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        
        self.canvas.draw()
        
    def update_selected_list(self):
        """更新选中原子列表"""
        self.selected_list.clear()
        symbols = self.substrate.get_chemical_symbols()
        for idx in sorted(self.selected_indices):
            self.selected_list.addItem(f"{symbols[idx]}{idx}")
            
    def populate_available_elements(self):
        """填充可用元素下拉框"""
        symbols = set(self.substrate.get_chemical_symbols())
        self.element_combo.clear()
        self.element_combo.addItems(sorted(symbols))
        
    def get_export_options(self):
        """获取导出选项"""
        return {
            'export_adsorption': self.export_adsorption.isChecked(),
            'export_structures': self.export_structures.isChecked(),
            'export_format': self.format_combo.currentText(),
            'export_images': self.export_images.isChecked(),
            'export_json': self.export_json.isChecked()
        }
        
    def accept(self):
        """确认选择"""
        if not self.selected_indices:
            QMessageBox.warning(self, "警告", "请至少选择一个吸附位点")
            return
        super().accept()

class BatchAnchorSelectionDialog(QDialog):
    """批量选择锚点原子的对话框"""
    
    def __init__(self, molecule, parent=None):
        super().__init__(parent)
        self.molecule = molecule
        self.parent = parent
        self.selected_indices = set()  # 存储选中的原子索引
        self.default_elements = ['O', 'N', 'F', 'Cl', 'Br', 'S', 'P']  # 默认选中的元素
        self.initUI()
        
    def initUI(self):
        """初始化UI"""
        self.setWindowTitle("批量选择锚点原子")
        self.setMinimumSize(1000, 800)
        
        # 创建主布局
        layout = QHBoxLayout()
        
        # 左侧控制面板
        control_panel = QWidget()
        control_layout = QVBoxLayout()
        control_panel.setMaximumWidth(300)
        
        # 元素筛选组
        filter_group = QGroupBox("元素筛选")
        filter_layout = QVBoxLayout()
        
        # 可用元素下拉框
        self.element_combo = QComboBox()
        self.populate_available_elements()
        filter_layout.addWidget(self.element_combo)
        
        filter_group.setLayout(filter_layout)
        control_layout.addWidget(filter_group)
        
        # 选中原子列表
        selected_group = QGroupBox("已选原子")
        selected_layout = QVBoxLayout()
        
        # 添加列表显示
        self.selected_list = QListWidget()
        selected_layout.addWidget(self.selected_list)
        
        # 添加操作按钮
        btn_layout = QHBoxLayout()
        self.add_element_btn = QPushButton("添加当前元素原子")
        self.remove_element_btn = QPushButton("删除当前元素原子")
        btn_layout.addWidget(self.add_element_btn)
        btn_layout.addWidget(self.remove_element_btn)
        selected_layout.addLayout(btn_layout)
        
        btn_layout2 = QHBoxLayout()
        self.add_all_btn = QPushButton("添加所有原子")
        self.clear_all_btn = QPushButton("清空所有原子")
        btn_layout2.addWidget(self.add_all_btn)
        btn_layout2.addWidget(self.clear_all_btn)
        selected_layout.addLayout(btn_layout2)
        
        selected_group.setLayout(selected_layout)
        control_layout.addWidget(selected_group)
        
        # 导出选项组
        export_group = QGroupBox("导出选项")
        export_layout = QVBoxLayout()
        
        # 导出吸附结构复选框
        self.export_adsorption = QCheckBox("导出吸附结构")
        self.export_adsorption.setChecked(True)
        export_layout.addWidget(self.export_adsorption)
        
        # 导出基底和分子结构复选框
        self.export_structures = QCheckBox("导出基底结构和分子结构")
        self.export_structures.setChecked(True)
        export_layout.addWidget(self.export_structures)
        
        # 导出格式选择
        format_layout = QHBoxLayout()
        format_layout.addWidget(QLabel("导出格式:"))
        self.format_combo = QComboBox()
        self.format_combo.addItems(["vasp_cartesian", "vasp_fractional", "xsf", "xyz", "cif"])
        format_layout.addWidget(self.format_combo)
        export_layout.addLayout(format_layout)
        
        # 导出图像复选框
        self.export_images = QCheckBox("导出图像")
        self.export_images.setChecked(True)
        export_layout.addWidget(self.export_images)
        
        # 导出JSON复选框
        self.export_json = QCheckBox("导出当前状态(JSON)")
        self.export_json.setChecked(True)
        export_layout.addWidget(self.export_json)
        
        export_group.setLayout(export_layout)
        control_layout.addWidget(export_group)
        
        # 确定和取消按钮
        button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        )
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        control_layout.addWidget(button_box)
        
        control_panel.setLayout(control_layout)
        
        # 右侧显示区域
        display_widget = QWidget()
        display_layout = QVBoxLayout()
        
        # 3D显示画布
        self.canvas = MplCanvas(self)
        display_layout.addWidget(self.canvas)
        
        # 导航工具栏
        toolbar = NavigationToolbar2QT(self.canvas, self)
        display_layout.addWidget(toolbar)
        
        display_widget.setLayout(display_layout)
        
        # 添加到主布局
        layout.addWidget(control_panel)
        layout.addWidget(display_widget, stretch=1)
        
        self.setLayout(layout)
        
        # 连接信号
        self.add_element_btn.clicked.connect(self.add_current_element)
        self.remove_element_btn.clicked.connect(self.remove_current_element)
        self.add_all_btn.clicked.connect(self.add_all_atoms)
        self.clear_all_btn.clicked.connect(self.clear_all_atoms)
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.selected_list.itemClicked.connect(self.on_list_item_clicked)
        
        # 绘制初始结构
        self.plot_molecule()
        
        # 默认选中指定元素
        self.select_default_elements()
        
    def select_default_elements(self):
        """选中默认元素的所有原子"""
        symbols = self.molecule.get_chemical_symbols()
        for i, symbol in enumerate(symbols):
            if symbol in self.default_elements:
                self.selected_indices.add(i)
        self.update_selected_list()
        self.plot_molecule()
        
    def on_list_item_clicked(self, item):
        """处理列表项点击事件"""
        # 从列表项文本中提取原子索引
        text = item.text()
        idx = int(''.join(filter(str.isdigit, text)))
        # 从选中集合中移除该原子
        if idx in self.selected_indices:
            self.selected_indices.remove(idx)
            # 更新列表和显示
            self.update_selected_list()
            self.plot_molecule()
            
    def on_click(self, event):
        """处理鼠标点击事件"""
        if event.inaxes != self.canvas.axes:
            return
            
        if not event.button == 1:  # 只处理左键点击
            return
            
        # 获取点击位置
        x, y = event.xdata, event.ydata
        
        # 找到最近的原子
        positions = self.molecule.positions
        distances = np.sqrt((positions[:, 0] - x)**2 + (positions[:, 1] - y)**2)
        nearest_idx = np.argmin(distances)
        
        # 切换选中状态
        if nearest_idx in self.selected_indices:
            self.selected_indices.remove(nearest_idx)
        else:
            self.selected_indices.add(nearest_idx)
            
        self.update_selected_list()
        self.plot_molecule()
        
    def add_current_element(self):
        """添加当前选中元素的所有原子"""
        element = self.element_combo.currentText()
        symbols = self.molecule.get_chemical_symbols()
        
        for i, symbol in enumerate(symbols):
            if symbol == element:
                self.selected_indices.add(i)
                
        self.update_selected_list()
        self.plot_molecule()
        
    def remove_current_element(self):
        """删除当前选中元素的所有原子"""
        element = self.element_combo.currentText()
        symbols = self.molecule.get_chemical_symbols()
        
        for i, symbol in enumerate(symbols):
            if symbol == element:
                self.selected_indices.discard(i)
                
        self.update_selected_list()
        self.plot_molecule()
        
    def add_all_atoms(self):
        """添加所有原子"""
        self.selected_indices = set(range(len(self.molecule)))
        self.update_selected_list()
        self.plot_molecule()
        
    def clear_all_atoms(self):
        """清空所有选中的原子"""
        self.selected_indices.clear()
        self.update_selected_list()
        self.plot_molecule()
        
    def plot_molecule(self):
        """绘制分子结构，标记选中的原子"""
        self.canvas.figure.clear()
        ax = self.canvas.figure.add_subplot(111, projection='3d')
        
        # 绘制所有原子
        positions = self.molecule.positions
        symbols = self.molecule.get_chemical_symbols()
        
        # 绘制未选中的原子
        for i in range(len(self.molecule)):
            if i not in self.selected_indices:
                pos = positions[i]
                symbol = symbols[i]
                color = ELEMENT_COLORS.get(symbol, DEFAULT_COLOR)
                ax.scatter(pos[0], pos[1], pos[2], c=color, s=100)
        
        # 绘制选中的原子（红色）
        for i in self.selected_indices:
            pos = positions[i]
            ax.scatter(pos[0], pos[1], pos[2], c='red', s=100)
            
        # 绘制分子键
        self.draw_molecular_bonds(self.molecule, ax)
        
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        
        self.canvas.draw()
        
    def update_selected_list(self):
        """更新选中原子列表"""
        self.selected_list.clear()
        symbols = self.molecule.get_chemical_symbols()
        for idx in sorted(self.selected_indices):
            self.selected_list.addItem(f"{symbols[idx]}{idx}")
            
    def populate_available_elements(self):
        """填充可用元素下拉框"""
        symbols = set(self.molecule.get_chemical_symbols())
        self.element_combo.clear()
        self.element_combo.addItems(sorted(symbols))
        
    def get_export_options(self):
        """获取导出选项"""
        return {
            'export_adsorption': self.export_adsorption.isChecked(),
            'export_structures': self.export_structures.isChecked(),
            'export_format': self.format_combo.currentText(),
            'export_images': self.export_images.isChecked(),
            'export_json': self.export_json.isChecked()
        }
        
    def accept(self):
        """确认选择"""
        if not self.selected_indices:
            QMessageBox.warning(self, "警告", "请至少选择一个锚点原子")
            return
        super().accept()

    def draw_molecular_bonds(self, atoms, ax):
        """绘制分子键"""
        # 获取原子位置和共价半径
        positions = atoms.positions
        symbols = atoms.get_chemical_symbols()
        
        # 计算所有原子对之间的距离
        n_atoms = len(atoms)
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                # 计算两个原子之间的距离
                dist = np.linalg.norm(positions[i] - positions[j])
                
                # 获取两个原子的共价半径
                r1 = COVALENT_RADII.get(symbols[i], 1.0)
                r2 = COVALENT_RADII.get(symbols[j], 1.0)
                
                # 如果距离小于共价半径之和的1.3倍，则认为存在化学键
                if dist < 1.3 * (r1 + r2):
                    # 绘制键（用线段表示）
                    ax.plot([positions[i][0], positions[j][0]],
                           [positions[i][1], positions[j][1]],
                           [positions[i][2], positions[j][2]], 'k-', alpha=0.5) 

class BatchPathDialog(QDialog):
    """批量路径处理对话框"""
    
    def __init__(self, molecule, current_anchor, parent=None):
        super().__init__(parent)
        self.molecule = molecule
        self.current_anchor = current_anchor
        self.parent = parent
        self.path_items = []  # 存储路径项 [(json_path, anchor_index, interpolation_steps), ...]
        self.initUI()
        
    def initUI(self):
        """初始化UI"""
        self.setWindowTitle("批量路径处理")
        self.setMinimumSize(800, 600)
        
        # 创建主布局
        layout = QVBoxLayout()
        
        # 路径列表
        list_group = QGroupBox("路径列表")
        list_layout = QVBoxLayout()
        
        # 创建表格
        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["JSON文件", "锚点原子", "插值步数", "操作"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)
        list_layout.addWidget(self.table)
        
        # 添加路径按钮
        add_btn = QPushButton("添加路径点")
        add_btn.clicked.connect(self.add_path_item)
        list_layout.addWidget(add_btn)
        
        list_group.setLayout(list_layout)
        layout.addWidget(list_group)
        
        # 导出选项组
        export_group = QGroupBox("导出选项")
        export_layout = QVBoxLayout()
        
        # 导出吸附结构复选框
        self.export_adsorption = QCheckBox("导出吸附结构")
        self.export_adsorption.setChecked(True)
        export_layout.addWidget(self.export_adsorption)
        
        # 导出基底和分子结构复选框
        self.export_structures = QCheckBox("导出基底结构和分子结构")
        self.export_structures.setChecked(True)
        export_layout.addWidget(self.export_structures)
        
        # 导出格式选择
        format_layout = QHBoxLayout()
        format_layout.addWidget(QLabel("导出格式:"))
        self.format_combo = QComboBox()
        self.format_combo.addItems(["vasp_cartesian", "vasp_fractional", "xsf", "xyz", "cif"])
        format_layout.addWidget(self.format_combo)
        export_layout.addLayout(format_layout)
        
        # 导出图像复选框
        self.export_images = QCheckBox("导出图像")
        self.export_images.setChecked(True)
        export_layout.addWidget(self.export_images)
        
        # 导出JSON复选框
        self.export_json = QCheckBox("导出当前状态(JSON)")
        self.export_json.setChecked(True)
        export_layout.addWidget(self.export_json)
        
        export_group.setLayout(export_layout)
        layout.addWidget(export_group)
        
        # 确定和取消按钮
        button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        )
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
        
        self.setLayout(layout)
        
    def add_path_item(self):
        """添加路径项"""
        # 打开文件对话框选择json文件
        json_file, _ = QFileDialog.getOpenFileName(
            self, "选择状态文件", "", "JSON文件 (*.json)"
        )
        
        if json_file:
            # 添加新行
            row = self.table.rowCount()
            self.table.insertRow(row)
            
            # 文件路径
            self.table.setItem(row, 0, QTableWidgetItem(json_file))
            
            # 锚点原子选择
            anchor_combo = QComboBox()
            symbols = self.molecule.get_chemical_symbols()
            for i, symbol in enumerate(symbols):
                anchor_combo.addItem(f"{symbol}{i}")
            anchor_combo.setCurrentIndex(self.current_anchor)  # 设置默认值为当前锚点
            self.table.setCellWidget(row, 1, anchor_combo)
            
            # 插值步数
            if row > 0:  # 只有第二个及以后的项才需要插值步数
                steps_spin = QSpinBox()
                steps_spin.setRange(1, 20)
                steps_spin.setValue(5)  # 默认5步
                self.table.setCellWidget(row, 2, steps_spin)
            
            # 删除按钮
            delete_btn = QPushButton("删除")
            delete_btn.clicked.connect(lambda: self.remove_path_item(row))
            self.table.setCellWidget(row, 3, delete_btn)
            
            # 更新所有行的插值步数显示
            self.update_interpolation_steps()
    
    def remove_path_item(self, row):
        """删除路径项"""
        self.table.removeRow(row)
        self.update_interpolation_steps()
    
    def update_interpolation_steps(self):
        """更新插值步数显示"""
        for row in range(self.table.rowCount()):
            if row == 0:  # 第一个项不需要插值步数
                self.table.setCellWidget(row, 2, None)
            else:
                # 如果还没有插值步数控件，添加一个
                if not isinstance(self.table.cellWidget(row, 2), QSpinBox):
                    steps_spin = QSpinBox()
                    steps_spin.setRange(1, 20)
                    steps_spin.setValue(5)  # 默认5步
                    self.table.setCellWidget(row, 2, steps_spin)
    
    def get_path_data(self):
        """获取路径数据"""
        path_data = []
        for row in range(self.table.rowCount()):
            json_file = self.table.item(row, 0).text()
            anchor_combo = self.table.cellWidget(row, 1)
            anchor_index = anchor_combo.currentIndex()
            
            if row == 0:
                steps = 0  # 第一个点不需要插值步数
            else:
                steps_spin = self.table.cellWidget(row, 2)
                steps = steps_spin.value()
            
            path_data.append((json_file, anchor_index, steps))
        
        return path_data
    
    def get_export_options(self):
        """获取导出选项"""
        return {
            'export_adsorption': self.export_adsorption.isChecked(),
            'export_structures': self.export_structures.isChecked(),
            'export_format': self.format_combo.currentText(),
            'export_images': self.export_images.isChecked(),
            'export_json': self.export_json.isChecked()
        }
        
    def accept(self):
        """确认选择"""
        if self.table.rowCount() < 2:
            QMessageBox.warning(self, "警告", "请至少添加两个路径点")
            return
        super().accept() 