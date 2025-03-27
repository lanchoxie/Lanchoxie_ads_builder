#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
将分子旋转、平移到基底(substrate)的指定位置，并检查分子是否在基底晶胞内
"""

import os
import numpy as np
from ase import Atoms
from ase.io import read, write
from scipy.spatial.distance import cdist
from scipy.spatial import KDTree
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from ase.data import covalent_radii, atomic_numbers

def find_surface_atoms(substrate, tolerance=0.5, surface_layers=1):
    """
    寻找基底表面原子
    
    参数:
        substrate (ase.Atoms): 基底结构
        tolerance (float): 容许误差(Å)
        surface_layers (int): 要考虑的表面层数
        
    返回:
        list: 表面原子的索引列表
    """
    # 获取原子坐标
    positions = substrate.get_positions()
    
    # 找出z方向的最大值
    max_z = np.max(positions[:, 2])
    
    if surface_layers == 1:
        # 传统方法：只选择最上层原子
        surface_indices = np.where(positions[:, 2] > max_z - tolerance)[0]
    else:
        # 多层方法
        # 1. 先找到最上层原子
        top_layer_indices = np.where(positions[:, 2] > max_z - tolerance)[0]
        
        # 2. 获取所有原子的z坐标并排序
        z_coords = np.unique(positions[:, 2])
        z_coords = np.sort(z_coords)[::-1]  # 降序排列
        
        # 3. 计算层间距（估计值）
        if len(z_coords) >= 2:
            layer_distance = np.mean(np.diff(z_coords[:5])) if len(z_coords) >= 5 else np.mean(np.diff(z_coords))
        else:
            layer_distance = tolerance  # 如果只有一层，使用tolerance作为默认层间距
        
        # 4. 选择多层原子
        surface_indices = []
        for layer in range(surface_layers):
            layer_min_z = max_z - (layer * layer_distance) - tolerance
            layer_max_z = max_z - (layer * layer_distance) + tolerance if layer > 0 else float('inf')
            layer_indices = np.where((positions[:, 2] > layer_min_z) & (positions[:, 2] < layer_max_z))[0]
            surface_indices.extend(layer_indices)
        
        # 去重并排序
        surface_indices = np.unique(surface_indices)
    
    return surface_indices.tolist()

def find_atom_by_element(atoms, element, surface_only=True, tolerance=0.5):
    """
    寻找指定元素的原子
    
    参数:
        atoms (ase.Atoms): 原子结构
        element (str): 元素符号
        surface_only (bool): 是否仅在表面搜索
        tolerance (float): 容许误差(Å)，用于表面原子判断
        
    返回:
        list: 匹配元素的原子索引列表
    """
    # 获取所有原子的元素符号
    symbols = atoms.get_chemical_symbols()
    
    # 寻找匹配元素的原子索引
    indices = [i for i, sym in enumerate(symbols) if sym == element]
    
    # 如果需要仅在表面搜索
    if surface_only and indices:
        surface_indices = find_surface_atoms(atoms, tolerance)
        indices = [i for i in indices if i in surface_indices]
    
    return indices

def is_inside_cell_xy(position, cell):
    """
    检查一个点在x-y平面上的投影是否在晶胞x-y平面内
    
    参数:
        position (numpy.ndarray): 点的坐标 [x, y, z]
        cell (numpy.ndarray): 晶胞矩阵 (3x3)
        
    返回:
        bool: 如果点的x-y投影在晶胞x-y平面内则返回True
    """
    # 获取晶胞的二维投影
    cell_2d = cell[:2, :2]
    
    # 计算晶胞在x-y平面上的四个顶点
    vertices = np.array([
        [0, 0],
        cell_2d[0],
        cell_2d[0] + cell_2d[1],
        cell_2d[1]
    ])
    
    # 获取点的x-y投影
    point_2d = position[:2]
    
    # 检查点是否在晶胞的x-y投影内部
    # 可以通过检查点相对于每条边的位置来判断
    n_vertices = len(vertices)
    inside = True
    
    for i in range(n_vertices):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % n_vertices]
        
        # 检查点是否在边的"右侧"
        # 使用叉积判断：(x2-x1)(y-y1) - (y2-y1)(x-x1) > 0 表示点在边的右侧
        if ((x2 - x1) * (point_2d[1] - y1) - (y2 - y1) * (point_2d[0] - x1)) < 0:
            inside = False
            break
    
    return inside

def is_molecule_inside_cell_xy(molecule, cell):
    """
    检查分子在x-y平面上的投影是否完全在晶胞x-y平面内
    
    参数:
        molecule (ase.Atoms): 分子结构
        cell (numpy.ndarray): 晶胞矩阵
        
    返回:
        tuple: (是否在晶胞内, 不在晶胞内的原子索引列表)
    """
    # 检查每个原子是否在晶胞x-y投影内
    atoms_outside = []
    for i, pos in enumerate(molecule.get_positions()):
        if not is_inside_cell_xy(pos, cell):
            atoms_outside.append(i)
    
    return len(atoms_outside) == 0, atoms_outside

def get_xy_convex_hull(substrate):
    """
    获取基底在x-y平面的凸包
    
    参数:
        substrate (ase.Atoms): 基底结构
        
    返回:
        ConvexHull: 基底在x-y平面的凸包
    """
    # 获取基底表面原子
    surface_indices = find_surface_atoms(substrate)
    surface_positions = substrate.positions[surface_indices]
    
    # 提取x-y坐标
    xy_positions = surface_positions[:, :2]
    
    # 计算凸包
    hull = ConvexHull(xy_positions)
    
    return hull

def is_point_in_convex_hull(point_2d, hull):
    """
    检查点是否在凸包内
    
    参数:
        point_2d (numpy.ndarray): 点的二维坐标 [x, y]
        hull (ConvexHull): 凸包
        
    返回:
        bool: 如果点在凸包内则返回True
    """
    from scipy.spatial import Delaunay
    
    # 获取凸包顶点
    hull_points = hull.points[hull.vertices]
    
    # 创建Delaunay三角剖分
    tri = Delaunay(hull_points)
    
    # 检查点是否在三角剖分内
    return tri.find_simplex(point_2d) >= 0

def is_molecule_inside_substrate_xy(molecule, substrate):
    """
    检查分子在x-y平面上的投影是否完全在基底表面凸包内
    
    参数:
        molecule (ase.Atoms): 分子结构
        substrate (ase.Atoms): 基底结构
        
    返回:
        tuple: (是否在基底内, 不在基底内的原子索引列表)
    """
    try:
        # 获取基底表面凸包
        hull = get_xy_convex_hull(substrate)
        
        # 检查每个原子是否在凸包内
        atoms_outside = []
        for i, pos in enumerate(molecule.get_positions()):
            if not is_point_in_convex_hull(pos[:2], hull):
                atoms_outside.append(i)
        
        return len(atoms_outside) == 0, atoms_outside
    except Exception as e:
        print(f"检查分子是否在基底内时出错: {str(e)}")
        print("回退到简单的晶胞检查...")
        return is_molecule_inside_cell_xy(molecule, substrate.get_cell())

def print_atoms_outside_info(molecule, atoms_outside):
    """
    打印在晶胞/基底外的原子信息
    
    参数:
        molecule (ase.Atoms): 分子结构
        atoms_outside (list): 在晶胞/基底外的原子索引列表
    """
    if not atoms_outside:
        print("所有原子都在基底内")
        return
    
    print(f"有 {len(atoms_outside)} 个原子在基底外:")
    for i in atoms_outside:
        symbol = molecule.get_chemical_symbols()[i]
        pos = molecule.get_positions()[i]
        print(f"  原子 {i}: {symbol} 位置: {pos}")

def check_and_adjust_position(molecule, substrate):
    """
    检查分子是否在基底内，如果不在则尝试调整
    
    参数:
        molecule (ase.Atoms): 分子结构
        substrate (ase.Atoms): 基底结构
        
    返回:
        tuple: (是否需要调整, 调整后的分子)
    """
    # 检查分子是否在基底内
    in_substrate, atoms_outside = is_molecule_inside_substrate_xy(molecule, substrate)
    
    # 如果有原子在基底外，打印信息
    if not in_substrate:
        print_atoms_outside_info(molecule, atoms_outside)
        
        # 尝试调整分子位置
        print("尝试将分子调整至基底内...")
        
        # 获取基底细胞
        cell = substrate.get_cell()
        
        # 找到基底表面原子
        surface_indices = find_surface_atoms(substrate)
        surface_positions = substrate.positions[surface_indices]
        
        # 计算基底表面中心
        surface_center = np.mean(surface_positions, axis=0)
        
        # 计算分子中心
        mol_center = np.mean(molecule.get_positions(), axis=0)
        
        # 将分子在x-y平面上移动到基底表面中心
        # 修复广播错误 - 正确计算x-y平面上的调整向量
        adjustment_vector = np.zeros_like(molecule.positions)
        adjustment_vector[:, 0] = surface_center[0] - mol_center[0]  # X轴调整
        adjustment_vector[:, 1] = surface_center[1] - mol_center[1]  # Y轴调整
        # Z轴保持不变
        
        # 应用调整
        molecule.positions += adjustment_vector
        
        # 再次检查
        in_substrate, atoms_outside = is_molecule_inside_substrate_xy(molecule, substrate)
        
        if in_substrate:
            print("已成功将分子调整至基底内")
        else:
            print("调整后仍有原子在基底外，可能需要手动调整分子位置")
            print_atoms_outside_info(molecule, atoms_outside)
    
    # 返回是否需要调整和调整后的分子
    return (not in_substrate), molecule

def rotate_molecule(molecule, angles=(0, 0, 0), anchor_atom=None):
    """
    旋转分子
    
    参数:
        molecule (ase.Atoms): 分子结构
        angles (tuple): 绕x, y, z轴的旋转角度(度)
        anchor_atom (int, optional): 指定锚点原子索引，如果提供，则绕此原子旋转，否则绕分子中心旋转
        
    返回:
        ase.Atoms: 旋转后的分子
    """
    # 创建副本避免修改原始对象
    rotated_mol = molecule.copy()
    
    # 确定旋转中心
    if anchor_atom is not None and anchor_atom < len(rotated_mol):
        # 使用指定的锚点原子作为旋转中心
        center = rotated_mol[anchor_atom].position
    else:
        # 使用分子的几何中心
        center = np.mean(rotated_mol.get_positions(), axis=0)
    
    # 将分子平移到旋转中心为原点
    rotated_mol.positions -= center
    
    # 依次绕各轴旋转
    for i, angle in enumerate(angles):
        if angle != 0:
            # 构建旋转角度向量
            rot_angles = np.zeros(3)
            rot_angles[i] = np.radians(angle)
            
            # 旋转分子
            rotated_mol.rotate(rot_angles[i], 'xyz'[i])
    
    # 将分子平移回原位置
    rotated_mol.positions += center
    
    return rotated_mol

def position_molecule_on_site(molecule, substrate, site_position, height=2.0,
                             rotation=(0, 0, 0), anchor_atom=0):
    """
    将分子定位到基底的指定位置
    
    参数:
        molecule (ase.Atoms): 分子结构
        substrate (ase.Atoms): 基底结构
        site_position (numpy.ndarray): 吸附位置的3D坐标
        height (float): 分子底部距离表面的高度(Å)
        rotation (tuple): 旋转角度
        anchor_atom (int): 分子中的锚定原子索引，该原子将位于吸附位置上方
        
    返回:
        ase.Atoms: 定位后的分子
    """
    # 创建分子的副本
    mol = molecule.copy()
    
    # 旋转分子
    if rotation != (0, 0, 0):
        mol = rotate_molecule(mol, rotation)
    
    # 获取锚定原子位置
    anchor_pos = mol.positions[anchor_atom]
    
    # 分子其他原子位置相对于锚定原子的偏移
    rel_positions = mol.positions - anchor_pos
    
    # 计算新的锚定原子位置
    new_anchor_pos = np.array([site_position[0], site_position[1], site_position[2] + height])
    
    # 平移所有原子
    mol.positions = rel_positions + new_anchor_pos
    
    return mol

def adjust_cell_height(system, min_vacuum_above=5.0):
    """
    调整晶胞高度，确保分子顶部与晶胞顶部有足够的空间
    
    参数:
        system (ase.Atoms): 原子系统
        min_vacuum_above (float): 分子顶部与晶胞顶部的最小距离(Å)
        
    返回:
        ase.Atoms: 调整后的系统
    """
    # 获取原子坐标
    positions = system.get_positions()
    
    # 获取晶胞尺寸
    cell = system.get_cell()
    
    # 获取z方向的最大原子坐标（分子顶部）
    max_z = np.max(positions[:, 2])
    
    # 获取当前晶胞在z方向的长度
    cell_z = cell[2, 2]
    
    # 计算分子顶部到晶胞顶部的距离
    distance_to_top = cell_z - max_z
    
    # 如果距离小于最小要求，则调整晶胞高度
    if distance_to_top < min_vacuum_above:
        # 计算需要增加的高度
        additional_height = min_vacuum_above - distance_to_top
        
        # 调整晶胞
        new_cell = cell.copy()
        new_cell[2, 2] += additional_height
        system.set_cell(new_cell)
        
        print(f"晶胞高度已调整: +{additional_height:.2f}Å (新高度: {new_cell[2, 2]:.2f}Å)")
    
    return system

def find_nearest_surface_atom(substrate, element=None, position=None):
    """
    寻找最接近指定位置的表面原子，或特定元素的表面原子
    
    参数:
        substrate (ase.Atoms): 基底结构
        element (str): 元素符号，如不指定则搜索所有元素
        position (numpy.ndarray): 3D坐标，寻找最接近该位置的原子
        
    返回:
        tuple: (原子索引, 原子位置)
    """
    # 寻找表面原子
    surface_indices = find_surface_atoms(substrate)
    
    # 如果指定了元素，过滤表面原子
    if element is not None:
        symbols = substrate.get_chemical_symbols()
        surface_indices = [i for i in surface_indices if symbols[i] == element]
    
    if not surface_indices:
        raise ValueError(f"未找到满足条件的表面原子: element={element}")
    
    # 获取表面原子位置
    positions = substrate.positions[surface_indices]
    
    # 如果指定了位置，寻找最接近的原子
    if position is not None:
        distances = np.linalg.norm(positions - position, axis=1)
        nearest_idx = np.argmin(distances)
        atom_idx = surface_indices[nearest_idx]
        return atom_idx, substrate.positions[atom_idx]
    else:
        # 否则返回第一个匹配的原子
        atom_idx = surface_indices[0]
        return atom_idx, substrate.positions[atom_idx]

def find_nearest_neighbor_pairs(positions, max_distance=3.5):
    """
    寻找最近邻原子对
    
    参数:
        positions (numpy.ndarray): 原子坐标
        max_distance (float): 最大考虑距离(Å)
        
    返回:
        list: 最近邻原子对索引列表 [(i, j), ...]
    """
    # 计算所有原子间的距离矩阵
    distances = cdist(positions, positions)
    
    # 将自身距离设为无穷大
    np.fill_diagonal(distances, np.inf)
    
    # 寻找所有距离小于max_distance的原子对
    pairs = []
    for i in range(len(positions)):
        for j in range(i+1, len(positions)):
            if distances[i, j] < max_distance:
                pairs.append((i, j))
    
    return pairs

def find_bridge_sites(substrate, max_distance=3.5, surface_layers=1):
    """
    寻找基底表面的桥位(两个原子之间的中点)
    
    参数:
        substrate (ase.Atoms): 基底结构
        max_distance (float): 最大考虑的原子间距(Å)
        surface_layers (int): 要考虑的表面层数
        
    返回:
        list: 桥位坐标列表 [(x, y, z), ...]
        list: 桥位元素列表 [['元素1', '元素2'], ...]
    """
    # 寻找表面原子
    surface_indices = find_surface_atoms(substrate, surface_layers=surface_layers)
    surface_positions = substrate.positions[surface_indices]
    
    # 获取原子元素符号
    symbols = substrate.get_chemical_symbols()
    surface_symbols = [symbols[i] for i in surface_indices]
    
    # 寻找最近邻原子对
    pairs = find_nearest_neighbor_pairs(surface_positions, max_distance)
    
    # 计算桥位坐标
    bridge_sites = []
    bridge_elements = []
    for i, j in pairs:
        pos1 = surface_positions[i]
        pos2 = surface_positions[j]
        bridge_pos = (pos1 + pos2) / 2
        bridge_sites.append(bridge_pos)
        
        # 记录形成桥位的两个原子的元素
        element1 = surface_symbols[i]
        element2 = surface_symbols[j]
        bridge_elements.append([element1, element2])
    
    return bridge_sites, bridge_elements

def are_atoms_in_plane(pos1, pos2, pos3, tolerance=0.1):
    """
    判断三个原子是否近似在同一平面
    
    参数:
        pos1, pos2, pos3 (numpy.ndarray): 三个原子的三维坐标
        tolerance (float): 容许误差(Å)
        
    返回:
        bool: 如果三个原子在同一平面则返回True
    """
    # 计算三个向量
    v1 = pos2 - pos1
    v2 = pos3 - pos1
    
    # 计算法向量
    normal = np.cross(v1, v2)
    
    # 检查法向量的z分量是否显著
    return abs(normal[2]) > tolerance

def find_hollow_sites(substrate, max_distance=5.0, surface_layers=1):
    """
    寻找基底表面的空位(三个原子围成的三角形中心)
    
    参数:
        substrate (ase.Atoms): 基底结构
        max_distance (float): 最大考虑的原子间距(Å)
        surface_layers (int): 要考虑的表面层数
        
    返回:
        list: 空位坐标列表 [(x, y, z), ...]
        list: 空位元素列表 [['元素1', '元素2', '元素3'], ...]
    """
    # 寻找表面原子
    surface_indices = find_surface_atoms(substrate, surface_layers=surface_layers)
    surface_positions = substrate.positions[surface_indices]
    
    # 获取原子元素符号
    symbols = substrate.get_chemical_symbols()
    surface_symbols = [symbols[i] for i in surface_indices]
    
    # 使用KDTree找到每个表面原子的近邻
    tree = KDTree(surface_positions)
    
    # 寻找所有可能的三角形
    hollow_sites = []
    hollow_elements = []
    
    for i in range(len(surface_positions)):
        # 寻找距离小于max_distance的邻居
        dists, neighbors = tree.query(surface_positions[i], k=7, distance_upper_bound=max_distance)
        
        # 移除自身和无效邻居
        valid_neighbors = [neighbors[j] for j in range(1, len(neighbors)) if neighbors[j] < len(surface_positions)]
        
        # 对于每对邻居，形成三角形
        for j in range(len(valid_neighbors)):
            for k in range(j+1, len(valid_neighbors)):
                # 三个原子的索引
                idx1 = i
                idx2 = valid_neighbors[j]
                idx3 = valid_neighbors[k]
                
                # 三个原子的坐标
                pos1 = surface_positions[idx1]
                pos2 = surface_positions[idx2]
                pos3 = surface_positions[idx3]
                
                # 检查三个原子是否形成一个合理的三角形
                side1 = np.linalg.norm(pos2 - pos1)
                side2 = np.linalg.norm(pos3 - pos1)
                side3 = np.linalg.norm(pos3 - pos2)
                
                # 检查三角形边长是否相近(判断是否为近似等边三角形)
                avg_side = (side1 + side2 + side3) / 3
                if (abs(side1 - avg_side) < 1.0 and
                    abs(side2 - avg_side) < 1.0 and
                    abs(side3 - avg_side) < 1.0):
                    
                    # 计算三角形的质心(重心)
                    hollow_pos = (pos1 + pos2 + pos3) / 3
                    
                    # 将空位高度调整为三个原子的平均高度
                    hollow_pos[2] = (pos1[2] + pos2[2] + pos3[2]) / 3
                    
                    hollow_sites.append(hollow_pos)
                    
                    # 记录形成空位的三个原子的元素
                    element1 = surface_symbols[idx1]
                    element2 = surface_symbols[idx2]
                    element3 = surface_symbols[idx3]
                    hollow_elements.append([element1, element2, element3])
    
    return hollow_sites, hollow_elements

def get_best_bridge_site(substrate, element=None, max_distance=3.5):
    """
    获取最佳桥位
    
    参数:
        substrate (ase.Atoms): 基底结构
        element (str): 如果指定，则寻找包含此元素的桥位
        max_distance (float): 最大考虑的原子间距(Å)
        
    返回:
        numpy.ndarray: 桥位坐标
    """
    # 寻找所有桥位
    bridge_sites, bridge_elements = find_bridge_sites(substrate, max_distance)
    
    if not bridge_sites:
        raise ValueError("未找到合适的桥位")
    
    # 如果指定了元素，筛选包含此元素的桥位
    if element is not None:
        filtered_sites = []
        filtered_elements = []
        
        for i, (site, elements) in enumerate(zip(bridge_sites, bridge_elements)):
            if element in elements:
                filtered_sites.append(site)
                filtered_elements.append(elements)
        
        if filtered_sites:
            bridge_sites = filtered_sites
            bridge_elements = filtered_elements
        else:
            print(f"警告：未找到包含元素 {element} 的桥位，使用所有桥位")
    
    # 选择桥位，优先选择z坐标最高的
    if bridge_sites:
        # 找到z坐标最高的桥位
        max_z_idx = np.argmax([site[2] for site in bridge_sites])
        return bridge_sites[max_z_idx]
    
    # 如果没有找到合适的桥位，返回第一个
    return bridge_sites[0]

def get_best_hollow_site(substrate, element=None, max_distance=5.0):
    """
    获取最佳空位
    
    参数:
        substrate (ase.Atoms): 基底结构
        element (str): 如果指定，则寻找包含此元素的空位
        max_distance (float): 最大考虑的原子间距(Å)
        
    返回:
        numpy.ndarray: 空位坐标
    """
    # 寻找所有空位
    hollow_sites, hollow_elements = find_hollow_sites(substrate, max_distance)
    
    if not hollow_sites:
        raise ValueError("未找到合适的空位")
    
    # 如果指定了元素，筛选包含此元素的空位
    if element is not None:
        filtered_sites = []
        filtered_elements = []
        
        for i, (site, elements) in enumerate(zip(hollow_sites, hollow_elements)):
            if element in elements:
                filtered_sites.append(site)
                filtered_elements.append(elements)
        
        if filtered_sites:
            hollow_sites = filtered_sites
            hollow_elements = filtered_elements
        else:
            print(f"警告：未找到包含元素 {element} 的空位，使用所有空位")
    
    # 选择空位，优先选择z坐标最高的
    if hollow_sites:
        # 找到z坐标最高的空位
        max_z_idx = np.argmax([site[2] for site in hollow_sites])
        return hollow_sites[max_z_idx]
    
    # 如果没有找到合适的空位，返回第一个
    return hollow_sites[0]

def get_atomic_radius(element):
    """
    获取元素的共价半径
    
    参数:
        element (str): 元素符号
        
    返回:
        float: 元素的共价半径
    """
    return covalent_radii[atomic_numbers[element]]

def check_atom_overlap(system, scale_factor=0.7):
    """
    检查系统中的原子是否存在重叠
    
    参数:
        system (ase.Atoms): 原子系统
        scale_factor (float): 缩放因子，用于调整重叠判断标准。
                             值越小，容忍的重叠程度越大。
        
    返回:
        list: 重叠原子对的列表 [(i, j, distance, min_distance), ...]
    """
    # 获取原子位置
    positions = system.get_positions()
    
    # 获取原子共价半径
    symbols = system.get_chemical_symbols()
    radii = np.array([get_atomic_radius(symbol) for symbol in symbols])
    
    # 使用KDTree加速距离计算
    tree = KDTree(positions)
    
    # 寻找重叠原子
    overlaps = []
    
    for i, (pos, radius) in enumerate(zip(positions, radii)):
        # 寻找可能与当前原子重叠的所有原子
        # 搜索半径为两个原子半径之和乘以缩放因子
        max_dist = radius * 2 * scale_factor
        indices = tree.query_ball_point(pos, max_dist)
        
        # 移除自身
        if i in indices:
            indices.remove(i)
        
        for j in indices:
            # 避免重复计算 (i,j) 和 (j,i)
            if j < i:
                continue
                
            # 计算两个原子之间的距离
            dist = np.linalg.norm(positions[i] - positions[j])
            
            # 计算两个原子允许的最小距离
            min_dist = (radii[i] + radii[j]) * scale_factor
            
            # 如果距离小于最小距离，则认为原子重叠
            if dist < min_dist:
                overlaps.append((i, j, dist, min_dist))
    
    return overlaps

def identify_molecule_substrate_atoms(system, n_substrate_atoms=None):
    """
    将系统中的原子分为基底原子和分子原子
    
    参数:
        system (ase.Atoms): 原子系统
        n_substrate_atoms (int): 基底原子的数量，如不指定则尝试自动识别
        
    返回:
        tuple: (基底原子索引列表, 分子原子索引列表)
    """
    # 如果未指定基底原子数量，尝试自动识别
    if n_substrate_atoms is None:
        # 一种简单的策略是寻找z坐标分布中的明显间隙
        z_coords = system.get_positions()[:, 2]
        z_sorted = np.sort(z_coords)
        z_diff = np.diff(z_sorted)
        
        # 查找最大间隙
        max_gap_idx = np.argmax(z_diff)
        n_substrate_atoms = max_gap_idx + 1
        
        print(f"自动识别的基底原子数量: {n_substrate_atoms}")
    
    substrate_indices = list(range(n_substrate_atoms))
    molecule_indices = list(range(n_substrate_atoms, len(system)))
    
    return substrate_indices, molecule_indices

def validate_system(system, n_substrate_atoms=None, scale_factor=0.7, min_mol_surf_dist=1.5, verbose=True):
    """
    验证吸附系统的合理性
    
    参数:
        system (ase.Atoms): 原子系统
        n_substrate_atoms (int): 基底原子的数量，如不指定则尝试自动识别
        scale_factor (float): 重叠判断的缩放因子
        min_mol_surf_dist (float): 分子到表面的最小距离
        verbose (bool): 是否输出详细信息
        
    返回:
        tuple: (是否合理, 问题列表)
    """
    issues = []
    
    # 识别基底和分子原子
    substrate_indices, molecule_indices = identify_molecule_substrate_atoms(system, n_substrate_atoms)
    
    if not molecule_indices:
        issues.append("未能识别到分子原子")
        return False, issues
    
    # 获取原子位置和符号
    positions = system.get_positions()
    symbols = system.get_chemical_symbols()
    
    # 检查1: 原子重叠
    if verbose:
        print("检查原子重叠...")
    overlaps = check_atom_overlap(system, scale_factor)
    
    # 按原子类型分类重叠
    mol_mol_overlaps = []  # 分子内部重叠
    sub_sub_overlaps = []  # 基底内部重叠
    mol_sub_overlaps = []  # 分子-基底重叠
    
    for i, j, dist, min_dist in overlaps:
        if i in molecule_indices and j in molecule_indices:
            mol_mol_overlaps.append((i, j, dist, min_dist))
        elif i in substrate_indices and j in substrate_indices:
            sub_sub_overlaps.append((i, j, dist, min_dist))
        else:
            mol_sub_overlaps.append((i, j, dist, min_dist))
    
    # 记录重叠信息
    if mol_mol_overlaps:
        issue_text = f"发现 {len(mol_mol_overlaps)} 个分子内部原子重叠"
        issues.append(issue_text)
        if verbose:
            print(f"警告: {issue_text}")
            for i, j, dist, min_dist in mol_mol_overlaps:
                print(f"  原子 {i}({symbols[i]}) 和 {j}({symbols[j]}) 距离 {dist:.3f}Å < {min_dist:.3f}Å")
    
    if mol_sub_overlaps:
        issue_text = f"发现 {len(mol_sub_overlaps)} 个分子-基底原子重叠"
        issues.append(issue_text)
        if verbose:
            print(f"警告: {issue_text}")
            for i, j, dist, min_dist in mol_sub_overlaps:
                print(f"  原子 {i}({symbols[i]}) 和 {j}({symbols[j]}) 距离 {dist:.3f}Å < {min_dist:.3f}Å")
    
    # 检查2: 分子到表面的距离
    if verbose:
        print("检查分子到表面的距离...")
    
    # 获取基底表面原子
    substrate_positions = positions[substrate_indices]
    max_z_substrate = np.max(substrate_positions[:, 2])
    
    # 获取分子底部原子
    molecule_positions = positions[molecule_indices]
    min_z_molecule = np.min(molecule_positions[:, 2])
    
    # 计算分子到表面的距离
    mol_surf_distance = min_z_molecule - max_z_substrate
    
    if verbose:
        print(f"分子到表面的距离: {mol_surf_distance:.3f}Å")
    
    if mol_surf_distance < min_mol_surf_dist:
        issue_text = f"分子到表面的距离 ({mol_surf_distance:.3f}Å) 小于最小要求 ({min_mol_surf_dist:.3f}Å)"
        issues.append(issue_text)
        if verbose:
            print(f"警告: {issue_text}")
    
    # 结构合理性判断
    is_valid = (len(mol_sub_overlaps) == 0) and (mol_surf_distance >= min_mol_surf_dist)
    
    if verbose:
        if is_valid:
            print("结构验证通过: 吸附结构合理")
        else:
            print("结构验证失败: 吸附结构存在问题")
    
    return is_valid, issues

def create_multiple_site_options(substrate, element=None, site_type=None, max_sites=5):
    """
    创建多个可能的吸附位点选项
    
    参数:
        substrate (ase.Atoms): 基底结构
        element (str): 目标元素
        site_type (str): 吸附位点类型 (top, bridge, hollow)
        max_sites (int): 最大位点数量
        
    返回:
        list: 吸附位点列表
    """
    sites = []
    
    try:
        if site_type == 'top':
            # 对于顶位，找到所有表面原子
            surface_indices = find_surface_atoms(substrate)
            
            # 如果指定了元素，筛选指定元素的表面原子
            if element is not None:
                symbols = substrate.get_chemical_symbols()
                element_indices = [i for i in surface_indices if symbols[i] == element]
                if element_indices:
                    surface_indices = element_indices
            
            # 按z坐标降序排序，获取最顶部的原子
            sorted_indices = sorted(surface_indices, 
                                   key=lambda i: substrate.positions[i][2], 
                                   reverse=True)
            
            # 取前max_sites个原子位置作为吸附位点
            sites = [substrate.positions[i] for i in sorted_indices[:max_sites]]
            
        elif site_type == 'bridge':
            # 对于桥位，找到所有可能的桥位
            bridge_sites, bridge_elements = find_bridge_sites(substrate)
            
            # 按z坐标降序排序
            sorted_indices = sorted(range(len(bridge_sites)), 
                                   key=lambda i: bridge_sites[i][2], 
                                   reverse=True)
            
            # 取前max_sites个桥位
            sites = [bridge_sites[i] for i in sorted_indices[:min(max_sites, len(bridge_sites))]]
            
        elif site_type == 'hollow':
            # 对于空位，找到所有可能的空位
            hollow_sites, hollow_elements = find_hollow_sites(substrate)
            
            # 按z坐标降序排序
            sorted_indices = sorted(range(len(hollow_sites)), 
                                   key=lambda i: hollow_sites[i][2], 
                                   reverse=True)
            
            # 取前max_sites个空位
            sites = [hollow_sites[i] for i in sorted_indices[:min(max_sites, len(hollow_sites))]]
            
        else:
            # 如果未指定位点类型，尝试所有类型
            top_sites = create_multiple_site_options(substrate, element, 'top', max(1, max_sites // 3))
            bridge_sites = create_multiple_site_options(substrate, element, 'bridge', max(1, max_sites // 3))
            hollow_sites = create_multiple_site_options(substrate, element, 'hollow', max(1, max_sites // 3))
            
            sites = top_sites + bridge_sites + hollow_sites
            
    except Exception as e:
        print(f"创建多个吸附位点选项时出错: {str(e)}")
        # 如果出错，至少返回一个表面中心点作为候选
        surface_indices = find_surface_atoms(substrate)
        surface_center = np.mean(substrate.positions[surface_indices], axis=0)
        sites = [surface_center]
    
    # 确保返回的位点数量不超过max_sites
    return sites[:max_sites]

def create_adsorption_system(molecule_file, substrate_file, site_spec=None, 
                            height=2.0, rotation=(0, 0, 0), anchor_atom=0,
                            min_vacuum_above=10.0, output_file=None,
                            check_position=True, validate=True, scale_factor=0.7,
                            min_mol_surf_dist=1.5, max_trial_sites=5):
    """
    创建分子吸附系统
    
    参数:
        molecule_file (str): 分子文件路径
        substrate_file (str): 基底文件路径
        site_spec (dict): 吸附位置规格，可以是:
                         {'element': 'Pt'} - 选择特定元素的表面原子
                         {'position': [x, y, z]} - 选择最接近指定位置的表面原子
                         {'element': 'Pt', 'position': [x, y, z]} - 同时限制元素和位置
                         {'site_type': 'top/bridge/hollow'} - 指定位点类型
        height (float): 分子底部距离表面的高度(Å)
        rotation (tuple): 旋转角度
        anchor_atom (int): 分子中的锚定原子索引
        min_vacuum_above (float): 分子顶部与晶胞顶部的最小距离(Å)
        output_file (str): 输出文件路径
        check_position (bool): 是否检查并调整分子位置，使其在基底内
        validate (bool): 是否验证吸附结构的合理性
        scale_factor (float): 原子重叠判断的缩放因子
        min_mol_surf_dist (float): 分子到表面的最小距离
        max_trial_sites (int): 尝试的最大位点数量
        
    返回:
        tuple: (吸附系统, 是否合理, 问题列表)
    """
    # 读取结构
    molecule = read(molecule_file)
    substrate = read(substrate_file)
    
    site_type = None
    if site_spec and 'site_type' in site_spec:
        site_type = site_spec['site_type']
    
    element = None
    if site_spec and 'element' in site_spec:
        element = site_spec['element']
    
    # 如果指定了精确位置，直接尝试该位置
    if site_spec and 'position' in site_spec:
        site_positions = [np.array(site_spec['position'])]
    else:
        # 否则生成多个位点选项
        site_positions = create_multiple_site_options(substrate, element, site_type, max_trial_sites)
    
    best_system = None
    best_issues = ["未能找到合适的吸附位置"]
    is_valid = False
    
    # 尝试所有候选位点
    print(f"尝试 {len(site_positions)} 个候选吸附位点...")
    
    for i, site_position in enumerate(site_positions):
        print(f"\n尝试位点 {i+1}/{len(site_positions)}: {site_position}")
        
        # 定位分子
        positioned_mol = position_molecule_on_site(
            molecule, substrate, site_position, height, rotation, anchor_atom
        )
        
        # 检查并可能调整分子位置
        if check_position:
            need_adjust, positioned_mol = check_and_adjust_position(positioned_mol, substrate)
        
        # 合并分子和基底
        trial_system = substrate.copy()
        trial_system += positioned_mol
        
        # 检查并调整晶胞高度
        if min_vacuum_above > 0:
            trial_system = adjust_cell_height(trial_system, min_vacuum_above)
        
        # 如果需要验证
        if validate:
            current_valid, current_issues = validate_system(
                trial_system, 
                n_substrate_atoms=len(substrate), 
                scale_factor=scale_factor, 
                min_mol_surf_dist=min_mol_surf_dist,
                verbose=False
            )
            
            # 输出验证结果
            if current_valid:
                print(f"位点 {i+1} 验证通过: 吸附结构合理")
                is_valid = True
                best_system = trial_system
                best_issues = current_issues
                break  # 找到合适位点，停止尝试
            else:
                print(f"位点 {i+1} 验证失败: {', '.join(current_issues)}")
                
                # 如果是第一个位点或比之前的情况更好（问题更少），则更新最佳系统
                if best_system is None or len(current_issues) < len(best_issues):
                    best_system = trial_system
                    best_issues = current_issues
        else:
            # 不需要验证，直接使用当前系统
            best_system = trial_system
            is_valid = True
            break
    
    # 如果所有位点都尝试完毕
    if best_system is None:
        print("所有候选位点均不符合要求，无法创建合理的吸附系统")
        return None, False, ["无法创建合理的吸附系统"]
    
    # 如果需要验证但没有找到合适的位点
    if validate and not is_valid:
        print(f"警告: 未找到完全满足要求的吸附位点，使用最佳候选位点 (问题: {', '.join(best_issues)})")
    
    # 保存结构
    if output_file is not None and best_system is not None:
        write(output_file, best_system, format='vasp')
        print(f"吸附系统已保存至: {output_file}")
    
    return best_system, is_valid, best_issues

def visualize_validation(system, n_substrate_atoms=None, output_file=None):
    """
    可视化验证结果
    
    参数:
        system (ase.Atoms): 原子系统
        n_substrate_atoms (int): 基底原子的数量
        output_file (str): 输出文件路径
    """
    if n_substrate_atoms is None:
        # 尝试自动识别
        substrate_indices, molecule_indices = identify_molecule_substrate_atoms(system)
    else:
        substrate_indices = list(range(n_substrate_atoms))
        molecule_indices = list(range(n_substrate_atoms, len(system)))
    
    positions = system.get_positions()
    
    plt.figure(figsize=(10, 6))
    
    # 绘制所有原子的z坐标分布
    plt.subplot(121)
    substrate_positions = positions[substrate_indices]
    molecule_positions = positions[molecule_indices]
    
    plt.hist([substrate_positions[:, 2], molecule_positions[:, 2]], 
            bins=30, label=['基底', '分子'])
    
    max_z_substrate = np.max(substrate_positions[:, 2])
    min_z_molecule = np.min(molecule_positions[:, 2])
    
    plt.axvline(x=max_z_substrate, color='r', linestyle='--', label='基底表面')
    plt.axvline(x=min_z_molecule, color='g', linestyle='--', label='分子底部')
    plt.xlabel('Z坐标 (Å)')
    plt.ylabel('原子数量')
    plt.legend()
    plt.title('Z坐标分布')
    
    # 绘制原子间距离分布
    plt.subplot(122)
    # 计算分子原子到基底原子的所有距离
    all_distances = []
    for i in molecule_indices:
        for j in substrate_indices:
            dist = np.linalg.norm(positions[i] - positions[j])
            all_distances.append(dist)
    
    plt.hist(all_distances, bins=30)
    plt.xlabel('距离 (Å)')
    plt.ylabel('原子对数量')
    plt.title('分子-基底原子间距离分布')
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file)
        print(f"可视化结果已保存至: {output_file}")
    else:
        plt.show()

def main():
    """
    命令行入口点
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='将分子定位到基底的指定位置')
    parser.add_argument('--molecule', '-m', required=True, help='分子文件路径')
    parser.add_argument('--substrate', '-s', required=True, help='基底文件路径')
    parser.add_argument('--site', choices=['top', 'bridge', 'hollow'], help='吸附位点类型')
    parser.add_argument('--element', help='目标吸附位置的元素')
    parser.add_argument('--position', type=float, nargs=3, help='目标吸附位置的坐标')
    parser.add_argument('--distance', '--height', type=float, default=2.0, 
                        help='分子距离表面的高度(Å)')
    parser.add_argument('--rotation', type=float, nargs=3, default=[0, 0, 0],
                       help='分子的旋转角度(度)')
    parser.add_argument('--anchor', type=int, default=0, 
                       help='分子中作为锚点的原子索引')
    parser.add_argument('--min-vacuum-above', type=float, default=10.0,
                       help='分子顶部与晶胞顶部的最小距离(Å)，设为0则不调整晶胞高度，默认10.0')
    parser.add_argument('--no-check', action='store_true',
                       help='不检查分子位置是否在基底内')
    parser.add_argument('--validate', action='store_true', default=True,
                       help='验证吸附结构的合理性')
    parser.add_argument('--no-validate', action='store_false', dest='validate',
                      help='不验证吸附结构的合理性')
    parser.add_argument('--scale-factor', type=float, default=0.7,
                       help='原子重叠判断的缩放因子')
    parser.add_argument('--min-distance', type=float, default=1.5,
                       help='分子到表面的最小距离(Å)')
    parser.add_argument('--max-trial-sites', type=int, default=5,
                       help='尝试的最大位点数量')
    parser.add_argument('--visualize', '-v', action='store_true',
                       help='可视化验证结果')
    parser.add_argument('--output', '-o', default=None, help='输出文件路径')
    
    args = parser.parse_args()
    
    # 构建吸附位置规格
    site_spec = {}
    if args.element is not None:
        site_spec['element'] = args.element
    if args.position is not None:
        site_spec['position'] = args.position
    if args.site is not None:
        site_spec['site_type'] = args.site
    
    # 设置默认输出文件名
    if args.output is None:
        mol_base = os.path.splitext(os.path.basename(args.molecule))[0]
        sub_base = os.path.splitext(os.path.basename(args.substrate))[0]
        args.output = f"{mol_base}_on_{sub_base}.vasp"
    
    # 设置可视化输出文件名
    vis_file = None
    if args.visualize:
        output_base = os.path.splitext(args.output)[0]
        vis_file = f"{output_base}_validation.png"
    
    # 创建吸附系统
    system, is_valid, issues = create_adsorption_system(
        args.molecule,
        args.substrate,
        site_spec,
        args.distance,
        tuple(args.rotation),
        args.anchor,
        args.min_vacuum_above,
        args.output,
        not args.no_check,
        args.validate,
        args.scale_factor,
        args.min_distance,
        args.max_trial_sites
    )
    
    # 如果需要可视化并且系统创建成功
    if args.visualize and system is not None:
        visualize_validation(system, len(read(args.substrate)), vis_file)
    
    # 返回验证结果作为退出状态
    if not is_valid:
        print(f"警告: 最终结构存在问题: {', '.join(issues)}")
        return 1
    
    return 0

if __name__ == "__main__":
    main() 