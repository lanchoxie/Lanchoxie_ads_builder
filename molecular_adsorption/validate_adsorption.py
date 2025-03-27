#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
判断分子吸附结构是否合理
"""

import os
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.data import covalent_radii, atomic_numbers
import matplotlib.pyplot as plt
from scipy.spatial import KDTree

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
        list: 重叠原子对的列表 [(i, j, distance), ...]
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

def identify_molecule_substrate_atoms(system, n_substrate_atoms):
    """
    将系统中的原子分为基底原子和分子原子
    
    参数:
        system (ase.Atoms): 原子系统
        n_substrate_atoms (int): 基底原子的数量
        
    返回:
        tuple: (基底原子索引列表, 分子原子索引列表)
    """
    substrate_indices = list(range(n_substrate_atoms))
    molecule_indices = list(range(n_substrate_atoms, len(system)))
    
    return substrate_indices, molecule_indices

def validate_adsorption_structure(system_file, n_substrate_atoms=None, scale_factor=0.7,
                       min_mol_surf_dist=1.5, visualization_file=None):
    """
    验证吸附结构的合理性
    
    参数:
        system_file (str): 吸附系统文件路径
        n_substrate_atoms (int): 基底原子的数量，如不指定则尝试自动识别
        scale_factor (float): 重叠判断的缩放因子
        min_mol_surf_dist (float): 分子到表面的最小距离
        visualization_file (str): 可视化输出文件路径
        
    返回:
        tuple: (是否合理(bool), 问题列表(list))
    """
    # 读取系统
    system = read(system_file)
    
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
    
    # 识别基底和分子原子
    substrate_indices, molecule_indices = identify_molecule_substrate_atoms(system, n_substrate_atoms)
    
    if not molecule_indices:
        print("错误: 未能识别到分子原子")
        return False, ["未能识别到分子原子"]
    
    # 获取原子位置和符号
    positions = system.get_positions()
    symbols = system.get_chemical_symbols()
    
    # 检查1: 原子重叠
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
    
    # 收集问题
    issues = []
    
    # 打印重叠信息
    if mol_mol_overlaps:
        msg = f"发现 {len(mol_mol_overlaps)} 个分子内部原子重叠"
        print(f"警告: {msg}")
        issues.append(msg)
        for i, j, dist, min_dist in mol_mol_overlaps:
            print(f"  原子 {i}({symbols[i]}) 和 {j}({symbols[j]}) 距离 {dist:.3f}Å < {min_dist:.3f}Å")
    
    if sub_sub_overlaps:
        msg = f"发现 {len(sub_sub_overlaps)} 个基底内部原子重叠"
        print(f"警告: {msg}")
        # 基底内部重叠通常不是主要关注点，不添加到issues中
    
    if mol_sub_overlaps:
        msg = f"发现 {len(mol_sub_overlaps)} 个分子-基底原子重叠"
        print(f"警告: {msg}")
        issues.append(msg)
        for i, j, dist, min_dist in mol_sub_overlaps:
            print(f"  原子 {i}({symbols[i]}) 和 {j}({symbols[j]}) 距离 {dist:.3f}Å < {min_dist:.3f}Å")
    
    # 检查2: 分子到表面的距离
    print("检查分子到表面的距离...")
    
    # 获取基底表面原子
    substrate_positions = positions[substrate_indices]
    max_z_substrate = np.max(substrate_positions[:, 2])
    
    # 获取分子底部原子
    molecule_positions = positions[molecule_indices]
    min_z_molecule = np.min(molecule_positions[:, 2])
    
    # 计算分子到表面的距离
    mol_surf_distance = min_z_molecule - max_z_substrate
    
    print(f"分子到表面的距离: {mol_surf_distance:.3f}Å")
    
    if mol_surf_distance < min_mol_surf_dist:
        msg = f"分子到表面的距离 ({mol_surf_distance:.3f}Å) 小于最小要求 ({min_mol_surf_dist:.3f}Å)"
        print(f"警告: {msg}")
        issues.append(msg)
    
    # 可视化
    if visualization_file:
        plt.figure(figsize=(10, 6))
        
        # 绘制所有原子的z坐标分布
        plt.subplot(121)
        plt.hist([substrate_positions[:, 2], molecule_positions[:, 2]], 
                bins=30, label=['基底', '分子'])
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
        plt.savefig(visualization_file)
        print(f"可视化结果已保存至: {visualization_file}")
    
    # 结构合理性判断
    is_valid = len(issues) == 0
    
    if is_valid:
        print("结构验证通过: 吸附结构合理")
    else:
        print("结构验证失败: 吸附结构存在问题")
    
    return is_valid, issues

def main():
    """
    命令行入口点
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='验证分子吸附结构的合理性')
    parser.add_argument('system_file', help='吸附系统文件路径')
    parser.add_argument('--n-substrate', type=int, default=None, 
                       help='基底原子的数量')
    parser.add_argument('--scale-factor', type=float, default=0.7, 
                       help='重叠判断的缩放因子')
    parser.add_argument('--min-distance', type=float, default=1.5, 
                       help='分子到表面的最小距离(Å)')
    parser.add_argument('--visualize', '-v', action='store_true', 
                       help='是否生成可视化图表')
    parser.add_argument('--output', '-o', default=None, 
                       help='可视化输出文件路径')
    
    args = parser.parse_args()
    
    # 设置默认可视化输出文件名
    if args.visualize and args.output is None:
        base_name = os.path.splitext(os.path.basename(args.system_file))[0]
        args.output = f"{base_name}_validation.png"
    elif not args.visualize:
        args.output = None
    
    # 验证吸附结构
    validate_adsorption_structure(
        args.system_file,
        args.n_substrate,
        args.scale_factor,
        args.min_distance,
        args.output
    )

if __name__ == "__main__":
    main() 