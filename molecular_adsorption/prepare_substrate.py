 #!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
对基底(substrate)切面并添加真空层
"""

import os
import numpy as np
from ase import Atoms
from ase.io import write, read
from ase.build import surface

def prepare_substrate(substrate_file, miller_indices=(1, 0, 0), layers=4, vacuum=15.0, output_file=None):
    """
    从VASP文件读取基底结构，根据Miller指数切面，并添加真空层
    
    参数:
        substrate_file (str): 基底VASP文件路径
        miller_indices (tuple): Miller指数，例如(1,0,0)
        layers (int): 表面层数
        vacuum (float): 真空层厚度(Å)
        output_file (str): 输出文件名，默认为None时生成自动文件名
        
    返回:
        ase.Atoms: 处理后的基底表面结构
    """
    try:
        # 读取基底结构
        print(f"正在读取基底文件: {substrate_file}")
        bulk = read(substrate_file)
        
        # 创建表面切片
        slab = surface(bulk, miller_indices, layers, vacuum=vacuum)
        print(f"已创建表面切片: Miller指数={miller_indices}, 层数={layers}, 真空层={vacuum}Å")
        
        # 保存结构
        if output_file is None:
            # 自动生成输出文件名
            base_name = os.path.splitext(os.path.basename(substrate_file))[0]
            indices_str = ''.join(map(str, miller_indices))
            output_file = f"{base_name}_{indices_str}_{layers}layers_{vacuum}A.vasp"
        
        write(output_file, slab, format='vasp')
        print(f"表面结构已保存至: {output_file}")
        
        return slab
        
    except Exception as e:
        print(f"处理基底表面时出错: {str(e)}")
        return None

def add_vacuum_only(substrate_file, vacuum=15.0, direction=2, output_file=None):
    """
    仅在特定方向添加真空层，不切表面
    
    参数:
        substrate_file (str): 基底VASP文件路径
        vacuum (float): 真空层厚度(Å)
        direction (int): 添加真空层的方向，0=x, 1=y, 2=z
        output_file (str): 输出文件名
        
    返回:
        ase.Atoms: 添加真空层后的结构
    """
    try:
        # 读取结构
        print(f"正在读取文件: {substrate_file}")
        atoms = read(substrate_file)
        
        # 获取原始晶胞尺寸
        cell = atoms.get_cell()
        
        # 获取原子坐标
        positions = atoms.get_positions()
        
        # 计算在指定方向上的最小和最大原子坐标
        min_pos = np.min(positions[:, direction])
        max_pos = np.max(positions[:, direction])
        
        # 原子在该方向上的跨度
        span = max_pos - min_pos
        
        # 计算总晶胞长度 = 原子跨度 + 真空层厚度
        cell_length = span + vacuum
        
        # 更新晶胞
        new_cell = cell.copy()
        new_cell[direction, direction] = cell_length
        atoms.set_cell(new_cell)
        
        # 将原子居中放置
        shift = np.zeros(3)
        shift[direction] = (cell_length - span) / 2 - min_pos
        atoms.positions += shift
        
        # 保存结构
        if output_file is None:
            # 自动生成输出文件名
            base_name = os.path.splitext(os.path.basename(substrate_file))[0]
            directions = ['x', 'y', 'z']
            output_file = f"{base_name}_{vacuum}A_vacuum_{directions[direction]}.vasp"
        
        write(output_file, atoms, format='vasp')
        print(f"添加真空层后的结构已保存至: {output_file}")
        
        return atoms
        
    except Exception as e:
        print(f"添加真空层时出错: {str(e)}")
        return None

def main():
    """
    命令行入口点
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='处理基底表面并添加真空层')
    parser.add_argument('substrate_file', help='基底VASP文件路径')
    parser.add_argument('--miller', type=int, nargs=3, default=[1, 0, 0], 
                        help='Miller指数，例如: 1 0 0')
    parser.add_argument('--layers', type=int, default=4, help='表面层数')
    parser.add_argument('--vacuum', type=float, default=15.0, help='真空层厚度(Å)')
    parser.add_argument('--no-surface', action='store_true', 
                        help='不切表面，只添加真空层')
    parser.add_argument('--direction', type=int, choices=[0, 1, 2], default=2,
                        help='添加真空层的方向，0=x, 1=y, 2=z')
    parser.add_argument('--output', '-o', default=None, help='输出文件名')
    
    args = parser.parse_args()
    
    if args.no_surface:
        # 仅添加真空层
        add_vacuum_only(args.substrate_file, args.vacuum, args.direction, args.output)
    else:
        # 切表面并添加真空层
        prepare_substrate(args.substrate_file, tuple(args.miller), args.layers, 
                         args.vacuum, args.output)

if __name__ == "__main__":
    main()