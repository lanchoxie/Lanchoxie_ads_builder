#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
分子吸附模拟主脚本
整合所有功能，提供完整的工作流程
"""

import os
import argparse
import sys
import textwrap
import numpy as np
from ase.io import read, write

# 导入子模块
try:
    from molecular_adsorption.download_molecule import download_molecule_by_name, create_molecular_crystal, save_molecule
    from molecular_adsorption.prepare_substrate import prepare_substrate, add_vacuum_only
    # 不再需要单独导入check_position模块
    from molecular_adsorption.position_molecule import create_adsorption_system, check_and_adjust_position
    from molecular_adsorption.validate_adsorption import validate_adsorption_structure
except ImportError:
    # 如果作为独立脚本运行，则尝试直接导入
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    try:
        # 先尝试使用相对导入
        from .download_molecule import download_molecule_by_name, create_molecular_crystal, save_molecule
        from .prepare_substrate import prepare_substrate, add_vacuum_only
        from .position_molecule import create_adsorption_system, check_and_adjust_position
        from .validate_adsorption import validate_adsorption_structure
    except (ImportError, ValueError):
        # 如果相对导入失败，回退到绝对导入
        from download_molecule import download_molecule_by_name, create_molecular_crystal, save_molecule
        from prepare_substrate import prepare_substrate, add_vacuum_only
        from position_molecule import create_adsorption_system, check_and_adjust_position
        from validate_adsorption import validate_adsorption_structure

def create_working_dirs(base_dir='adsorption_work'):
    """创建工作目录结构"""
    dirs = {
        'molecules': os.path.join(base_dir, 'molecules'),
        'substrates': os.path.join(base_dir, 'substrates'),
        'adsorption': os.path.join(base_dir, 'adsorption'),
        'analysis': os.path.join(base_dir, 'analysis')
    }
    
    for name, dir_path in dirs.items():
        os.makedirs(dir_path, exist_ok=True)
        print(f"创建目录: {dir_path}")
    
    return dirs

def download_molecule_workflow(args, dirs):
    """分子下载工作流"""
    print(f"\n=== 从PubChem下载分子: {args.molecule} ===")
    
    # 下载分子
    mol = download_molecule_by_name(args.molecule)
    
    if mol is None:
        print(f"下载分子失败: {args.molecule}")
        return None
    
    # 创建晶格结构
    crystal = create_molecular_crystal(mol, vacuum=args.vacuum, pbc=True)
    
    # 生成输出文件名
    output_file = os.path.join(dirs['molecules'], 
                              f"{args.molecule.replace(' ', '_')}.vasp")
    
    # 保存结构
    save_molecule(crystal, output_file, format='vasp')
    
    return output_file

def prepare_substrate_workflow(args, dirs):
    """基底处理工作流"""
    print(f"\n=== 处理基底: {args.substrate} ===")
    
    input_file = args.substrate
    
    # 生成输出文件名
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    
    if args.no_surface:
        # 仅添加真空层
        output_file = os.path.join(dirs['substrates'], 
                                  f"{base_name}_vacuum{args.vacuum}A.vasp")
        
        substrate = add_vacuum_only(input_file, vacuum=args.vacuum, 
                                   direction=args.direction, 
                                   output_file=output_file)
    else:
        # 切表面并添加真空层
        if args.surface:
            # 如果提供的是表面编号，如111、100等
            surface = args.surface
        else:
            # 否则使用Miller指数
            indices_str = ''.join(map(str, args.miller))
            surface = indices_str
            
        output_file = os.path.join(dirs['substrates'], 
                                  f"{base_name}_{surface}_{args.layers}l_vac{args.vacuum}A.vasp")
        
        substrate = prepare_substrate(input_file, 
                                     surface=surface, 
                                     layers=args.layers, 
                                     vacuum=args.vacuum, 
                                     output_file=output_file)
    
    return output_file

def position_molecule_workflow(args, dirs):
    """分子定位工作流"""
    print(f"\n=== 将分子定位到基底上 ===")
    
    # 构建吸附位置规格
    site_spec = {}
    if args.element is not None:
        site_spec['element'] = args.element
    if args.position is not None:
        site_spec['position'] = args.position
    
    # 设置输出文件名
    mol_base = os.path.splitext(os.path.basename(args.molecule_file))[0]
    sub_base = os.path.splitext(os.path.basename(args.substrate_file))[0]
    output_file = os.path.join(dirs['adsorption'], 
                              f"{mol_base}_on_{sub_base}.vasp")
    
    # 创建吸附系统
    system = create_adsorption_system(
        args.molecule_file,
        args.substrate_file,
        site_spec,
        args.height,
        tuple(args.rotation),
        args.anchor,
        args.min_vacuum_above,
        output_file,
        check_position=not args.no_check
    )
    
    return output_file

def validate_workflow(args, dirs):
    """吸附结构验证工作流"""
    print(f"\n=== 验证吸附结构 ===")
    
    # 设置可视化输出文件名
    if args.visualize:
        base_name = os.path.splitext(os.path.basename(args.system_file))[0]
        vis_file = os.path.join(dirs['analysis'], f"{base_name}_validation.png")
    else:
        vis_file = None
    
    # 验证吸附结构
    is_valid, issues = validate_adsorption_structure(
        args.system_file,
        args.n_substrate,
        args.tolerance,
        vis_file
    )
    
    if is_valid:
        print("吸附结构验证通过")
    else:
        print("吸附结构存在以下问题:")
        for issue in issues:
            print(f"  - {issue}")
    
    return is_valid

def full_workflow(args):
    """完整工作流程"""
    # 创建工作目录
    dirs = create_working_dirs(args.work_dir)
    
    # 1. 下载分子
    if args.molecule and not args.molecule_file:
        molecule_file = download_molecule_workflow(args, dirs)
        if molecule_file is None:
            return
        args.molecule_file = molecule_file
    
    # 2. 处理基底
    if args.substrate and not args.substrate_file:
        substrate_file = prepare_substrate_workflow(args, dirs)
        if substrate_file is None:
            return
        args.substrate_file = substrate_file
    
    # 3. 将分子定位到基底上
    if args.molecule_file and args.substrate_file and not args.system_file:
        system_file = position_molecule_workflow(args, dirs)
        if system_file is None:
            return
        args.system_file = system_file
    
    # 4. 验证吸附结构
    if args.validate and args.system_file:
        is_valid = validate_workflow(args, dirs)
        
        if not is_valid and args.auto_adjust:
            print("\n=== 自动调整吸附结构 ===")
            # 这里可以添加自动调整逻辑
            # ...
            print("自动调整功能尚未实现")
    
    print("\n=== 工作流程完成 ===")

def main():
    parser = argparse.ArgumentParser(
        description='分子吸附模拟工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''
        示例用法:
        -------------------------------------------
        # 完整工作流程
        python main_adsorption.py --molecule CO --substrate Pt.cif --surface 111 --site top --distance 2.0 --min-vacuum-above 7.0
        
        # 仅下载分子
        python main_adsorption.py --molecule CO --vacuum 10 --only-download
        
        # 仅处理基底
        python main_adsorption.py --substrate Pt.cif --surface 111 --layers 4 --vacuum 15 --only-substrate
        
        # 将已有分子定位到基底上
        python main_adsorption.py --molecule-file CO.vasp --substrate-file Pt111.vasp --site top --distance 2.0 --min-vacuum-above 10.0
        
        # 验证现有的吸附结构
        python main_adsorption.py --system-file CO_on_Pt111.vasp --validate --visualize
        '''))
    
    # 工作目录
    parser.add_argument('--work-dir', default='adsorption_work', 
                       help='工作目录')
    
    # 分子相关参数
    group_mol = parser.add_argument_group('分子参数')
    group_mol.add_argument('--molecule', help='分子名称或CID')
    group_mol.add_argument('--molecule-file', help='分子文件路径')
    group_mol.add_argument('--vacuum', type=float, default=10.0, 
                          help='晶胞中的真空层厚度(Å)')
    group_mol.add_argument('--only-download', action='store_true',
                          help='仅下载分子，不执行后续步骤')
    
    # 基底相关参数
    group_sub = parser.add_argument_group('基底参数')
    group_sub.add_argument('--substrate', help='基底文件路径')
    group_sub.add_argument('--substrate-file', help='处理后的基底文件路径')
    group_sub.add_argument('--miller', type=int, nargs=3, default=[1, 0, 0], 
                          help='Miller指数，例如: 1 0 0')
    group_sub.add_argument('--surface', 
                          help='表面编号，例如: 111, 100等')
    group_sub.add_argument('--layers', type=int, default=4, 
                          help='表面层数')
    group_sub.add_argument('--no-surface', action='store_true',
                          help='不切割表面，仅添加真空层')
    group_sub.add_argument('--direction', type=int, default=2, 
                          help='添加真空层的方向，0=x，1=y，2=z')
    group_sub.add_argument('--only-substrate', action='store_true',
                          help='仅处理基底，不执行后续步骤')
    
    # 分子定位参数
    group_pos = parser.add_argument_group('分子定位参数')
    group_pos.add_argument('--site', choices=['top', 'bridge', 'hollow'], 
                          help='吸附位点类型')
    group_pos.add_argument('--element', help='目标吸附位置的元素')
    group_pos.add_argument('--position', type=float, nargs=3, 
                          help='目标吸附位置的坐标')
    group_pos.add_argument('--distance', '--height', dest='height', type=float, default=2.0, 
                          help='分子距离表面的高度(Å)')
    group_pos.add_argument('--rotation', type=float, nargs=3, default=[0, 0, 0], 
                          help='分子的旋转角度(度)')
    group_pos.add_argument('--anchor', type=int, default=0, 
                          help='分子中作为锚点的原子索引')
    group_pos.add_argument('--min-vacuum-above', type=float, default=10.0, 
                          help='分子顶部与晶胞顶部的最小距离(Å)，设为0则不调整晶胞高度')
    group_pos.add_argument('--no-check', action='store_true',
                          help='不检查分子位置是否在基底内')
    
    # 验证参数
    group_val = parser.add_argument_group('验证参数')
    group_val.add_argument('--system-file', help='吸附系统文件路径')
    group_val.add_argument('--validate', action='store_true', 
                          help='验证吸附结构')
    group_val.add_argument('--n-substrate', type=int, 
                          help='基底原子数量，用于区分基底和分子')
    group_val.add_argument('--tolerance', type=float, default=0.7, 
                          help='原子间距容许误差(倍数)')
    group_val.add_argument('--visualize', action='store_true', 
                          help='可视化验证结果')
    group_val.add_argument('--auto-adjust', action='store_true', 
                          help='自动调整问题结构')
    
    args = parser.parse_args()
    
    # 根据命令判断执行哪个工作流
    if args.only_download and args.molecule:
        # 仅下载分子
        dirs = create_working_dirs(args.work_dir)
        download_molecule_workflow(args, dirs)
    elif args.only_substrate and args.substrate:
        # 仅处理基底
        dirs = create_working_dirs(args.work_dir)
        prepare_substrate_workflow(args, dirs)
    elif args.system_file and args.validate:
        # 仅验证结构
        dirs = create_working_dirs(args.work_dir)
        validate_workflow(args, dirs)
    else:
        # 执行完整工作流
        full_workflow(args)

if __name__ == "__main__":
    main() 