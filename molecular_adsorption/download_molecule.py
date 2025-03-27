#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
从PubChem下载分子并构建包含分子的晶格结构
"""

import os
import numpy as np
import pubchempy as pcp
from ase import Atoms
from ase.io import write, read
from ase.visualize import view
from ase.data import covalent_radii
from ase.build import molecule

# 尝试导入RDKit，如果不可用则显示警告
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("警告: 未检测到RDKit模块。如果PubChem 3D结构不可用，将无法使用SMILES构建3D结构。")
    print("      建议安装RDKit: pip install rdkit")

def build_molecule_from_smiles(smiles):
    """
    使用RDKit从SMILES构建3D分子结构
    
    参数:
        smiles (str): 分子的SMILES表示
        
    返回:
        ase.Atoms: 分子结构，如果失败则返回None
    """
    if not RDKIT_AVAILABLE:
        print("无法从SMILES构建分子: RDKit未安装")
        return None
    
    try:
        # 创建RDKit分子对象
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"无法从SMILES创建分子: {smiles}")
            return None
        
        # 添加氢原子
        mol = Chem.AddHs(mol)
        
        # 生成3D构象
        success = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if success == -1:
            print(f"无法生成3D构象: {smiles}")
            return None
        
        # 优化几何结构
        AllChem.MMFFOptimizeMolecule(mol)
        
        # 转换为ASE Atoms对象
        atoms = []
        positions = []
        
        conf = mol.GetConformer()
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            idx = atom.GetIdx()
            pos = conf.GetAtomPosition(idx)
            atoms.append(symbol)
            positions.append([pos.x, pos.y, pos.z])
        
        # 创建ASE Atoms对象
        atoms_obj = Atoms(symbols=atoms, positions=positions)
        
        print(f"已从SMILES构建3D分子结构: {smiles}")
        return atoms_obj
        
    except Exception as e:
        print(f"从SMILES构建分子时出错: {str(e)}")
        return None

def try_get_3d_structure_from_pubchem(compound_id):
    """
    尝试从PubChem获取3D结构
    
    参数:
        compound_id (int): PubChem CID
        
    返回:
        ase.Atoms: 分子结构，如果失败则返回None
    """
    try:
        # 尝试获取3D结构
        properties = pcp.get_properties('IsomericSMILES,CanonicalSMILES,MolecularFormula,InChI', compound_id, 'cid')
        if not properties:
            print(f"无法获取分子属性: CID={compound_id}")
            return None
        
        # 尝试使用ASE内置的分子库
        formula = properties[0].get('MolecularFormula', '').lower()
        try:
            mol = molecule(formula)
            print(f"使用ASE内置分子库创建: {formula}")
            return mol
        except:
            print(f"ASE不支持该分子: {formula}")
        
        # 尝试直接从PubChem获取3D结构
        try:
            # 尝试获取3D结构
            mol3d = pcp.get_cids(compound_id, record_type='3d')
            if mol3d:
                sdf = pcp.download('SDF', mol3d[0], '3d')
                if sdf:
                    # 保存临时SDF文件并使用ASE读取
                    with open('temp_3d.sdf', 'w') as f:
                        f.write(sdf)
                    try:
                        mol = read('temp_3d.sdf')
                        os.remove('temp_3d.sdf')
                        print(f"从PubChem下载的3D结构已解析")
                        return mol
                    except Exception as e:
                        print(f"解析SDF文件时出错: {str(e)}")
                        if os.path.exists('temp_3d.sdf'):
                            os.remove('temp_3d.sdf')
        except Exception as e:
            print(f"获取3D结构时出错: {str(e)}")
        
        # 如果以上方法都失败，尝试使用SMILES构建
        smiles = properties[0].get('IsomericSMILES') or properties[0].get('CanonicalSMILES')
        if smiles:
            print(f"尝试从SMILES构建: {smiles}")
            return build_molecule_from_smiles(smiles)
            
        return None
        
    except Exception as e:
        print(f"从PubChem获取3D结构时出错: {str(e)}")
        return None

def download_molecule_by_name(name, output_dir='.'):
    """
    通过分子名称从PubChem下载分子结构
    
    参数:
        name (str): 分子名称或CID
        output_dir (str): 输出目录
        
    返回:
        ase.Atoms: 分子结构
    """
    try:
        print(f"正在查找分子: {name}")
        
        # 尝试将输入解析为CID
        try:
            cid = int(name)
            compounds = pcp.get_compounds(cid, 'cid')
        except ValueError:
            # 如果无法转换为整数，则按名称搜索
            compounds = pcp.get_compounds(name, 'name')
        
        if not compounds:
            print(f"未找到分子: {name}")
            return None
        
        # 获取第一个匹配的化合物
        compound = compounds[0]
        print(f"找到分子: {compound.iupac_name} (CID: {compound.cid})")
        
        # 尝试从内置数据库获取常见分子
        try:
            common_name = name.lower().strip()
            if common_name in ['h2o', 'water', '水']:
                mol = molecule('H2O')
                print(f"使用ASE内置分子库创建: 水分子(H2O)")
                return mol
            elif common_name in ['co', 'carbon monoxide', '一氧化碳']:
                mol = molecule('CO')
                print(f"使用ASE内置分子库创建: 一氧化碳(CO)")
                return mol
            elif common_name in ['co2', 'carbon dioxide', '二氧化碳']:
                mol = molecule('CO2')
                print(f"使用ASE内置分子库创建: 二氧化碳(CO2)")
                return mol
            elif common_name in ['ch4', 'methane', '甲烷']:
                mol = molecule('CH4')
                print(f"使用ASE内置分子库创建: 甲烷(CH4)")
                return mol
            elif common_name in ['nh3', 'ammonia', '氨']:
                mol = molecule('NH3')
                print(f"使用ASE内置分子库创建: 氨(NH3)")
                return mol
            # 可以添加更多常见分子...
        except Exception as e:
            print(f"尝试使用内置分子库时出错: {str(e)}")
        
        # 尝试获取3D结构
        mol = try_get_3d_structure_from_pubchem(compound.cid)
        if mol is not None:
            return mol
        
        # 如果以上方法都失败，则获取SMILES表示并使用RDKit构建
        smiles = compound.canonical_smiles
        print(f"SMILES: {smiles}")
        return build_molecule_from_smiles(smiles)
    
    except Exception as e:
        print(f"下载分子时出错: {str(e)}")
        return None

def create_molecular_crystal(molecule, vacuum=10.0, pbc=True):
    """
    将分子放入周期性晶格中
    
    参数:
        molecule (ase.Atoms): 分子结构
        vacuum (float): 真空层厚度(Å)
        pbc (bool/list): 周期性边界条件
        
    返回:
        ase.Atoms: 周期性晶格中的分子
    """
    if molecule is None:
        return None
        
    # 计算分子的包围盒
    positions = molecule.get_positions()
    min_pos = np.min(positions, axis=0)
    max_pos = np.max(positions, axis=0)
    
    # 计算分子尺寸
    mol_size = max_pos - min_pos
    
    # 创建足够大的晶胞，四周留有vacuum指定的空间
    cell = mol_size + 2 * vacuum
    
    # 将分子居中
    center = (min_pos + max_pos) / 2
    molecule.positions += (cell / 2) - center
    
    # 设置晶胞和周期性边界条件
    molecule.set_cell(cell)
    if isinstance(pbc, bool):
        molecule.set_pbc([pbc, pbc, pbc])
    else:
        molecule.set_pbc(pbc)
    
    return molecule

def save_molecule(molecule, filename, format='vasp'):
    """
    保存分子结构到文件
    
    参数:
        molecule (ase.Atoms): 分子结构
        filename (str): 输出文件名
        format (str): 输出格式
    """
    if molecule is None:
        print("没有分子结构可保存")
        return
        
    try:
        write(filename, molecule, format=format)
        print(f"分子结构已保存至: {filename}")
    except Exception as e:
        print(f"保存分子结构时出错: {str(e)}")

def main():
    """
    命令行入口点
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='从PubChem下载分子并创建晶格结构')
    parser.add_argument('molecule', help='分子名称或CID')
    parser.add_argument('--vacuum', type=float, default=10.0, help='真空层厚度(Å)')
    parser.add_argument('--pbc', type=str, default='TTT', help='周期性边界条件 (例如: TTT, FFF, TFT)')
    parser.add_argument('--output', '-o', default=None, help='输出文件名')
    parser.add_argument('--format', '-f', default='vasp', help='输出格式')
    
    args = parser.parse_args()
    
    # 解析周期性边界条件
    if args.pbc.lower() in ('t', 'true'):
        pbc = True
    elif args.pbc.lower() in ('f', 'false'):
        pbc = False
    else:
        pbc = [c.lower() == 't' for c in args.pbc[:3]]
    
    # 设置输出文件名
    if args.output is None:
        args.output = f"{args.molecule.replace(' ', '_')}.{args.format}"
    
    # 下载分子
    mol = download_molecule_by_name(args.molecule)
    
    if mol is not None:
        # 创建晶格结构
        crystal = create_molecular_crystal(mol, vacuum=args.vacuum, pbc=pbc)
        
        # 保存结构
        save_molecule(crystal, args.output, format=args.format)
    
if __name__ == "__main__":
    main()