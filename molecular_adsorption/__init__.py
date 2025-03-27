"""
分子吸附模拟工具包
用于从PubChem下载分子，构建吸附结构，以及进行吸附结构的验证
"""

__version__ = "0.1.0"

# 确保相对导入可以正常工作
try:
    # 导出主要模块
    from . import download_molecule
    from . import prepare_substrate
    from . import position_molecule
    from . import validate_adsorption
    from . import main_adsorption
    from . import adsorption_gui

    # 导出主要函数
    from .download_molecule import download_molecule_by_name, create_molecular_crystal, save_molecule
    from .prepare_substrate import prepare_substrate, add_vacuum_only
    from .position_molecule import (
        create_adsorption_system, 
        find_surface_atoms,
        rotate_molecule,
        check_and_adjust_position,
        validate_system
    )
    from .validate_adsorption import validate_adsorption_structure
    from .adsorption_gui import main as run_gui, AdsorptionGUI, MplCanvas, ELEMENT_COLORS, DEFAULT_COLOR
    from .batch_dialogs import BatchSiteSelectionDialog, BatchAnchorSelectionDialog, BatchPathDialog
    
except (ImportError, ValueError):
    # 当不能使用相对导入时不做任何事
    pass

__all__ = ['AdsorptionGUI', 'BatchSiteSelectionDialog', 'BatchAnchorSelectionDialog', 'BatchPathDialog']