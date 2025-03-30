#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
分子吸附模拟GUI最简启动脚本 (英文版)
"""

import os
import sys
import traceback

# 获取当前脚本目录
current_dir = os.path.dirname(os.path.abspath(__file__))

# 确保可以找到模块
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

try:
    # 直接运行脚本文件
    script_path = os.path.join(current_dir, "molecular_adsorption", "adsorption_gui_eng.py")
    
    # 不修改工作目录，保持在项目根目录
    # original_dir = os.getcwd()
    # os.chdir(os.path.dirname(script_path))
    
    # 执行脚本
    with open(script_path, 'r', encoding='utf-8') as f:
        script_code = f.read()
    
    # 添加执行名字为主脚本
    namespace = {
        "__name__": "__main__",
        "__file__": script_path,
    }
    
    exec(script_code, namespace)
    
    # 不需要恢复工作目录
    # os.chdir(original_dir)
    
except ImportError as e:
    print(f"错误: 无法导入GUI模块: {e}")
    print("详细错误信息:")
    traceback.print_exc()
    sys.exit(1)
except Exception as e:
    print(f"启动GUI时发生错误: {e}")
    traceback.print_exc()
    sys.exit(1) 