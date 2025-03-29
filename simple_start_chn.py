#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
分子吸附模拟GUI最简启动脚本
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
    # 直接导入GUI模块并执行main函数
    from molecular_adsorption.adsorption_gui_chn import main
    main()
except ImportError as e:
    print(f"错误: 无法导入GUI模块: {e}")
    print("详细错误信息:")
    traceback.print_exc()
    sys.exit(1)
except Exception as e:
    print(f"启动GUI时发生错误: {e}")
    traceback.print_exc()
    sys.exit(1) 