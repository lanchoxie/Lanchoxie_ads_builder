#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
分子吸附模拟GUI启动器 - 语言选择
Molecular Adsorption GUI Launcher - Language Selection
"""

import os
import sys
import traceback
import tkinter as tk
from tkinter import messagebox
import subprocess

def start_chinese_version():
    """启动中文版GUI"""
    try:
        script_path = os.path.join(current_dir, "simple_start_chn.py")
        root.destroy()  # 关闭选择窗口
        subprocess.call([sys.executable, script_path])
    except Exception as e:
        messagebox.showerror("错误 / Error", f"启动失败: {str(e)}\nLaunch failed: {str(e)}")
        traceback.print_exc()

def start_english_version():
    """启动英文版GUI"""
    try:
        script_path = os.path.join(current_dir, "simple_start_eng.py")
        root.destroy()  # 关闭选择窗口
        subprocess.call([sys.executable, script_path])
    except Exception as e:
        messagebox.showerror("错误 / Error", f"启动失败: {str(e)}\nLaunch failed: {str(e)}")
        traceback.print_exc()

# 获取当前脚本目录
current_dir = os.path.dirname(os.path.abspath(__file__))

# 确保可以找到模块
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

# 创建选择窗口
root = tk.Tk()
root.title("语言选择 / Language Selection")
root.geometry("400x200")
root.resizable(False, False)

# 设置窗口图标
try:
    icon_path = os.path.join(current_dir, "molecular_adsorption", "icon.ico")
    if os.path.exists(icon_path):
        root.iconbitmap(icon_path)
except Exception:
    pass  # 忽略图标设置错误

# 顶部文本
header_label = tk.Label(
    root, 
    text="分子吸附模拟工具\nMolecular Adsorption Simulation Tool", 
    font=("Arial", 14, "bold"),
    pady=10
)
header_label.pack()

# 提示文本
prompt_label = tk.Label(
    root, 
    text="请选择界面语言 / Please select interface language:", 
    font=("Arial", 10),
    pady=5
)
prompt_label.pack()

# 按钮框架
button_frame = tk.Frame(root)
button_frame.pack(pady=10)

# 中文按钮
chinese_button = tk.Button(
    button_frame, 
    text="中文", 
    font=("Arial", 12, "bold"),
    width=10, 
    height=2,
    command=start_chinese_version
)
chinese_button.grid(row=0, column=0, padx=20)

# 英文按钮
english_button = tk.Button(
    button_frame, 
    text="English", 
    font=("Arial", 12, "bold"),
    width=10, 
    height=2,
    command=start_english_version
)
english_button.grid(row=0, column=1, padx=20)

# 启动主循环
if __name__ == "__main__":
    try:
        root.mainloop()
    except Exception as e:
        print(f"启动选择器时发生错误: {e}")
        print("Launcher error occurred: {e}")
        traceback.print_exc()
        sys.exit(1) 