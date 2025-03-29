"""Language configuration file for the Molecular Adsorption GUI"""

class LanguageManager:
    """全局语言管理器"""
    _instance = None
    _current_language = 'en'
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    @classmethod
    def get_current_language(cls):
        """获取当前语言"""
        return cls._current_language
    
    @classmethod
    def set_language(cls, language):
        """设置当前语言"""
        if language in ['en', 'zh']:
            cls._current_language = language
    
    @classmethod
    def get_text(cls, key):
        """获取指定key的翻译文本"""
        try:
            return TRANSLATIONS[cls._current_language][key]
        except KeyError:
            # 如果找不到翻译，返回英文版本，如果英文版本也没有，返回key
            return TRANSLATIONS['en'].get(key, key)

# 翻译字典
TRANSLATIONS = {
    'en': {
        # Window title
        'window_title': "Molecular Adsorption Simulation Tool",
        
        # Menu items
        'file': 'File',
        'view': 'View',
        'language': "Language",
        
        # File operations
        'file_operations': "File Operations",
        'molecule_file': "Molecule File:",
        'substrate_file': "Substrate File:",
        'browse': "Browse...",
        'import_molecule': "Import Molecule",
        'import_substrate': "Import Substrate",
        'import_state': "Import State",
        'export_state': "Export State",
        'export_structure': "Export Structure",
        'display_style': "Display Style:",
        'display_styles': ["Ball and Stick", "Space Filling", "Wireframe"],
        
        # Font related
        'chinese_font_setup_completed': "Chinese font setup completed",
        'error_setting_chinese_font': "Error setting Chinese font: {}",
        
        # Adsorption parameters
        'adsorption_parameters': "Adsorption Parameters",
        'adsorption_site': "Adsorption Site:",
        'site_types': ["Top", "Bridge", "Hollow"],
        'select_site': "Select Site",
        'target_element': "Target Element:",
        'distance_from_surface': "Distance from Surface (Å):",
        'minimum_vacuum': "Minimum Vacuum Height (Å):",
        
        # Rotation controls
        'rotation_controls': "Rotation Controls",
        'x_axis': "X-axis:",
        'y_axis': "Y-axis:",
        'z_axis': "Z-axis:",
        'reset_pose': "Reset Pose",
        
        # Translation adjustment
        'translation_adjustment': "Translation Adjustment",
        'check_xy_plane': "Check if molecule is within substrate XY plane",
        'auto_adjust': "Auto Adjust Molecule Position",
        'xy_offset': "XY Offset",
        'x_offset': "X Offset:",
        'y_offset': "Y Offset:",
        'apply_offset': "Apply Offset",
        
        # Action buttons
        'reset_system': "Reset Adsorption System",
        'validate_structure': "Validate Structure",
        
        # Status messages
        'status_ready': "Status: Ready",
        'status_loading': "Status: Loading...",
        'status_processing': "Status: Processing...",
        'status_completed': "Status: Completed",
        'status_error': "Status: Error",
        'status_selected_site': 'Selected {} site at {}',
        'error_occurred_when_creating_initial_molecule_structure': 'Error occurred when creating initial molecule structure',
        'status_applied_vertical_adjustment_msg': 'Applied vertical adjustment: {}',
        'vertical_adjustment_error_msg': 'Vertical adjustment error: {}',
        'status_vertical_adjustment_failed_msg': 'Vertical adjustment failed: {}',
        'status_no_vertical_adjustment_applied_msg': 'No vertical adjustment applied',
        'status_error_occurred_when_creating_initial_molecule_structure': 'Error occurred when creating initial molecule structure: {}',
        'status_created_adsorption_system_msg': 'Adsorption system created successfully',
        'error_occurred_when_creating_adsorption_system': 'Error occurred when creating adsorption system: {}',
        'error_occurred_when_importing_adsorption_state': 'Error occurred when importing adsorption state: {}',
        'no_valid_adsorption_coordinates_found': 'No valid adsorption coordinates found',
        'please_load_substrate_and_molecule_first': 'Please load substrate and molecule first',
        'please_select_adsorption_site_first': 'Please select adsorption site first',
        'no_structure_to_save': 'No structure to save',
        
        # View controls
        'view_controls': "View Controls",
        'zoom_in': "Zoom In",
        'zoom_out': "Zoom Out",
        'reset_view': "Reset View",
        'xy_projection': "XY Projection View",
        
        # Error messages
        'loading_error': "Loading Error",
        'error_loading_molecule': "Error loading molecule file: {}",
        'error_loading_substrate': "Error loading substrate file: {}",
        'structure_issues': "Structure Issues",
        'system_issues': "The adsorption system has the following issues:\n{}",
        'validation_error': "Validation Error",
        'error_validating_system': "Error validating adsorption system: {}",
        'validation_result': "Validation Result",
        'validation_passed': "Adsorption structure validation passed: Structure is reasonable",
        'validation_failed': "Adsorption structure validation failed:\n{}",
        'positioning_error': "Positioning Error",
        'error_positioning_molecule': "Error positioning molecule: {}",
        'plotting_structure_msg': "Plotting structure with {} atoms...",
        'first_few_positions_msg': "First few atom positions: {}",
        
        # Export related
        'warning': "Warning",
        'error': "Error",
        'create_system_first': "Please create adsorption system first",
        'select_adsorption_site_first': "Please select an adsorption site first",
        'export_success': "Export Success",
        'export_success_msg': "Structure has been successfully exported to directory:\n{}",
        'export_error': "Export Error",
        'export_error_msg': "Error occurred during export:\n{}",
        'saved_system_to': "Saved adsorption system to: {}",
        'saved_substrate_to': "Saved substrate structure to: {}",
        'saved_molecule_to': "Saved molecule structure to: {}",
        'saved_state_to': "Saved state to: {}",
        'images_saved_to_dir': "Images have been saved to {} directory",
        'completed': "Completed",
        'batch_site_processing_completed_msg': "Batch site processing completed, results saved in {} directory",
        'batch_anchor_processing_completed_msg': "Batch anchor processing completed, results saved in {} directory",
        'batch_path_processing_completed_msg': "Batch path processing completed, results saved in {} directory",
        
        # Additional UI elements
        'molecule_adjustment': "Molecule Adjustment",
        'anchor_atom': "Anchor Atom:",
        'select_by_structure': "Select by Structure",
        'keep_vertical': "Keep anchor point perpendicular to surface",
        'batch_processing': "Batch Processing",
        'site_batch': "Site Batch Processing",
        'anchor_batch': "Anchor Batch Processing",
        'path_batch': "Path Batch Processing",
        'save_current_view': "Save Current View",
        'view_direction': "View Direction:",
        'view_directions': ["Default", "A Axis", "B Axis", "C Axis", 
                          "A-B Plane", "B-C Plane", "A-C Plane"],
        
        # Projection view related
        'substrate': "Substrate",
        'molecule': "Molecule",
        'mirror': "Mirror",
        'anchor_point': "Anchor Point",
        'adsorption_site': "Adsorption Site",
        'lattice_height': "Lattice Height",
        'molecule_highest_point': "Molecule Highest Point",
        
        # Export options dialog
        'export_options': "Export Options",
        'base_name': "Base Name:",
        'export_content': "Export Content",
        'export_adsorption_system': "Export Adsorption System",
        'export_substrate_and_molecule_structures': "Export Substrate and Molecule Structures",
        'export_picture': "Export Picture",
        'export_current_state': "Export Current State",
        'file_format': "File Format",
        
        # Site selection dialog
        'select': "Select",
        'site_filter': "Site Filter",
        'enable_element_filter': "Enable Element Filter",
        'surface_layer': "Surface Layer:",
        'display_mode': "Display Mode:",
        'top_layer': "Top Layer Only",
        'multi_layer': "Multiple Layers",
        'available_elements': "Available Elements:",
        'add': "Add",
        'remove': "Remove",
        'clear': "Clear",
        'filter_elements': "Filter Elements:",
        'ready': "Ready",
        'cancel': "Cancel",
        'current_anchor_position': "Current Anchor Position",
        'site_type': "Site Type",
        'click_position': "Click Position",
        'available_sites_number': "Available Sites Number",
        'current_view': "Current View",
        'site': "Site",
        'position': "Position",
        'distance': "Distance",
        'selected_site_index': "Selected Site Index",
        'selected_site': "Selected Site",
        'element': "Element",
        
        # Atom selection dialog
        'select_anchor_atom': "Select Anchor Atom",
        'atom_index': "Atom Index:",
        'select_an_atom_to_view_information': "Select an atom to view information",
        'direction_adjustment_options': "Direction Adjustment Options",
        'anchor_perpendicular_to_molecule_center': "Anchor Perpendicular to Molecule Center",
        'anchor_perpendicular_to_C_atom': "Anchor Perpendicular to C Atom",
        'keep_current_molecule_state': "Keep Current Molecule State",
        'adjust_molecule_position': "Adjust Molecule Position",
        'ok': "OK",
        
        # System messages
        'cell_height_adjusted_msg': 'Cell height adjusted: +{:.2f}Å (New height: {:.2f}Å)',
        'error_occurred_when_adjusting_vacuum_height': 'Error occurred when adjusting vacuum height: {}',
        'position_error': 'Position Error',
        'adsorption_system_invalid_msg': 'The adsorption system is invalid for the following reasons:',
        
        # Reset related
        'reset_error': 'Reset Error',
        'reset_molecule_pose_error_msg': 'Error occurred when resetting molecule pose: {}',
        'status_reset_molecule_pose_msg': 'Reset molecule pose successfully',
        
        # Calculator related
        'calculator_title': "Simple Calculator",
        'calculator_description': "This script provides basic calculator functionality, supporting addition, subtraction, multiplication, division, power and modulo operations.",
        'calculator_operator_not_found': "No supported operator found (+, -, *, /, **, %)",
        'calculator_expression_error': "Expression format error: {}",
        'calculator_number_conversion_error': "Cannot convert expression parts to numbers: {}",
        'calculator_division_by_zero': "Divisor cannot be zero",
        'calculator_operation_names': {
            '+': "Addition",
            '-': "Subtraction",
            '*': "Multiplication",
            '/': "Division",
            '**': "Power",
            '%': "Modulo"
        },
        
        # Signal related
        'rotation_angle': "{}°",
        'anchor_element': "({})",
        'file_dialog_title': "Open {} File",
        'file_dialog_filter': "VASP Files (*.vasp);;XYZ Files (*.xyz);;CIF Files (*.cif);;All Files (*.*)",
        'status_loaded_molecule': "Status: Loaded molecule with {} atoms",
        'status_loaded_substrate': "Status: Loaded substrate with {} atoms",
        'warning_load_substrate_first': "Please load substrate structure first",
        
        # Import/Export state related
        'export_adsorption_state': "Export Adsorption State",
        'import_adsorption_state': "Import Adsorption State",
        'json_file_filter': "JSON Files (*.json);;All Files (*.*)",
        'state_imported_successfully': "State imported successfully",
        'state_exported_successfully': "State exported successfully to:\n{}",
        'status_state_exported': "Status: Adsorption state exported to {}",
        'status_state_imported': "Status: Adsorption state imported from {}",
        'error_occurred_when_exporting_adsorption_state': "Error occurred when exporting adsorption state: {}",
        
        # Rotation related
        'rotation_error': "Rotation Error",
        'error_applying_rotation': "Error applying rotation: {}",
        'status_rotation_applied': "Status: Rotation applied",
        'status_rotation_applied_with_vertical_adjustment': "Status: Rotation applied (keeping vertical adjustment)",
        
        # Auto adjust position related
        'surface_atoms_highest_z': "Surface atoms highest Z position",
        'molecule_lowest_z_atom_index': "Molecule lowest Z atom index",
        'current_molecule_distance_to_surface': "Current molecule distance to surface",
        'need_upward_shift': "Need upward shift",
        'auto_adjusted_molecule_height': "Auto adjusted molecule height, shift amount",
        'status_molecule_adjusted_to_safe_distance': "Status: Molecule lowest point adjusted to 1.5Å above surface",
        'molecule_surface_distance_reasonable': "Molecule distance to surface is reasonable, no adjustment needed",
        'error_auto_adjusting_position': "Error auto adjusting position",
        
        # Distance and vacuum related
        'status_distance_adjusted': "Status: Adjusted molecule distance from surface to {} Å",
        'error_updating_distance': "Error updating distance from surface",
        'status_vacuum_height_adjusted': "Status: Adjusted cell height to meet minimum vacuum height {} Å",
        'status_current_vacuum_height': "Status: Current vacuum height {:.2f} Å already meets minimum requirement {} Å",
        
        # Orientation related
        'orientation_error': "Orientation Error",
        'error_adjusting_molecule_orientation': "Error adjusting molecule orientation: {}",
        'error_updating_projection_views': "Error updating projection views",
        
        # Custom adsorption system creation related
        'positioning_info': "Positioning Information",
        'anchor_atom_index': "Anchor atom index",
        'selected_site_coordinates': "Selected site coordinates",
        'current_anchor_position': "Current anchor position",
        'height': "Height",
        'offset': "Offset",
        'positioning_error_msg': "Positioning error",
        'expected_position': "Expected position",
        'actual_position': "Actual position",
        'positioning_success_msg': "Positioning successful, anchor precisely placed at target position",
        'molecule_outside_substrate_xy_msg': "Molecule is outside substrate XY range, applying periodic adjustment",
        'anchor_deviation_after_adjustment_msg': "After adjustment, anchor has shifted by",
        'adjustment_completed_anchor_unchanged_msg': "Adjustment completed, anchor position unchanged",
        
        # Plotting related
        'no_atoms_to_plot': "No atoms provided for plot_structure",
        'saving_current_view_msg': "Saving current view: elevation={}, azimuth={}",
        
        # Anchor selection related
        'please_load_molecule_first': "Please load molecule structure first",
        'selection_error': "Selection Error",
        'error_selecting_anchor': "Error selecting anchor point: {}",
        
        # View direction
        'status_view_direction_changed': "Status: View direction changed to {}"
    },
    
    'zh': {
        # Window title
        'window_title': "分子吸附模拟工具",
        
        # Menu items
        'file': '文件',
        'view': '视图',
        'language': "语言",
        
        # File operations
        'file_operations': "文件操作",
        'molecule_file': "分子文件:",
        'substrate_file': "基底文件:",
        'browse': "浏览...",
        'import_molecule': "导入分子",
        'import_substrate': "导入基底",
        'import_state': "导入状态",
        'export_state': "导出状态",
        'export_structure': "导出结构",
        'display_style': "显示样式:",
        'display_styles': ["球棍模型", "空间填充", "线框模型"],
        
        # Font related
        'chinese_font_setup_completed': "设置中文字体完成",
        'error_setting_chinese_font': "设置中文字体时出错: {}",
        
        # Adsorption parameters
        'adsorption_parameters': "吸附参数",
        'adsorption_site': "吸附位点:",
        'site_types': ["顶位", "桥位", "空位"],
        'select_site': "选择位点",
        'target_element': "目标元素:",
        'distance_from_surface': "距表面距离 (Å):",
        'minimum_vacuum': "最小真空层高度 (Å):",
        
        # Rotation controls
        'rotation_controls': "旋转控制",
        'x_axis': "X轴:",
        'y_axis': "Y轴:",
        'z_axis': "Z轴:",
        'reset_pose': "重置姿态",
        
        # Translation adjustment
        'translation_adjustment': "平移调整",
        'check_xy_plane': "检查分子是否在基底XY平面内",
        'auto_adjust': "自动调整分子位置",
        'xy_offset': "XY偏移",
        'x_offset': "X偏移:",
        'y_offset': "Y偏移:",
        'apply_offset': "应用偏移",
        
        # Action buttons
        'reset_system': "重置吸附系统",
        'validate_structure': "验证结构",
        
        # Status messages
        'status_ready': "状态: 就绪",
        'status_loading': "状态：加载中...",
        'status_processing': "状态：处理中...",
        'status_completed': "状态：已完成",
        'status_error': "状态：错误",
        'status_selected_site': '已选择{}位点，位置：{}',
        'error_occurred_when_creating_initial_molecule_structure': '创建初始分子结构时发生错误',
        'status_applied_vertical_adjustment_msg': '已应用垂直调整：{}',
        'vertical_adjustment_error_msg': '垂直调整错误：{}',
        'status_vertical_adjustment_failed_msg': '垂直调整失败：{}',
        'status_no_vertical_adjustment_applied_msg': '未应用垂直调整',
        'status_error_occurred_when_creating_initial_molecule_structure': '创建初始分子结构时发生错误：{}',
        'status_created_adsorption_system_msg': '吸附系统创建成功',
        'error_occurred_when_creating_adsorption_system': '创建吸附系统时发生错误：{}',
        'error_occurred_when_importing_adsorption_state': '导入吸附状态时出错：{}',
        'no_valid_adsorption_coordinates_found': '未找到有效的吸附坐标',
        'please_load_substrate_and_molecule_first': '请先加载基底和分子',
        'please_select_adsorption_site_first': '请先选择吸附位点',
        'no_structure_to_save': '没有可保存的结构',
        
        # View controls
        'view_controls': "视图控制",
        'zoom_in': "放大",
        'zoom_out': "缩小",
        'reset_view': "重置视图",
        'xy_projection': "XY投影视图",
        
        # Error messages
        'loading_error': "加载错误",
        'error_loading_molecule': "加载分子文件时出错: {}",
        'error_loading_substrate': "加载基底文件时出错: {}",
        'structure_issues': "结构问题",
        'system_issues': "吸附系统存在以下问题:\n{}",
        'validation_error': "验证错误",
        'error_validating_system': "验证吸附系统时出错: {}",
        'validation_result': "验证结果",
        'validation_passed': "吸附结构验证通过: 结构合理",
        'validation_failed': "吸附结构验证失败:\n{}",
        'positioning_error': "定位错误",
        'error_positioning_molecule': "定位分子时出错: {}",
        'plotting_structure_msg': "正在绘制包含 {} 个原子的结构...",
        'first_few_positions_msg': "前几个原子位置: {}",
        
        # Export related
        'warning': "警告",
        'error': "错误",
        'create_system_first': "请先创建吸附系统",
        'select_adsorption_site_first': "请先选择吸附位点",
        'export_success': "导出成功",
        'export_success_msg': "结构已成功导出到目录：\n{}",
        'export_error': "导出错误",
        'export_error_msg': "导出过程中发生错误：\n{}",
        'saved_system_to': "已保存吸附系统结构到: {}",
        'saved_substrate_to': "已保存基底结构到: {}",
        'saved_molecule_to': "已保存分子结构到: {}",
        'saved_state_to': "已保存状态到: {}",
        'images_saved_to_dir': "图片已保存到 {} 目录",
        'completed': "完成",
        'batch_site_processing_completed_msg': "批量位点处理完成，结果保存在 {} 目录下",
        'batch_anchor_processing_completed_msg': "批量锚点处理完成，结果保存在 {} 目录下",
        'batch_path_processing_completed_msg': "批量路径处理完成，结果保存在 {} 目录下",
        
        # Additional UI elements
        'molecule_adjustment': "分子调整",
        'anchor_atom': "锚点原子:",
        'select_by_structure': "从结构中选择",
        'keep_vertical': "保持锚点垂直于表面",
        'batch_processing': "批处理",
        'site_batch': "位点批处理",
        'anchor_batch': "锚点批处理",
        'path_batch': "路径批处理",
        'save_current_view': "保存当前视图",
        'view_direction': "视角方向:",
        'view_directions': ["默认", "A轴", "B轴", "C轴", 
                          "A-B平面", "B-C平面", "A-C平面"],
        
        # Projection view related
        'substrate': "基底",
        'molecule': "分子",
        'mirror': "镜像",
        'anchor_point': "锚点",
        'adsorption_site': "吸附位点",
        'lattice_height': "晶格高度",
        'molecule_highest_point': "分子最高点",
        
        # Export options dialog
        'export_options': "导出选项",
        'base_name': "基础名称:",
        'export_content': "导出内容",
        'export_adsorption_system': "导出吸附系统",
        'export_substrate_and_molecule_structures': "导出基底和分子结构",
        'export_picture': "导出图片",
        'export_current_state': "导出当前状态",
        'file_format': "文件格式",
        
        # Site selection dialog
        'select': "选择",
        'site_filter': "位点过滤",
        'enable_element_filter': "启用元素过滤",
        'surface_layer': "表面层:",
        'display_mode': "显示模式:",
        'top_layer': "仅顶层",
        'multi_layer': "多层",
        'available_elements': "可用元素:",
        'add': "添加",
        'remove': "移除",
        'clear': "清除",
        'filter_elements': "过滤元素:",
        'ready': "就绪",
        'cancel': "取消",
        'current_anchor_position': "当前锚点位置",
        'site_type': "位点类型",
        'click_position': "点击位置",
        'available_sites_number': "可用位点数量",
        'current_view': "当前视图",
        'site': "位点",
        'position': "位置",
        'distance': "距离",
        'selected_site_index': "选中位点索引",
        'selected_site': "选中位点",
        'element': "元素",
        
        # Atom selection dialog
        'select_anchor_atom': "选择锚点原子",
        'atom_index': "原子索引:",
        'select_an_atom_to_view_information': "选择一个原子以查看信息",
        'direction_adjustment_options': "方向调整选项",
        'anchor_perpendicular_to_molecule_center': "锚点垂直于分子中心",
        'anchor_perpendicular_to_C_atom': "锚点垂直于C原子",
        'keep_current_molecule_state': "保持当前分子状态",
        'adjust_molecule_position': "调整分子位置",
        'ok': "确定",
        
        # System messages
        'cell_height_adjusted_msg': '晶胞高度已调整：+{:.2f}Å（新高度：{:.2f}Å）',
        'error_occurred_when_adjusting_vacuum_height': '调整真空层高度时发生错误：{}',
        'position_error': '位置错误',
        'adsorption_system_invalid_msg': '吸附系统无效，原因如下：',
        
        # Reset related
        'reset_error': '重置错误',
        'reset_molecule_pose_error_msg': '重置分子姿态时发生错误：{}',
        'status_reset_molecule_pose_msg': '已成功重置分子姿态',
        
        # Calculator related
        'calculator_title': "简单计算器",
        'calculator_description': "这个脚本提供基本的计算器功能，支持加减乘除、幂运算和模运算。",
        'calculator_operator_not_found': "未找到支持的操作符（+, -, *, /, **, %）",
        'calculator_expression_error': "表达式格式错误: {}",
        'calculator_number_conversion_error': "无法将表达式部分转换为数字: {}",
        'calculator_division_by_zero': "除数不能为零",
        'calculator_operation_names': {
            '+': "加法",
            '-': "减法",
            '*': "乘法",
            '/': "除法",
            '**': "幂运算",
            '%': "取模"
        },
        
        # Signal related
        'rotation_angle': "{}°",
        'anchor_element': "({})",
        'file_dialog_title': "打开{}文件",
        'file_dialog_filter': "VASP文件 (*.vasp);;XYZ文件 (*.xyz);;CIF文件 (*.cif);;所有文件 (*.*)",
        'status_loaded_molecule': "状态：已加载包含 {} 个原子的分子",
        'status_loaded_substrate': "状态：已加载包含 {} 个原子的基底",
        'warning_load_substrate_first': "请先加载基底结构",
        
        # Import/Export state related
        'export_adsorption_state': "导出吸附状态",
        'import_adsorption_state': "导入吸附状态",
        'json_file_filter': "JSON文件 (*.json);;所有文件 (*.*)",
        'state_imported_successfully': "状态导入成功",
        'state_exported_successfully': "吸附状态已成功导出到:\n{}",
        'status_state_exported': "状态：吸附状态已导出到 {}",
        'status_state_imported': "状态：吸附状态已从 {} 导入",
        'error_occurred_when_exporting_adsorption_state': "导出吸附状态时发生错误: {}",
        
        # Rotation related
        'rotation_error': "旋转错误",
        'error_applying_rotation': "应用旋转时出错: {}",
        'status_rotation_applied': "状态：已应用旋转",
        'status_rotation_applied_with_vertical_adjustment': "状态：已应用旋转（保持垂直调整）",
        
        # Auto adjust position related
        'surface_atoms_highest_z': "表面原子z轴最高位置",
        'molecule_lowest_z_atom_index': "分子z轴最低原子索引",
        'current_molecule_distance_to_surface': "当前分子最低点到表面的距离",
        'need_upward_shift': "需要向上移动",
        'auto_adjusted_molecule_height': "已自动调整分子高度，移动量",
        'status_molecule_adjusted_to_safe_distance': "状态：已将分子最低点调整到表面上方1.5Å处",
        'molecule_surface_distance_reasonable': "分子与表面距离合理，无需调整",
        'error_auto_adjusting_position': "自动调整位置时出错",
        
        # Distance and vacuum related
        'status_distance_adjusted': "状态：已调整分子到表面的距离为 {} Å",
        'error_updating_distance': "更新表面距离时出错",
        'status_vacuum_height_adjusted': "状态：已调整晶胞高度以满足最小真空高度 {} Å",
        'status_current_vacuum_height': "状态：当前真空高度 {:.2f} Å 已满足最小要求 {} Å",
        
        # Orientation related
        'orientation_error': "方向调整错误",
        'error_adjusting_molecule_orientation': "调整分子方向时出错: {}",
        'error_updating_projection_views': "更新投影视图时出错",
        
        # Custom adsorption system creation related
        'positioning_info': "定位信息",
        'anchor_atom_index': "锚点原子索引",
        'selected_site_coordinates': "选定位点坐标",
        'current_anchor_position': "当前锚点位置",
        'height': "高度",
        'offset': "偏移量",
        'positioning_error_msg': "定位误差",
        'expected_position': "预期位置",
        'actual_position': "实际位置",
        'positioning_success_msg': "定位成功，锚点精确放置在目标位置",
        'molecule_outside_substrate_xy_msg': "分子不在基底XY范围内，进行周期性调整",
        'anchor_deviation_after_adjustment_msg': "调整后锚点偏移了",
        'adjustment_completed_anchor_unchanged_msg': "调整完成，锚点位置保持不变",
        
        # Plotting related
        'no_atoms_to_plot': "没有提供原子进行plot_structure",
        'saving_current_view_msg': "正在保存当前视图：elevation={}, azimuth={}",
        
        # Anchor selection related
        'please_load_molecule_first': "请先加载分子结构",
        'selection_error': "选择错误",
        'error_selecting_anchor': "选择锚点时出错: {}",
        
        # View direction
        'status_view_direction_changed': "状态：视角方向已改变为 {}"
    }
} 