"""
AutoDock Vina受体准备模块
生成PDBQT格式受体文件和对接盒子配置
"""

import os
import subprocess
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass


@dataclass
class VinaBox:
    """Vina对接盒子"""
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    
    def to_config_string(self) -> str:
        """生成Vina配置字符串"""
        return f"""center_x = {self.center_x:.3f}
center_y = {self.center_y:.3f}
center_z = {self.center_z:.3f}
size_x = {self.size_x:.1f}
size_y = {self.size_y:.1f}
size_z = {self.size_z:.1f}
"""


class VinaPreparator:
    """Vina受体准备器"""
    
    # 环肽尺寸估算 (用于确定盒子大小)
    # 21个氨基酸的环肽，直径约25-35Å
    PEPTIDE_DIAMETER = 30.0  # Å
    
    def __init__(self, config: Dict[str, Any], logger=None):
        self.config = config
        self.logger = logger
        
        # Vina可执行文件路径
        self.vina_executable = config.get('vina', {}).get('executable', 'vina')
    
    def prepare_receptor(self, clean_pdb: str, output_pdbqt: str) -> bool:
        """
        准备Vina受体文件 (PDBQT格式)
        
        步骤:
        1. 合并非极性氢
        2. 分配Gasteiger电荷
        3. 转换原子类型为AutoDock类型
        
        Args:
            clean_pdb: 清洁后的PDB文件
            output_pdbqt: 输出PDBQT文件路径
        
        Returns:
            是否成功
        """
        if not os.path.exists(clean_pdb):
            if self.logger:
                self.logger.error(f"输入文件不存在: {clean_pdb}")
            return False
        
        if self.logger:
            self.logger.info(f"准备Vina受体: {clean_pdb}")
        
        os.makedirs(os.path.dirname(output_pdbqt), exist_ok=True)
        
        # 方法1: 使用OpenBabel (推荐)
        if self._prepare_with_openbabel(clean_pdb, output_pdbqt):
            return True
        
        # 方法2: 使用AutoDockTools的prepare_receptor4.py
        if self._prepare_with_adt(clean_pdb, output_pdbqt):
            return True
        
        # 方法3: 简化转换 (可能不准确)
        if self.logger:
            self.logger.warning("【VINA_PREP_FALLBACK】")
            self.logger.warning("占位符含义: 使用简化PDBQT转换")
        
        return self._prepare_fallback(clean_pdb, output_pdbqt)
    
    def _prepare_with_openbabel(self, input_pdb: str, output_pdbqt: str) -> bool:
        """使用OpenBabel准备受体"""
        try:
            # OpenBabel命令: 合并非极性氢，分配部分电荷
            cmd = [
                'obabel', 
                '-i', 'pdb', input_pdb,
                '-o', 'pdbqt',
                '-O', output_pdbqt,
                '-xr'  # 删除残基类型记录 (受体模式)
            ]
            
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                timeout=60
            )
            
            if result.returncode == 0 and os.path.exists(output_pdbqt):
                if self.logger:
                    self.logger.info(f"OpenBabel转换完成: {output_pdbqt}")
                return True
            else:
                if self.logger:
                    self.logger.debug(f"OpenBabel失败: {result.stderr}")
                return False
                
        except FileNotFoundError:
            if self.logger:
                self.logger.debug("OpenBabel未安装")
            return False
        except Exception as e:
            if self.logger:
                self.logger.debug(f"OpenBabel异常: {e}")
            return False
    
    def _prepare_with_adt(self, input_pdb: str, output_pdbqt: str) -> bool:
        """使用AutoDockTools准备受体"""
        try:
            # 查找prepare_receptor4.py
            adt_paths = [
                '/usr/bin/prepare_receptor4.py',
                '/usr/local/bin/prepare_receptor4.py',
                '/opt/mgltools/bin/prepare_receptor4.py',
            ]
            
            prepare_script = None
            for path in adt_paths:
                if os.path.exists(path):
                    prepare_script = path
                    break
            
            if not prepare_script:
                return False
            
            cmd = [
                'python2',  # AutoDockTools通常需要Python2
                prepare_script,
                '-r', input_pdb,
                '-o', output_pdbqt,
                '-A', 'hydrogens',  # 添加氢
                '-U', 'nphs',       # 合并非极性氢
                '-v'                # 详细输出
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=120
            )
            
            if result.returncode == 0 and os.path.exists(output_pdbqt):
                if self.logger:
                    self.logger.info(f"ADT转换完成: {output_pdbqt}")
                return True
            else:
                return False
                
        except Exception as e:
            if self.logger:
                self.logger.debug(f"ADT异常: {e}")
            return False
    
    def _prepare_fallback(self, input_pdb: str, output_pdbqt: str) -> bool:
        """简化PDBQT转换 (备用方法)"""
        try:
            with open(input_pdb, 'r') as f:
                pdb_lines = f.readlines()
            
            pdbqt_lines = []
            atom_count = 0
            
            # AutoDock原子类型映射
            atom_type_map = {
                'C': 'C',    # 碳
                'N': 'N',    # 氮
                'O': 'O',    # 氧
                'S': 'S',    # 硫
                'H': 'HD',   # 氢 (给体)
                'P': 'P',    # 磷
                'F': 'F',    # 氟
                'CL': 'Cl',  # 氯
                'BR': 'Br',  # 溴
                'I': 'I',    # 碘
                'FE': 'Fe',  # 铁
                'ZN': 'Zn',  # 锌
                'CA': 'Ca',  # 钙
                'MG': 'Mg',  # 镁
                'NA': 'Na',  # 钠
                'K': 'K',    # 钾
            }
            
            for line in pdb_lines:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom_count += 1
                    
                    # 解析原子信息
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    
                    # 确定AutoDock原子类型
                    element = atom_name[0] if atom_name else 'C'
                    ad_type = atom_type_map.get(element, 'A')  # 默认芳香碳
                    
                    # 构建PDBQT行
                    # PDB格式: ATOM/HETATM + 序号 + 名称 + 残基 + 链 + 序号 + 坐标 + 占据 + B因子
                    # PDBQT添加: 电荷 + 原子类型
                    
                    # 简化: 假设中性电荷
                    charge = 0.0
                    
                    # 构建新行
                    pdbqt_line = line[:66] + f"{charge:6.3f} {ad_type:2s}\n"
                    pdbqt_lines.append(pdbqt_line)
                    
                elif line.startswith('TER') or line.startswith('END'):
                    pdbqt_lines.append(line)
            
            # 写入输出文件
            with open(output_pdbqt, 'w') as f:
                f.writelines(pdbqt_lines)
            
            if self.logger:
                self.logger.info(f"简化PDBQT完成: {output_pdbqt} ({atom_count} atoms)")
            
            return True
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"简化转换失败: {e}")
            return False
    
    def calculate_box(self, pocket_center: List[float], 
                     peptide_diameter: float = None) -> VinaBox:
        """
        计算对接盒子
        
        盒子大小 = 肽直径 + 缓冲空间
        
        Args:
            pocket_center: 口袋中心坐标 [x, y, z]
            peptide_diameter: 环肽直径 (默认使用类变量)
        
        Returns:
            VinaBox对象
        """
        if peptide_diameter is None:
            peptide_diameter = self.PEPTIDE_DIAMETER
        
        # 缓冲空间 (给肽链活动和对接搜索)
        buffer = 10.0  # Å
        
        box_size = peptide_diameter + 2 * buffer
        
        # Vina盒子大小限制 (最大126Å)
        box_size = min(box_size, 120.0)
        
        box = VinaBox(
            center_x=pocket_center[0],
            center_y=pocket_center[1],
            center_z=pocket_center[2],
            size_x=box_size,
            size_y=box_size,
            size_z=box_size
        )
        
        if self.logger:
            self.logger.info(f"对接盒子: center={pocket_center}, size={box_size:.1f}Å")
        
        return box
    
    def generate_config(self, receptor_pdbqt: str, box: VinaBox, 
                       output_config: str,
                       exhaustiveness: int = None,
                       num_modes: int = None) -> bool:
        """
        生成Vina配置文件
        
        Args:
            receptor_pdbqt: 受体PDBQT文件路径
            box: 对接盒子
            output_config: 输出配置文件路径
            exhaustiveness: 搜索详尽程度
            num_modes: 输出构象数
        
        Returns:
            是否成功
        """
        if exhaustiveness is None:
            exhaustiveness = self.config.get('vina', {}).get('exhaustiveness', 4)
        if num_modes is None:
            num_modes = self.config.get('vina', {}).get('num_modes', 9)
        
        config_content = f"""# AutoDock Vina Configuration
# Generated by MCTS Peptide Design Platform

# Receptor
receptor = {receptor_pdbqt}

# Search space
center_x = {box.center_x:.3f}
center_y = {box.center_y:.3f}
center_z = {box.center_z:.3f}
size_x = {box.size_x:.1f}
size_y = {box.size_y:.1f}
size_z = {box.size_z:.1f}

# Search parameters
exhaustiveness = {exhaustiveness}
num_modes = {num_modes}
energy_range = 3

# Output (will be set by docking script)
# out = output.pdbqt
# log = vina.log
"""
        
        try:
            os.makedirs(os.path.dirname(output_config), exist_ok=True)
            with open(output_config, 'w') as f:
                f.write(config_content)
            
            if self.logger:
                self.logger.info(f"Vina配置生成: {output_config}")
            
            return True
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"配置生成失败: {e}")
            return False


def prepare_vina_files(clean_pdb: str, pocket_center: List[float],
                      output_dir: str, config: Dict[str, Any],
                      logger=None) -> Dict[str, str]:
    """
    一站式准备Vina所需的所有文件
    
    Args:
        clean_pdb: 清洁后的PDB文件
        pocket_center: 口袋中心坐标
        output_dir: 输出目录
        config: 配置字典
        logger: 日志对象
    
    Returns:
        文件路径字典 {'receptor': ..., 'config': ..., 'box': ...}
    """
    preparator = VinaPreparator(config, logger)
    
    # 准备受体文件
    receptor_pdbqt = os.path.join(output_dir, 'receptor.pdbqt')
    preparator.prepare_receptor(clean_pdb, receptor_pdbqt)
    
    # 计算对接盒子
    box = preparator.calculate_box(pocket_center)
    
    # 生成配置文件
    config_file = os.path.join(output_dir, 'vina_config.txt')
    preparator.generate_config(receptor_pdbqt, box, config_file)
    
    return {
        'receptor': receptor_pdbqt,
        'config': config_file,
        'box': box
    }
