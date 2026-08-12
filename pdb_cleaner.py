"""
PDB清洗模块
去除结晶水、溶剂，保留金属离子和辅因子，补全缺失残基，pH7.4加氢
"""

import os
import re
from typing import Dict, List, Set, Tuple, Optional
from dataclasses import dataclass


# 水分子残基名
WATER_RESIDUES = {'HOH', 'WAT', 'H2O', 'DOD', 'T3P', 'SPC', 'TIP'}

# 常见溶剂分子 (需要去除)
SOLVENT_RESIDUES = {
    'DMS', 'DMSO',          # 二甲基亚砜
    'EDO', 'EOH', 'ETA',    # 乙醇
    'IPA', 'ISO',           # 异丙醇
    'ACN', 'ACT',           # 乙腈
    'DME', 'DIO',           # 二甲醚/二氧六环
    'GOL', 'GLY',           # 甘油
    'PEG', 'PG4', 'PE4',    # 聚乙二醇
    'SO4', 'SUL',           # 硫酸根
    'PO4', 'PHO',           # 磷酸根
    'CL', 'CLA', 'NA', 'SOD', 'K', 'POT', 'CA', 'CAL',  # 简单离子 (去除)
}

# 保留的金属离子和辅因子
KEEP_RESIDUES = {
    'ZN', 'ZN2',            # 锌
    'MG', 'MG2',            # 镁
    'CA', 'CA2',            # 钙
    'FE', 'FE2', 'FE3',     # 铁
    'MN', 'MN2', 'MN3',     # 锰
    'CU', 'CU1', 'CU2',     # 铜
    'CO', 'CO2',            # 钴
    'NI', 'NI2',            # 镍
    'HEM',                  # 血红素
    'FAD',                  # 黄素腺嘌呤二核苷酸
    'NAD', 'NAP',           # 烟酰胺腺嘌呤二核苷酸
    'ATP', 'ADP', 'AMP',    # 核苷酸
    'FMN',                  # 黄素单核苷酸
}

# 标准氨基酸
STANDARD_AA = {
    'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
    'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYX',
}


@dataclass
class CleanPDBResult:
    """PDB清洗结果"""
    output_file: str
    chain_id: str
    num_atoms: int
    num_residues: int
    kept_hetatms: List[str]  # 保留的HETATM记录
    removed_waters: int
    removed_solvents: int


class PDBCleaner:
    """PDB文件清洗器"""
    
    def __init__(self, logger=None):
        self.logger = logger
    
    def clean(self, input_pdb: str, output_pdb: str, 
              ph: float = 7.4) -> Optional[CleanPDBResult]:
        """
        清洗PDB文件
        
        Args:
            input_pdb: 输入PDB文件路径
            output_pdb: 输出清洗后的PDB文件路径
            ph: pH值 (用于决定质子化状态)
        
        Returns:
            CleanPDBResult对象，失败返回None
        """
        if not os.path.exists(input_pdb):
            if self.logger:
                self.logger.error(f"输入文件不存在: {input_pdb}")
            return None
        
        if self.logger:
            self.logger.info(f"清洗PDB: {input_pdb}")
        
        try:
            # 读取PDB
            with open(input_pdb, 'r') as f:
                lines = f.readlines()
            
            # 分析并清洗
            cleaned_lines, stats = self._process_lines(lines, ph)
            
            # 写入输出文件
            os.makedirs(os.path.dirname(output_pdb), exist_ok=True)
            with open(output_pdb, 'w') as f:
                f.writelines(cleaned_lines)
            
            if self.logger:
                self.logger.info(f"清洗完成: {output_pdb}")
                self.logger.info(f"  去除水分子: {stats['waters']}")
                self.logger.info(f"  去除溶剂: {stats['solvents']}")
                self.logger.info(f"  保留辅因子: {len(stats['kept_hetatms'])}")
            
            # 获取最长链信息
            chain_id = self._get_longest_chain(cleaned_lines)
            num_atoms, num_residues = self._count_atoms_residues(cleaned_lines)
            
            return CleanPDBResult(
                output_file=output_pdb,
                chain_id=chain_id,
                num_atoms=num_atoms,
                num_residues=num_residues,
                kept_hetatms=stats['kept_hetatms'],
                removed_waters=stats['waters'],
                removed_solvents=stats['solvents']
            )
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"PDB清洗失败: {e}")
            return None
    
    def _process_lines(self, lines: List[str], ph: float) -> Tuple[List[str], Dict]:
        """处理PDB行"""
        cleaned = []
        stats = {
            'waters': 0,
            'solvents': 0,
            'kept_hetatms': []
        }
        
        # 记录已处理的残基 (用于去重)
        processed_residues = set()
        
        for line in lines:
            record_type = line[:6].strip()
            
            if record_type == 'ATOM':
                # 标准原子记录，保留
                cleaned.append(line)
                
            elif record_type == 'HETATM':
                # 非标准原子，需要判断
                res_name = line[17:20].strip()
                
                # 去除水分子
                if res_name in WATER_RESIDUES:
                    stats['waters'] += 1
                    continue
                
                # 去除溶剂
                if res_name in SOLVENT_RESIDUES and res_name not in KEEP_RESIDUES:
                    stats['solvents'] += 1
                    continue
                
                # 保留金属离子和辅因子
                if res_name in KEEP_RESIDUES or res_name in STANDARD_AA:
                    cleaned.append(line)
                    res_id = line[17:27]
                    if res_id not in processed_residues:
                        stats['kept_hetatms'].append(res_name)
                        processed_residues.add(res_id)
                else:
                    # 未知分子，保守策略：去除
                    stats['solvents'] += 1
                    if self.logger:
                        self.logger.debug(f"去除未知分子: {res_name}")
                    
            elif record_type in ['CONECT', 'MASTER', 'END', 'ENDMDL']:
                # 保留连接信息和结束标记
                cleaned.append(line)
                
            elif record_type in ['REMARK', 'HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS']:
                # 保留元数据
                cleaned.append(line)
            
            elif record_type == 'SEQRES':
                # 保留序列信息
                cleaned.append(line)
            
            elif record_type == 'CRYST1':
                # 保留晶体信息
                cleaned.append(line)
        
        return cleaned, stats
    
    def _get_longest_chain(self, lines: List[str]) -> str:
        """获取最长的氨基酸链ID"""
        chain_lengths = {}
        
        for line in lines:
            if line.startswith('ATOM'):
                chain_id = line[21].strip()
                res_num = line[22:26].strip()
                
                if chain_id not in chain_lengths:
                    chain_lengths[chain_id] = set()
                chain_lengths[chain_id].add(res_num)
        
        if not chain_lengths:
            return 'A'  # 默认链
        
        # 返回最长的链
        longest_chain = max(chain_lengths.keys(), 
                          key=lambda x: len(chain_lengths[x]))
        return longest_chain
    
    def _count_atoms_residues(self, lines: List[str]) -> Tuple[int, int]:
        """统计原子和残基数量"""
        atoms = 0
        residues = set()
        
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atoms += 1
                res_id = line[17:27]
                residues.add(res_id)
        
        return atoms, len(residues)
    
    def add_hydrogens(self, input_pdb: str, output_pdb: str, ph: float = 7.4) -> bool:
        """
        在指定pH条件下加氢
        
        【需要PDBFixer或OpenBabel】
        占位符含义: 需要外部工具进行质子化
        
        Args:
            input_pdb: 输入PDB文件
            output_pdb: 输出加氢后的PDB文件
            ph: pH值
        
        Returns:
            是否成功
        """
        if self.logger:
            self.logger.info(f"加氢处理 (pH={ph}): {input_pdb}")
        
        # 尝试使用OpenBabel
        try:
            import subprocess
            
            # OpenBabel加氢命令
            cmd = ['obabel', '-i', 'pdb', input_pdb, '-o', 'pdb', '-O', output_pdb, '-h']
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            
            if result.returncode == 0 and os.path.exists(output_pdb):
                if self.logger:
                    self.logger.info(f"OpenBabel加氢完成: {output_pdb}")
                return True
            else:
                if self.logger:
                    self.logger.warning(f"OpenBabel加氢失败: {result.stderr}")
        except Exception as e:
            if self.logger:
                self.logger.warning(f"OpenBabel调用失败: {e}")
        
        # 尝试使用PDBFixer
        try:
            from pdbfixer import PDBFixer
            from simtk.openmm.app import PDBFile
            
            fixer = PDBFixer(filename=input_pdb)
            fixer.addMissingHydrogens(ph)
            
            with open(output_pdb, 'w') as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)
            
            if self.logger:
                self.logger.info(f"PDBFixer加氢完成: {output_pdb}")
            return True
            
        except ImportError:
            if self.logger:
                self.logger.warning("PDBFixer未安装，跳过加氢")
        except Exception as e:
            if self.logger:
                self.logger.warning(f"PDBFixer加氢失败: {e}")
        
        # 如果都失败，复制原文件
        if self.logger:
            self.logger.warning("【ADD_HYDROGEN_FAILED】")
            self.logger.warning("占位符含义: 加氢失败，需要安装OpenBabel或PDBFixer")
            self.logger.warning("命令: sudo apt install openbabel 或 pip install pdbfixer")
        
        import shutil
        shutil.copy(input_pdb, output_pdb)
        return False
    
    def repair_structure(self, input_pdb: str, output_pdb: str) -> bool:
        """
        修复PDB结构 (补全缺失残基和原子)
        
        【需要PDBFixer】
        占位符含义: 需要PDBFixer库修复缺失残基和原子
        
        Args:
            input_pdb: 输入PDB文件
            output_pdb: 输出修复后的PDB文件
        
        Returns:
            是否成功
        """
        if self.logger:
            self.logger.info(f"修复结构: {input_pdb}")
        
        try:
            from pdbfixer import PDBFixer
            from simtk.openmm.app import PDBFile
            
            fixer = PDBFixer(filename=input_pdb)
            
            # 查找并补全缺失残基
            fixer.findMissingResidues()
            if fixer.missingResidues:
                if self.logger:
                    self.logger.info(f"发现缺失残基: {len(fixer.missingResidues)}处")
                fixer.addMissingResidues()
            
            # 查找并补全缺失原子
            fixer.findMissingAtoms()
            if fixer.missingAtoms:
                if self.logger:
                    self.logger.info(f"发现缺失原子: {len(fixer.missingAtoms)}个")
                fixer.addMissingAtoms()
            
            # 写入输出文件
            with open(output_pdb, 'w') as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)
            
            if self.logger:
                self.logger.info(f"结构修复完成: {output_pdb}")
            return True
            
        except ImportError:
            if self.logger:
                self.logger.warning("【PDBFIXER_NOT_INSTALLED】")
                self.logger.warning("占位符含义: PDBFixer未安装，无法修复缺失残基")
                self.logger.warning("安装命令: pip install pdbfixer")
        except Exception as e:
            if self.logger:
                self.logger.error(f"结构修复失败: {e}")
        
        # 复制原文件
        import shutil
        shutil.copy(input_pdb, output_pdb)
        return False


def process_all_pdbs(input_dir: str, output_dir: str, logger=None) -> List[CleanPDBResult]:
    """
    处理目录中的所有PDB文件
    
    Args:
        input_dir: 输入目录
        output_dir: 输出目录
        logger: 日志对象
    
    Returns:
        处理结果列表
    """
    cleaner = PDBCleaner(logger)
    results = []
    
    # 支持的PDB文件扩展名
    pdb_extensions = ['.pdb', '.ent', '.pdb1']
    
    for filename in os.listdir(input_dir):
        if any(filename.endswith(ext) for ext in pdb_extensions):
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, f"{filename}_clean.pdb")
            
            result = cleaner.clean(input_file, output_file)
            if result:
                results.append(result)
    
    return results
