"""
AlphaFold3输入准备模块
准备终局精修所需的输入文件
"""

import os
import json
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict


@dataclass
class AlphaFold3Input:
    """AlphaFold3输入数据结构"""
    name: str
    sequences: List[Dict]  # 蛋白序列和配体序列
    modelSeeds: List[int]  # 随机种子
    
    def to_json(self) -> str:
        """转换为JSON字符串"""
        return json.dumps(asdict(self), indent=2)


class AlphaFold3Preparator:
    """AlphaFold3输入准备器"""
    
    def __init__(self, config: Dict[str, Any], logger=None):
        self.config = config
        self.logger = logger
        
        # AlphaFold3服务器配置
        self.server_url = config.get('alphafold3', {}).get('server_url', '')
    
    def prepare_input(self, receptor_pdb: str, peptide_sequence: str,
                     output_file: str) -> bool:
        """
        准备AlphaFold3多体预测输入
        
        Args:
            receptor_pdb: 受体蛋白PDB文件 (含辅因子和金属离子)
            peptide_sequence: 肽序列
            output_file: 输出JSON文件路径
        
        Returns:
            是否成功
        """
        if not os.path.exists(receptor_pdb):
            if self.logger:
                self.logger.error(f"受体文件不存在: {receptor_pdb}")
            return False
        
        if self.logger:
            self.logger.info(f"准备AlphaFold3输入: {output_file}")
        
        try:
            # 读取受体序列
            receptor_sequence = self._extract_sequence_from_pdb(receptor_pdb)
            
            if not receptor_sequence:
                if self.logger:
                    self.logger.error("无法提取受体序列")
                return False
            
            # 构建AlphaFold3输入
            af3_input = self._build_input(receptor_sequence, peptide_sequence)
            
            # 写入JSON文件
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, 'w') as f:
                f.write(af3_input.to_json())
            
            if self.logger:
                self.logger.info(f"AlphaFold3输入生成: {output_file}")
            
            return True
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"准备AlphaFold3输入失败: {e}")
            return False
    
    def _extract_sequence_from_pdb(self, pdb_file: str) -> Optional[str]:
        """从PDB文件提取序列"""
        try:
            with open(pdb_file, 'r') as f:
                lines = f.readlines()
            
            # 氨基酸三字母到单字母的映射
            aa_map = {
                'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
                'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
                'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
                'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
            }
            
            # 按链和残基号收集残基
            chain_residues = {}
            
            for line in lines:
                if line.startswith('ATOM'):
                    chain_id = line[21].strip() or 'A'
                    res_num = line[22:26].strip()
                    res_name = line[17:20].strip()
                    
                    if chain_id not in chain_residues:
                        chain_residues[chain_id] = {}
                    
                    # 只记录每个残基一次
                    if res_num not in chain_residues[chain_id]:
                        aa = aa_map.get(res_name, 'X')
                        chain_residues[chain_id][res_num] = aa
            
            # 选择最长的链
            longest_chain = max(chain_residues.keys(), 
                              key=lambda x: len(chain_residues[x]))
            
            # 按残基号排序并构建序列
            sorted_residues = sorted(chain_residues[longest_chain].items(),
                                   key=lambda x: int(x[0]))
            sequence = ''.join([aa for _, aa in sorted_residues])
            
            return sequence
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"提取序列失败: {e}")
            return None
    
    def _build_input(self, receptor_sequence: str, 
                    peptide_sequence: str) -> AlphaFold3Input:
        """构建AlphaFold3输入结构"""
        sequences = [
            {
                "proteinChain": {
                    "sequence": receptor_sequence,
                    "count": 1
                }
            },
            {
                "proteinChain": {
                    "sequence": peptide_sequence,
                    "count": 1
                }
            }
        ]
        
        return AlphaFold3Input(
            name="peptide_docking",
            sequences=sequences,
            modelSeeds=[1]  # 可设置多个种子进行ensemble预测
        )
    
    def prepare_batch_input(self, receptor_pdb: str, 
                           peptide_sequences: List[str],
                           output_dir: str) -> List[str]:
        """
        批量准备AlphaFold3输入
        
        Args:
            receptor_pdb: 受体蛋白PDB文件
            peptide_sequences: 肽序列列表
            output_dir: 输出目录
        
        Returns:
            生成的输入文件路径列表
        """
        os.makedirs(output_dir, exist_ok=True)
        
        input_files = []
        
        for i, seq in enumerate(peptide_sequences):
            output_file = os.path.join(output_dir, f'af3_input_{i+1}.json')
            
            if self.prepare_input(receptor_pdb, seq, output_file):
                input_files.append(output_file)
        
        if self.logger:
            self.logger.info(f"批量准备完成: {len(input_files)}个输入文件")
        
        return input_files
    
    def copy_for_alphafold3(self, clean_pdb: str, output_pdb: str,
                           keep_hetatms: bool = True) -> bool:
        """
        复制清洁PDB供AlphaFold3使用
        
        保留辅因子和金属离子
        
        Args:
            clean_pdb: 清洁后的PDB文件
            output_pdb: 输出PDB文件
            keep_hetatms: 是否保留HETATM记录
        
        Returns:
            是否成功
        """
        if not os.path.exists(clean_pdb):
            return False
        
        try:
            with open(clean_pdb, 'r') as f:
                lines = f.readlines()
            
            output_lines = []
            
            for line in lines:
                if line.startswith('ATOM'):
                    output_lines.append(line)
                elif line.startswith('HETATM'):
                    if keep_hetatms:
                        output_lines.append(line)
                elif line.startswith(('TER', 'END', 'CONECT')):
                    output_lines.append(line)
            
            os.makedirs(os.path.dirname(output_pdb), exist_ok=True)
            with open(output_pdb, 'w') as f:
                f.writelines(output_lines)
            
            if self.logger:
                self.logger.info(f"AlphaFold3输入PDB: {output_pdb}")
            
            return True
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"复制PDB失败: {e}")
            return False


def prepare_alphafold3_inputs(receptor_pdb: str,
                              peptide_sequences: List[str],
                              output_dir: str,
                              config: Dict[str, Any],
                              logger=None) -> Dict[str, Any]:
    """
    一站式准备AlphaFold3所需的所有文件
    
    Args:
        receptor_pdb: 受体蛋白PDB文件
        peptide_sequences: 肽序列列表 (Top-N候选)
        output_dir: 输出目录
        config: 配置字典
        logger: 日志对象
    
    Returns:
        结果字典
    """
    preparator = AlphaFold3Preparator(config, logger)
    
    # 复制受体PDB
    receptor_out = os.path.join(output_dir, 'receptor_for_af3.pdb')
    preparator.copy_for_alphafold3(receptor_pdb, receptor_out)
    
    # 批量生成输入文件
    input_files = preparator.prepare_batch_input(
        receptor_pdb, peptide_sequences, output_dir
    )
    
    return {
        'receptor': receptor_out,
        'inputs': input_files,
        'count': len(input_files)
    }
