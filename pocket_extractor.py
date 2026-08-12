"""
口袋提取模块
利用fpocket提取口袋特征，输出口袋PDB和序列文件
"""

import os
import sys
# 添加项目路径
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(project_root, 'src'))

import subprocess
import re
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, asdict
import json

from fpocket_interface import FpocketInterface, Pocket


@dataclass
class PocketFeatures:
    """口袋特征数据结构"""
    pocket_id: int
    score: float
    drug_score: float
    volume: float
    center: List[float]
    residues: List[str]
    residue_count: int
    hydrophobic_ratio: float  # 疏水残基比例
    polar_ratio: float        # 极性残基比例
    charged_ratio: float      # 带电残基比例
    
    # 几何特征
    diameter: float           # 口袋直径
    depth: float              # 口袋深度
    
    # 文件路径
    pdb_file: str
    sequence_file: str


class PocketExtractor:
    """口袋提取器"""
    
    # 氨基酸性质分类
    HYDROPHOBIC = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'GLY'}
    POLAR = {'SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN'}
    CHARGED_POS = {'LYS', 'ARG', 'HIS'}
    CHARGED_NEG = {'ASP', 'GLU'}
    
    def __init__(self, config: Dict[str, Any], logger=None):
        self.config = config
        self.logger = logger
        self.fpocket = FpocketInterface(config, logger)
    
    def extract_pockets(self, clean_pdb: str, output_dir: str) -> List[PocketFeatures]:
        """
        从清洁PDB中提取所有口袋特征
        
        优先使用fpocket，如果未安装则从配体推断口袋中心
        
        Args:
            clean_pdb: 清洁后的PDB文件路径
            output_dir: 输出目录
        
        Returns:
            口袋特征列表
        """
        if not os.path.exists(clean_pdb):
            if self.logger:
                self.logger.error(f"PDB文件不存在: {clean_pdb}")
            return []
        
        os.makedirs(output_dir, exist_ok=True)
        
        if self.logger:
            self.logger.info(f"提取口袋: {clean_pdb}")
        
        # 检查fpocket是否可用
        if self.fpocket.fpocket_path:
            # 使用fpocket检测口袋
            fpocket_outdir = os.path.join(output_dir, 'fpocket_raw')
            pockets = self.fpocket.get_all_pockets(clean_pdb, fpocket_outdir)
            
            if pockets:
                # 处理每个口袋
                pocket_features = []
                for pocket in pockets:
                    features = self._analyze_pocket(pocket, output_dir)
                    if features:
                        pocket_features.append(features)
                
                # 保存汇总信息
                self._save_pocket_summary(pocket_features, output_dir)
                
                if self.logger:
                    self.logger.info(f"提取完成: {len(pocket_features)}个口袋")
                
                return pocket_features
        
        # fpocket不可用或未检测到口袋，使用配体推断
        if self.logger:
            self.logger.warning("fpocket不可用，使用配体推断口袋中心")
        
        return self._extract_from_ligand(clean_pdb, output_dir)
    
    def _extract_from_ligand(self, clean_pdb: str, output_dir: str) -> List[PocketFeatures]:
        """
        从PDB中的配体推断口袋中心
        
        计算所有HETATM（非水、非金属离子）的质心作为口袋中心
        """
        try:
            with open(clean_pdb, 'r') as f:
                lines = f.readlines()
            
            # 收集配体原子坐标
            ligand_atoms = []
            ligand_residues = set()
            
            # 排除的残基（水、金属离子）
            exclude = {'HOH', 'WAT', 'H2O', 'ZN', 'MG', 'CA', 'FE', 'NA', 'K', 'CL'}
            
            for line in lines:
                if line.startswith('HETATM'):
                    res_name = line[17:20].strip()
                    if res_name in exclude:
                        continue
                    
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        ligand_atoms.append([x, y, z])
                        
                        chain = line[21].strip() or 'A'
                        res_num = line[22:26].strip()
                        ligand_residues.add(f"{res_name}_{chain}_{res_num}")
                    except:
                        continue
            
            if ligand_atoms:
                # 有配体，计算质心
                center = [
                    sum(a[0] for a in ligand_atoms) / len(ligand_atoms),
                    sum(a[1] for a in ligand_atoms) / len(ligand_atoms),
                    sum(a[2] for a in ligand_atoms) / len(ligand_atoms),
                ]
                
                # 估算体积
                volume = len(ligand_atoms) * 20.0
                
                # 创建口袋特征
                pocket = Pocket(
                    id=1,
                    score=0.5,
                    drug_score=0.5,
                    volume=volume,
                    center=center,
                    residues=list(ligand_residues),
                    pdb_file=clean_pdb
                )
                
                features = self._analyze_pocket(pocket, output_dir)
                
                if features:
                    if self.logger:
                        self.logger.info(f"从配体推断口袋中心: {center}")
                        self.logger.info(f"配体残基数: {len(ligand_residues)}")
                    return [features]
            
            # 没有配体，提取最长链并用fpocket
            if self.logger:
                self.logger.warning("PDB中未找到配体，提取最长链进行fpocket分析")
            
            return self._extract_from_longest_chain(clean_pdb, output_dir)
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"配体推断失败: {e}")
            return []
    
    def _extract_from_longest_chain(self, clean_pdb: str, output_dir: str) -> List[PocketFeatures]:
        """
        提取最长氨基酸链，用fpocket分析潜在口袋
        """
        try:
            with open(clean_pdb, 'r') as f:
                lines = f.readlines()
            
            # 统计每条链的残基数
            chain_residues = {}
            for line in lines:
                if line.startswith('ATOM'):
                    chain_id = line[21].strip() or 'A'
                    res_num = line[22:26].strip()
                    
                    if chain_id not in chain_residues:
                        chain_residues[chain_id] = set()
                    chain_residues[chain_id].add(res_num)
            
            if not chain_residues:
                if self.logger:
                    self.logger.error("PDB中没有ATOM记录")
                return []
            
            # 找到最长链
            longest_chain = max(chain_residues.keys(), 
                              key=lambda x: len(chain_residues[x]))
            chain_length = len(chain_residues[longest_chain])
            
            if self.logger:
                self.logger.info(f"最长链: {longest_chain}, 残基数: {chain_length}")
            
            # 提取最长链的PDB
            chain_pdb = os.path.join(output_dir, f'chain_{longest_chain}.pdb')
            with open(chain_pdb, 'w') as f:
                for line in lines:
                    if line.startswith('ATOM'):
                        chain_id = line[21].strip() or 'A'
                        if chain_id == longest_chain:
                            f.write(line)
                    elif line.startswith('HETATM'):
                        # 保留该链的辅因子
                        chain_id = line[21].strip() or 'A'
                        if chain_id == longest_chain:
                            f.write(line)
                f.write("END\n")
            
            if self.logger:
                self.logger.info(f"最长链PDB: {chain_pdb}")
            
            # 用fpocket分析最长链
            if self.fpocket.fpocket_path:
                fpocket_outdir = os.path.join(output_dir, f'fpocket_chain_{longest_chain}')
                pockets = self.fpocket.get_all_pockets(chain_pdb, fpocket_outdir)
                
                if pockets:
                    pocket_features = []
                    for pocket in pockets:
                        features = self._analyze_pocket(pocket, output_dir)
                        if features:
                            pocket_features.append(features)
                    
                    if pocket_features:
                        self._save_pocket_summary(pocket_features, output_dir)
                        if self.logger:
                            self.logger.info(f"从最长链提取 {len(pocket_features)} 个口袋")
                        return pocket_features
            
            # fpocket不可用或未找到口袋，使用几何中心
            if self.logger:
                self.logger.warning("fpocket未找到口袋，使用链几何中心")
            
            # 计算链的几何中心
            chain_atoms = []
            for line in lines:
                if line.startswith('ATOM'):
                    chain_id = line[21].strip() or 'A'
                    if chain_id == longest_chain:
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            chain_atoms.append([x, y, z])
                        except:
                            continue
            
            if chain_atoms:
                center = [
                    sum(a[0] for a in chain_atoms) / len(chain_atoms),
                    sum(a[1] for a in chain_atoms) / len(chain_atoms),
                    sum(a[2] for a in chain_atoms) / len(chain_atoms),
                ]
                
                # 创建虚拟口袋
                pocket = Pocket(
                    id=1,
                    score=0.3,
                    drug_score=0.3,
                    volume=chain_length * 100,  # 粗略估计
                    center=center,
                    residues=[f"{longest_chain}_{i}" for i in range(min(10, chain_length))],
                    pdb_file=chain_pdb
                )
                
                features = self._analyze_pocket(pocket, output_dir)
                if features:
                    if self.logger:
                        self.logger.info(f"使用链几何中心: {center}")
                    return [features]
            
            return []
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"最长链提取失败: {e}")
            return []
    
    def _analyze_pocket(self, pocket: Pocket, output_dir: str) -> Optional[PocketFeatures]:
        """分析单个口袋的特征"""
        try:
            # 计算残基组成
            hydrophobic_count = 0
            polar_count = 0
            charged_count = 0
            
            for res in pocket.residues:
                res_name = res.split('_')[0]
                if res_name in self.HYDROPHOBIC:
                    hydrophobic_count += 1
                elif res_name in self.POLAR:
                    polar_count += 1
                elif res_name in self.CHARGED_POS or res_name in self.CHARGED_NEG:
                    charged_count += 1
            
            total = len(pocket.residues)
            if total > 0:
                hydrophobic_ratio = hydrophobic_count / total
                polar_ratio = polar_count / total
                charged_ratio = charged_count / total
            else:
                hydrophobic_ratio = polar_ratio = charged_ratio = 0
            
            # 估算口袋直径 (基于体积)
            # V = 4/3 * π * r³ → r = (3V/4π)^(1/3)
            import math
            radius = ((3 * pocket.volume) / (4 * math.pi)) ** (1/3)
            diameter = 2 * radius
            
            # 估算深度 (简化: 基于残基数量)
            depth = total * 3.5  # 粗略估计，每个残基约3.5Å
            
            # 复制口袋PDB到输出目录
            pocket_pdb_out = os.path.join(output_dir, f'pocket_{pocket.id}.pdb')
            if os.path.exists(pocket.pdb_file):
                import shutil
                shutil.copy(pocket.pdb_file, pocket_pdb_out)
            
            # 生成序列文件
            sequence_file = os.path.join(output_dir, f'pocket_{pocket.id}_sequence.txt')
            self._write_sequence_file(pocket, sequence_file)
            
            features = PocketFeatures(
                pocket_id=pocket.id,
                score=pocket.score,
                drug_score=pocket.drug_score,
                volume=pocket.volume,
                center=pocket.center,
                residues=pocket.residues,
                residue_count=len(pocket.residues),
                hydrophobic_ratio=hydrophobic_ratio,
                polar_ratio=polar_ratio,
                charged_ratio=charged_ratio,
                diameter=diameter,
                depth=depth,
                pdb_file=pocket_pdb_out,
                sequence_file=sequence_file
            )
            
            return features
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"分析口袋{pocket.id}失败: {e}")
            return None
    
    def _write_sequence_file(self, pocket: Pocket, output_file: str):
        """写入口袋序列文件"""
        try:
            with open(output_file, 'w') as f:
                f.write(f">Pocket_{pocket.id}\n")
                f.write(f"Score: {pocket.score:.4f}\n")
                f.write(f"Drug_Score: {pocket.drug_score:.4f}\n")
                f.write(f"Volume: {pocket.volume:.2f} A^3\n")
                f.write(f"Center: {pocket.center[0]:.3f}, {pocket.center[1]:.3f}, {pocket.center[2]:.3f}\n")
                f.write(f"Residues ({len(pocket.residues)}):\n")
                for res in pocket.residues:
                    f.write(f"  {res}\n")
        except Exception as e:
            if self.logger:
                self.logger.error(f"写入序列文件失败: {e}")
    
    def _save_pocket_summary(self, pockets: List[PocketFeatures], output_dir: str):
        """保存口袋汇总信息"""
        summary_file = os.path.join(output_dir, 'pocket_summary.json')
        
        summary = {
            'total_pockets': len(pockets),
            'pockets': [asdict(p) for p in pockets]
        }
        
        try:
            with open(summary_file, 'w') as f:
                json.dump(summary, f, indent=2)
            
            if self.logger:
                self.logger.info(f"口袋汇总保存: {summary_file}")
        except Exception as e:
            if self.logger:
                self.logger.error(f"保存汇总失败: {e}")
    
    def get_best_pocket_for_docking(self, pockets: List[PocketFeatures]) -> Optional[PocketFeatures]:
        """
        选择最适合对接的口袋
        
        选择标准:
        1. 评分高
        2. 体积适中 (100-2000 Å³)
        3. 成药性好 (drug_score > 0.5)
        
        Args:
            pockets: 口袋列表
        
        Returns:
            最佳口袋
        """
        if not pockets:
            return None
        
        # 过滤并排序
        valid_pockets = [
            p for p in pockets 
            if 100 <= p.volume <= 2000 and p.drug_score > 0.3
        ]
        
        if not valid_pockets:
            # 如果没有完全符合的，返回评分最高的
            return max(pockets, key=lambda x: x.score)
        
        # 按综合评分排序 (fpocket评分 + 成药性)
        best = max(valid_pockets, key=lambda x: x.score + x.drug_score)
        
        return best


def extract_pocket_for_esmif(pocket_features: PocketFeatures, 
                              clean_pdb: str, 
                              output_file: str) -> bool:
    """
    为ESM-IF准备口袋PDB文件
    
    从清洁PDB中提取口袋区域的残基，保存为新的PDB文件
    
    Args:
        pocket_features: 口袋特征
        clean_pdb: 清洁后的PDB文件
        output_file: 输出PDB文件
    
    Returns:
        是否成功
    """
    if not os.path.exists(clean_pdb):
        return False
    
    try:
        # 读取清洁PDB
        with open(clean_pdb, 'r') as f:
            lines = f.readlines()
        
        # 解析口袋残基ID
        pocket_residue_ids = set()
        for res in pocket_features.residues:
            parts = res.split('_')
            if len(parts) >= 3:
                res_name, chain, res_num = parts[0], parts[1], parts[2]
                pocket_residue_ids.add((chain, res_num))
        
        # 提取口袋残基的原子
        pocket_lines = []
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain = line[21].strip()
                res_num = line[22:26].strip()
                
                if (chain, res_num) in pocket_residue_ids:
                    pocket_lines.append(line)
        
        # 写入输出文件
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with open(output_file, 'w') as f:
            f.writelines(pocket_lines)
            f.write("END\n")
        
        return True
        
    except Exception as e:
        return False
