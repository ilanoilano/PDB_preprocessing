"""
PDB预处理主入口

处理流程:
1. 清洗原始PDB (去水、去溶剂、保留辅因子)
2. 提取口袋 (fpocket检测，输出特征)
3. 准备Vina受体 (PDBQT格式，对接盒子)
4. 准备AlphaFold3输入 (保留辅因子的PDB)

用法:
    python main.py --input_dir /mnt/d/code/my_MCTS_based_PJ/PDB --output_dir /mnt/d/code/my_MCTS_based_PJ/results
"""

import os
import sys
import argparse
import logging
from typing import Dict, List, Any

# 添加项目路径 (当前文件在 preprocessing/ 目录)
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src_path = os.path.join(project_root, 'src')
# 使用append而不是insert，确保系统路径优先
if src_path not in sys.path:
    sys.path.append(src_path)

# 验证yaml可用
try:
    import yaml
except ImportError:
    print("Error: PyYAML not installed. Run: pip3 install pyyaml")
    sys.exit(1)

from utils import load_config, setup_logging
from pdb_cleaner import PDBCleaner, process_all_pdbs
from pocket_extractor import PocketExtractor, extract_pocket_for_esmif
from vina_preparator import prepare_vina_files
from alphafold3_preparator import prepare_alphafold3_inputs


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='PDB Preprocessing Pipeline for MCTS Peptide Design'
    )
    
    parser.add_argument(
        '--input_dir', '-i',
        type=str,
        default='/mnt/d/code/my_MCTS_based_PJ/PDB',
        help='输入PDB文件目录'
    )
    
    parser.add_argument(
        '--output_dir', '-o',
        type=str,
        default='/mnt/d/code/my_MCTS_based_PJ/results',
        help='输出目录'
    )
    
    parser.add_argument(
        '--config', '-c',
        type=str,
        default='../config/config.yaml',
        help='配置文件路径'
    )
    
    parser.add_argument(
        '--steps',
        type=str,
        default='all',
        choices=['all', 'clean', 'pocket', 'vina', 'alphafold3'],
        help='执行哪些步骤'
    )
    
    return parser.parse_args()


def run_preprocessing(input_dir: str, output_dir: str, 
                     config: Dict[str, Any], logger: logging.Logger) -> bool:
    """
    运行完整预处理流程
    
    Args:
        input_dir: 输入PDB文件目录
        output_dir: 输出目录
        config: 配置字典
        logger: 日志对象
    
    Returns:
        是否成功
    """
    logger.info("=" * 80)
    logger.info("PDB预处理流程启动")
    logger.info("=" * 80)
    logger.info(f"输入目录: {input_dir}")
    logger.info(f"输出目录: {output_dir}")
    
    # 检查输入目录
    if not os.path.exists(input_dir):
        logger.error(f"输入目录不存在: {input_dir}")
        return False
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 步骤1: 清洗PDB
    logger.info("\n[步骤1/4] 清洗PDB文件")
    clean_dir = os.path.join(output_dir, '1_cleaned')
    os.makedirs(clean_dir, exist_ok=True)
    
    cleaner = PDBCleaner(logger)
    
    # 处理所有PDB文件
    pdb_files = []
    for f in os.listdir(input_dir):
        if f.endswith(('.pdb', '.ent', '.pdb1')):
            pdb_files.append(f)
    
    if not pdb_files:
        logger.error(f"未找到PDB文件: {input_dir}")
        return False
    
    logger.info(f"发现 {len(pdb_files)} 个PDB文件")
    
    clean_results = []
    for pdb_file in pdb_files:
        input_path = os.path.join(input_dir, pdb_file)
        output_path = os.path.join(clean_dir, f"{os.path.splitext(pdb_file)[0]}_clean.pdb")
        
        result = cleaner.clean(input_path, output_path)
        if result:
            clean_results.append(result)
            logger.info(f"  清洗完成: {pdb_file} -> {result.chain_id}链, "
                       f"{result.num_residues}残基, {result.num_atoms}原子")
    
    if not clean_results:
        logger.error("PDB清洗失败")
        return False
    
    # 使用第一个成功的结果继续处理
    main_result = clean_results[0]
    clean_pdb = main_result.output_file
    
    # 步骤2: 提取口袋
    logger.info("\n[步骤2/4] 提取口袋")
    pocket_dir = os.path.join(output_dir, '2_pockets')
    os.makedirs(pocket_dir, exist_ok=True)
    
    extractor = PocketExtractor(config, logger)
    pockets = extractor.extract_pockets(clean_pdb, pocket_dir)
    
    if not pockets:
        logger.error("口袋提取失败")
        return False
    
    # 选择最佳口袋用于对接
    best_pocket = extractor.get_best_pocket_for_docking(pockets)
    if best_pocket:
        logger.info(f"  最佳口袋: ID={best_pocket.pocket_id}, "
                   f"Score={best_pocket.score:.4f}, "
                   f"Volume={best_pocket.volume:.1f}Å³")
        logger.info(f"  口袋中心: {best_pocket.center}")
        logger.info(f"  残基组成: 疏水{best_pocket.hydrophobic_ratio:.1%}, "
                   f"极性{best_pocket.polar_ratio:.1%}, "
                   f"带电{best_pocket.charged_ratio:.1%}")
        
        # 为ESM-IF准备口袋PDB
        esmif_pocket = os.path.join(pocket_dir, 'pocket_for_esmif.pdb')
        extract_pocket_for_esmif(best_pocket, clean_pdb, esmif_pocket)
        logger.info(f"  ESM-IF口袋PDB: {esmif_pocket}")
    
    # 步骤3: 准备Vina受体
    logger.info("\n[步骤3/4] 准备Vina受体")
    vina_dir = os.path.join(output_dir, '3_vina')
    os.makedirs(vina_dir, exist_ok=True)
    
    if best_pocket:
        vina_files = prepare_vina_files(
            clean_pdb, 
            best_pocket.center,
            vina_dir,
            config,
            logger
        )
        logger.info(f"  受体PDBQT: {vina_files['receptor']}")
        logger.info(f"  Vina配置: {vina_files['config']}")
        logger.info(f"  盒子大小: {vina_files['box'].size_x:.1f}Å")
    
    # 步骤4: 准备AlphaFold3输入
    logger.info("\n[步骤4/4] 准备AlphaFold3输入")
    af3_dir = os.path.join(output_dir, '4_alphafold3')
    os.makedirs(af3_dir, exist_ok=True)
    
    # 复制受体PDB供AlphaFold3使用
    af3_receptor = os.path.join(af3_dir, 'receptor.pdb')
    import shutil
    shutil.copy(clean_pdb, af3_receptor)
    logger.info(f"  AlphaFold3受体: {af3_receptor}")
    
    # 生成示例肽序列输入 (实际使用时从MCTS结果传入)
    sample_peptides = [
        "ACDEFGHICLMNPACQRSTACG",  # 示例序列1
        "ACMASSSCGSSGTCSSSSTCG",   # 示例序列2
    ]
    
    af3_inputs = prepare_alphafold3_inputs(
        clean_pdb,
        sample_peptides,
        af3_dir,
        config,
        logger
    )
    logger.info(f"  AlphaFold3输入: {af3_inputs['count']}个")
    
    # 保存处理汇总
    logger.info("\n" + "=" * 80)
    logger.info("预处理完成!")
    logger.info("=" * 80)
    logger.info(f"输出目录: {output_dir}")
    logger.info(f"  1. 清洁PDB: {clean_dir}")
    logger.info(f"  2. 口袋特征: {pocket_dir}")
    logger.info(f"  3. Vina文件: {vina_dir}")
    logger.info(f"  4. AlphaFold3输入: {af3_dir}")
    logger.info("=" * 80)
    
    return True


def main():
    """主函数"""
    args = parse_args()
    
    # 加载配置
    if not os.path.exists(args.config):
        print(f"错误: 配置文件不存在: {args.config}")
        sys.exit(1)
    
    config = load_config(args.config)
    
    # 设置输出目录
    config['paths']['output_dir'] = args.output_dir
    
    # 设置日志
    logger = setup_logging(config)
    
    # 运行预处理
    success = run_preprocessing(
        args.input_dir,
        args.output_dir,
        config,
        logger
    )
    
    if success:
        logger.info("程序正常结束")
        sys.exit(0)
    else:
        logger.error("程序执行失败")
        sys.exit(1)


if __name__ == '__main__':
    main()
