"""测试导入"""
import sys
import os

# 添加项目路径
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(project_root, 'src'))
sys.path.insert(0, project_root)

print(f"Python path: {sys.path[:3]}")

try:
    from utils import load_config, setup_logging
    print("Import OK: load_config, setup_logging")
except Exception as e:
    print(f"Import failed: {e}")

try:
    from fpocket_interface import FpocketInterface
    print("Import OK: FpocketInterface")
except Exception as e:
    print(f"Import failed: {e}")
