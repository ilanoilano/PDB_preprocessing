import sys
import os

# 测试导入utils
sys.path.insert(0, '../src')

try:
    from utils import load_config
    print("Import utils OK")
except Exception as e:
    print(f"Import utils failed: {e}")
    import traceback
    traceback.print_exc()
