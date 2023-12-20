import os
from setuptools import setup, find_packages
import glob
import shutil

# 定义你的包名
package_name = 'mimo'

# 当前目录
here = os.path.abspath(os.path.dirname(__file__))
build_dir = os.path.join(here, 'build')

# 找到所有的.so文件
so_files = glob.glob(os.path.join(build_dir, '*.so'))

# 如果找不到.so文件则抛出异常
if not so_files:
    raise FileNotFoundError("No .so files found in ./build directory. Please build them first.")

# 创建目标目录路径
target_dir = os.path.join(here, package_name)
os.makedirs(target_dir, exist_ok=True)

# 复制.so文件到目标目录
for so_file in so_files:
    shutil.copy(so_file, target_dir)

setup(
    name=package_name,
    version="0.0.1",
    author="Dean Moldovan",
    author_email="dean0x7d@gmail.com",
    description="A test project using pybind11 and CMake",
    packages=find_packages(),
    package_data={
        package_name: ['*.so'],  # 包括.so文件
    },
    include_package_data=True,  # 缺省为False，不包括MANIFEST.in中的数据文件
    zip_safe=False,
    python_requires=">=3.7",
)