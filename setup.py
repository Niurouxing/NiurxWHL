import os
from setuptools import setup, find_packages
import glob
import shutil

# 定义你的包名
package_name = 'mimo'

# 当前目录
here = os.path.abspath(os.path.dirname(__file__))
build_dir = os.path.join(here, 'build')

# 找到所有的.so和.pyd文件
so_files = glob.glob(os.path.join(build_dir, '*.so'))
pyd_files = glob.glob(os.path.join(build_dir, '*.pyd'))

# 合并.so和.pyd文件列表
shared_object_files = so_files + pyd_files

# 如果找不到文件则抛出异常
if not shared_object_files:
    raise FileNotFoundError("No .so or .pyd files found in ./build directory. Please build them first.")

# 创建目标目录路径
target_dir = os.path.join(here, package_name)
os.makedirs(target_dir, exist_ok=True)

# 复制文件到目标目录
for file in shared_object_files:
    shutil.copy(file, target_dir)

# 修改setup函数中的package_data参数，包括.pyd文件
setup(
    name=package_name,
    version="0.0.1",
    author="Dean Moldovan",
    author_email="dean0x7d@gmail.com",
    description="A test project using pybind11 and CMake",
    packages=find_packages(),
    package_data={
        package_name: ['*.so', '*.pyd'],  # 包括.so和.pyd文件
    },
    include_package_data=True,  # 缺省为False，不包括MANIFEST.in中的数据文件
    zip_safe=False,
    python_requires=">=3.7",
)