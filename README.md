
# Niurx的奇怪轮子

> **警告**：使用这些轮子过程中可能会导致一些问题，包括但不限于头晕、恶心、呕吐、高血压、精神失常、心源性休克。请自备医疗支持。

## 🤡 使用须知

`include/det.h`里面是C艹部分的调用，大概改一下就行

`test/test_basic.py`和`test/test_parallel.py`是python部分的调用，分别为多进程和单进程版本。

`det.h`里面设置`openblas`线程数和外面python设置的进程数会影响性能，但我不知道怎么影响的。

似乎进程太多会撞到内存瓶颈。

如果在windows上，需要在python文件开头手动设置dll目录。对此我很绝望。

ubuntu上需要手动去下openblas源码编译并在CMakeLists里面设置好路径。apt的openblas没有lapacke.h，天坑。

macos上需要brew安装openblas，然后不用改CMakeLists。

## 🛠️ 安装

cmake一下，编译器用gcc或者clang，不要用msvc，因为我懒得写兼容了。

上一行由copilot自动生成

然后`pip install ./`
