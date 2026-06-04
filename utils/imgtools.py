import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

class PlotOverlay:
    """
    用于在同一张图上叠加不同数据集的类 (简洁版)
    """

    def __init__(self, figsize=(8, 6), xlim=(0.1, 450), ylim=(1, 100),
                 xlabel='Shear Rate (1/s)', ylabel='$\\eta$ (Pa.s)', ft=16,
                 tick_width=3, is_xlog=True, is_ylog=True):
        """
        初始化绘图。
        """
        self.fig, self.ax = plt.subplots(figsize=figsize)
        if is_xlog:
            self.ax.set_xscale('log')
        if is_ylog:
            self.ax.set_yscale('log')
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        self.ax.set_xlabel(xlabel, fontsize=ft)
        self.ax.set_ylabel(ylabel, fontsize=ft)
        self.ax.tick_params(axis='both', direction='in', width=tick_width)

    def add_data(self, x, y, marker=None, color='black', linestyle='-', label=None, linewidth=1.5):
        """
        将数据添加到绘图中。
        """
        self.ax.plot(x, y, marker=marker, color=color, linestyle=linestyle, label=label, linewidth=linewidth)
        if label:
            self.ax.legend()

    def show(self):
        """
        显示绘图。
        """
        plt.show()

# 示例用法
if __name__ == '__main__':

    # 创建一些示例数据
    x1 = np.logspace(-1, 2, 50)
    y1 = 10 * np.exp(-x1 / 10) + 1

    x2 = np.logspace(-1, 2, 50)
    y2 = 5 * np.exp(-x2 / 20) + 0.5

    # 创建 PlotOverlay 对象
    plotter = PlotOverlay(ylabel='$\\eta$ (Pa.s)')

    # 添加第一组数据
    plotter.add_data(x1, y1, marker='o', color='blue', linestyle='-', label='Data Set 1')

    # 添加第二组数据
    plotter.add_data(x2, y2, marker='s', color='red', linestyle='--', label='Data Set 2')

    # 显示绘图
    plotter.show()
