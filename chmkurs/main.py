import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def read_points(filename):
    with open(filename, 'r') as file:
        num_points = int(file.readline())
        points = []
        for _ in range(num_points):
            coords = list(map(float, file.readline().split()))
            if len(coords) == 3:
                points.append(coords)
    return np.array(points)

def plot_selected_region(points, x_lim=None, y_lim=None, z_lim=None):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # ?????????? ????? ?? ???????? ????????
    mask = np.ones(len(points), dtype=bool)
    if x_lim:
        mask &= (points[:,0] >= x_lim[0]) & (points[:,0] <= x_lim[1])
    if y_lim:
        mask &= (points[:,1] >= y_lim[0]) & (points[:,1] <= y_lim[1])
    if z_lim:
        mask &= (points[:,2] >= z_lim[0]) & (points[:,2] <= z_lim[1])
    
    filtered_points = points[mask]
    
    # ????????????
    sc = ax.scatter(filtered_points[:,0], filtered_points[:,1], filtered_points[:,2],
                   c=filtered_points[:,2], cmap='viridis', marker='o', s=50, alpha=0.8)
    
    # ????????? ????
    ax.set_xlabel('X ')
    ax.set_ylabel('Y ')
    ax.set_zlabel('Z ')
    
    # ????????? ???????? ????
    if x_lim: ax.set_xlim(x_lim)
    if y_lim: ax.set_ylim(y_lim)
    if z_lim: ax.set_zlim(z_lim)
    
    # ???????? ?????
    cbar = fig.colorbar(sc, ax=ax, shrink=0.7)
    cbar.set_label('Ось (Z)')
    
    plt.title('Визуализация сетки')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    filename = "cross.txt"
    points = read_points(filename)
    
    # ??????: ???????? ??????? 0.5 ? x ? 1.5, 0 ? y ? 1, ??? z
    plot_selected_region(points,
                       x_lim=(-150.0, 150.0),
                       y_lim=(-150.0, 150.0),
                       z_lim=None)  # None ???????? "??? ???????????"