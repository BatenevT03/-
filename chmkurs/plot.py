import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np
from matplotlib.patches import Patch

def read_points(filename, read_values=False, y_zero=False):
    """Чтение точек из файла"""
    points = []
    values = []
    with open(filename, 'r') as file:
        for line in file:
            if line.strip():
                parts = list(map(float, line.strip().split()))
                if y_zero:
                    if len(parts) >= 2:
                        points.append([parts[0], 0.0, parts[1]])
                        if read_values and len(parts) > 2:
                            values.append(parts[2])
                else:
                    points.append(parts[:3])
                    if read_values and len(parts) > 3:
                        values.append(parts[3])
    return np.array(points), (np.array(values) if read_values and values else None)

def generate_quads(points):
    """Генерация четырёхугольников"""
    quads = []
    n = len(points)
    i = 0
    while i < n:
        p1 = points[i]
        if i + 1 >= n or not np.isclose(points[i+1][1], p1[1]):
            i += 1
            continue
        p2 = points[i+1]
        k = -1
        for j in range(i+1, n):
            if np.isclose(points[j][0], p2[0]) and not np.isclose(points[j][1], p2[1]):
                k = j
                break
        if k == -1 or k + 1 > n or not np.isclose(points[k][1], points[k-1][1]):
            i += 1
            continue
        p4 = points[k-1]
        p3 = points[k]
        quads.append([p1, p2, p3, p4])
        i += 1
    return quads

def draw_sphere(ax, x, y, z, radius, color, alpha=1.0):
    """Рисование сферы"""
    u = np.linspace(0, 2 * np.pi, 30)
    v = np.linspace(0, np.pi, 30)
    x_sphere = x + radius * np.outer(np.cos(u), np.sin(v))
    y_sphere = y + radius * np.outer(np.sin(u), np.sin(v))
    z_sphere = z + radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x_sphere, y_sphere, z_sphere, 
                   color=color, alpha=alpha, shade=True, edgecolor='none')

def plot_all_elements(quads, line_points, receiver_points, sphere_radius=15):
    """Отрисовка всех элементов без осей и фона"""
    fig = plt.figure(figsize=(16, 12), facecolor='white')
    ax = fig.add_subplot(111, projection='3d')
    
    # Полностью убираем фоновые плоскости
    ax.xaxis.pane.set_visible(False)
    ax.yaxis.pane.set_visible(False)
    ax.zaxis.pane.set_visible(False)
    
    # Убираем оси (линии X,Y,Z)
    ax.set_axis_off()
    
    # Сбор всех точек для правильного масштабирования
    all_points = []
    if quads:
        all_points.extend(np.vstack(quads))
    if line_points is not None:
        all_points.extend(line_points)
    if receiver_points is not None:
        all_points.extend(receiver_points)
    
    if all_points:
        all_points = np.array(all_points)
        x_min, x_max = all_points[:, 0].min(), all_points[:, 0].max()
        y_min, y_max = all_points[:, 1].min(), all_points[:, 1].max()
        z_min, z_max = all_points[:, 2].min(), all_points[:, 2].max()
        
        # Добавляем 10% отступ
        padding = 0.1
        x_pad = (x_max - x_min) * padding
        y_pad = (y_max - y_min) * padding
        z_pad = (z_max - z_min) * padding
        
        ax.set_xlim(x_min - x_pad, x_max + x_pad)
        ax.set_ylim(y_min - y_pad, y_max + y_pad)
        ax.set_zlim(z_min - z_pad, z_max + z_pad)
    
    # Отрисовка элементов
    if quads:
        quad_lines = [quad + [quad[0]] for quad in quads]
        ax.add_collection3d(Line3DCollection(quad_lines, colors='black',
                                          linewidths=1.5, linestyle='-', alpha=0.7))
    
    if line_points is not None and len(line_points) >= 2:
        draw_sphere(ax, *line_points[0], sphere_radius, 'blue')
        draw_sphere(ax, *line_points[1], sphere_radius, 'red')
    
    if receiver_points is not None:
        for point in receiver_points:
            draw_sphere(ax, *point, sphere_radius, 'green', alpha=0.7)
    
    # Легенда (если нужна)
    legend_elements = []
   
    
    if legend_elements:
        ax.legend(handles=legend_elements, fontsize=12, loc='upper right')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Загрузка данных
    points, _ = read_points("hill.txt")
    quads = generate_quads(points)
    
    try:
        line_points, _ = read_points("el.txt")
        if len(line_points) < 2:
            print("В el.txt должно быть минимум 2 точки")
            line_points = None
    except FileNotFoundError:
        print("Файл el.txt не найден")
        line_points = None
    
    try:
        receiver_points, _ = read_points("value_from_receiver.txt", y_zero=True)
    except FileNotFoundError:
        print("Файл value_from_receiver.txt не найден")
        receiver_points = None
    
    plot_all_elements(quads, line_points, receiver_points, sphere_radius=7)