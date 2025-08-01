import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np

def read_points(filename):
    """Чтение точек из файла"""
    points = []
    with open(filename, 'r') as file:
        for line in file:
            if line.strip():
                x, y, z = map(float, line.strip().split())
                points.append([x, y, z])
    return np.array(points)

def generate_quads(points):
    """Генерация четырёхугольников по строгому алгоритму"""
    quads = []
    n = len(points)
    
    i = 0
    while i < n:
        p1 = points[i]
        
        # Ищем вторую точку (i+1 с таким же y)
        if i + 1 >= n or not np.isclose(points[i+1][1], p1[1]):
            i += 1
            continue
        p2 = points[i+1]
        
        # Ищем четвертую точку (с x как у p2 и другим y)
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

def draw_sphere(ax, x, y, z, radius, color):
    """Рисует идеальную сферу"""
    u = np.linspace(0, 2 * np.pi, 30)
    v = np.linspace(0, np.pi, 30)
    
    x_sphere = x + radius * np.outer(np.cos(u), np.sin(v))
    y_sphere = y + radius * np.outer(np.sin(u), np.sin(v))
    z_sphere = z + radius * np.outer(np.ones(np.size(u)), np.cos(v))
    
    ax.plot_surface(x_sphere, y_sphere, z_sphere, 
                   color=color, alpha=1.0, shade=True, edgecolor='none')

def plot_quads_and_points(quads, points, sphere_radius=15):
    """Отрисовка без фоновых плоскостей"""
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Убираем полупрозрачные плоскости
    #ax.xaxis.pane.fill = False
    #ax.yaxis.pane.fill = False
    #ax.zaxis.pane.fill = False
    
    # Убираем серые линии сетки
    ax.xaxis._axinfo["grid"].update({"visible": False})
    ax.yaxis._axinfo["grid"].update({"visible": False})
    ax.zaxis._axinfo["grid"].update({"visible": False})
    

    # 1. Отрисовка четырёхугольников
    quad_lines = []
    for quad in quads:
        closed_line = quad + [quad[0]]
        quad_lines.append(closed_line)
    
    ax.add_collection3d(Line3DCollection(quad_lines, 
                                      colors='black',
                                      linewidths=1.5,
                                      linestyle='-',
                                      alpha=0.7))
    
    # 2. Отрисовка сфер
    if points is not None and len(points) >= 2:
        draw_sphere(ax, points[0][0], points[0][1], points[0][2], 
                   sphere_radius, 'blue')
        draw_sphere(ax, points[1][0], points[1][1], points[1][2], 
                   sphere_radius, 'red')
        
        # Легенда
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='blue', label='Начальная точка'),
            Patch(facecolor='red', label='Конечная точка')
        ]
        ax.legend(handles=legend_elements, fontsize=12)

    # Настройка осей
    ax.set_xlim(-300, 300)
    ax.set_ylim(-300, 300)

    
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_zlabel('Z', fontsize=12)
    ax.set_title('Визуализация без фоновых плоскостей', fontsize=14, pad=20)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    points = read_points("hill.txt")
    quads = generate_quads(points)
    
    try:
        sphere_points = read_points("el.txt")
        if len(sphere_points) < 2:
            print("Ошибка: в el.txt нужно минимум 2 точки")
            sphere_points = None
    except FileNotFoundError:
        print("Файл el.txt не найден")
        sphere_points = None
    
    SPHERE_RADIUS = 20
    plot_quads_and_points(quads, 
                         sphere_points[:2] if sphere_points is not None else None,
                         sphere_radius=SPHERE_RADIUS)