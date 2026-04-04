import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def read_component(filename):
    """
    Читает файл формата:
        x  y  displacement_component
    """
    data = np.loadtxt(filename)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 3:
        raise ValueError(f"Файл {filename} должен содержать минимум 3 столбца: x y value")

    xy = data[:, :2]
    values = data[:, 2]
    return xy, values


def point_key(x, y, ndigits=12):
    return (round(float(x), ndigits), round(float(y), ndigits))


def build_body_mesh(file_x, file_y):
    """
    Собирает регулярную сетку тела из двух файлов:
      - file_x: x y u_x
      - file_y: x y u_y

    Возвращает:
      X0, Y0  - исходные координаты узлов в виде матриц
      Ux, Uy  - смещения узлов в виде матриц
    """
    xy_x, ux_vals = read_component(file_x)
    xy_y, uy_vals = read_component(file_y)

    ux_map = {point_key(x, y): ux for (x, y), ux in zip(xy_x, ux_vals)}
    uy_map = {point_key(x, y): uy for (x, y), uy in zip(xy_y, uy_vals)}

    all_keys = sorted(set(ux_map.keys()) | set(uy_map.keys()), key=lambda p: (p[1], p[0]))
    if not all_keys:
        raise ValueError("Не найдено ни одного узла.")

    xs = np.array(sorted({p[0] for p in all_keys}), dtype=float)
    ys = np.array(sorted({p[1] for p in all_keys}), dtype=float)

    X0, Y0 = np.meshgrid(xs, ys)
    Ux = np.zeros_like(X0, dtype=float)
    Uy = np.zeros_like(Y0, dtype=float)

    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            key = point_key(x, y)

            if key not in ux_map:
                raise ValueError(f"В файле {file_x} отсутствует узел ({x}, {y})")
            if key not in uy_map:
                raise ValueError(f"В файле {file_y} отсутствует узел ({x}, {y})")

            Ux[j, i] = ux_map[key]
            Uy[j, i] = uy_map[key]

    return X0, Y0, Ux, Uy


def plot_grid(ax, X, Y, color, linewidth=1.2, alpha=1.0, linestyle='-'):
    """
    Рисует сетку по матрицам координат X, Y:
    - линии по строкам
    - линии по столбцам
    """
    nrows, ncols = X.shape

    for j in range(nrows):
        ax.plot(X[j, :], Y[j, :],
                color=color, linewidth=linewidth, alpha=alpha, linestyle=linestyle)

    for i in range(ncols):
        ax.plot(X[:, i], Y[:, i],
                color=color, linewidth=linewidth, alpha=alpha, linestyle=linestyle)


def main():
    parser = argparse.ArgumentParser(
        description="Построение исходной и деформированной сетки двух тел по файлам смещений."
    )

    # Пути к файлам
    parser.add_argument("--bottom-x", default="bottom_displacement_x.txt",
                        help="Файл смещений u_x для нижнего тела")
    parser.add_argument("--bottom-y", default="bottom_displacement_y.txt",
                        help="Файл смещений u_y для нижнего тела")
    parser.add_argument("--top-x", default="top_displacement_x.txt",
                        help="Файл смещений u_x для верхнего тела")
    parser.add_argument("--top-y", default="top_displacement_y.txt",
                        help="Файл смещений u_y для верхнего тела")

    # Общие цвета
    parser.add_argument("--before-color", default="0.65",
                        help="Цвет исходной конфигурации (по умолчанию серый)")
    parser.add_argument("--after-color", default="tab:red",
                        help="Цвет деформированной конфигурации (по умолчанию красный)")

    # При необходимости — отдельные цвета для каждого тела
    parser.add_argument("--bottom-before-color", default=None,
                        help="Цвет исходной конфигурации нижнего тела")
    parser.add_argument("--bottom-after-color", default=None,
                        help="Цвет деформированной конфигурации нижнего тела")
    parser.add_argument("--top-before-color", default=None,
                        help="Цвет исходной конфигурации верхнего тела")
    parser.add_argument("--top-after-color", default=None,
                        help="Цвет деформированной конфигурации верхнего тела")

    # Прочие параметры
    parser.add_argument("--scale", type=float, default=1.0,
                        help="Масштаб смещений")
    parser.add_argument("--lw-before", type=float, default=1.0,
                        help="Толщина линий исходной конфигурации")
    parser.add_argument("--lw-after", type=float, default=1.4,
                        help="Толщина линий деформированной конфигурации")
    parser.add_argument("--alpha-before", type=float, default=0.9,
                        help="Прозрачность исходной конфигурации")
    parser.add_argument("--alpha-after", type=float, default=1.0,
                        help="Прозрачность деформированной конфигурации")
    parser.add_argument("--figsize", type=float, nargs=2, default=(8.0, 6.0),
                        metavar=("W", "H"),
                        help="Размер фигуры")
    parser.add_argument("--title", default="",
                        help="Заголовок графика")
    parser.add_argument("--output", default=None,
                        help="Если указан путь, рисунок сохраняется в файл")
    parser.add_argument("--dpi", type=int, default=200,
                        help="DPI при сохранении")
    parser.add_argument("--show-legend", action="store_true",
                        help="Показать легенду")

    args = parser.parse_args()

    # Цвета с учетом возможных переопределений
    bottom_before_color = args.bottom_before_color or args.before_color
    bottom_after_color = args.bottom_after_color or args.after_color
    top_before_color = args.top_before_color or args.before_color
    top_after_color = args.top_after_color or args.after_color

    # Чтение данных
    Xb0, Yb0, Ubx, Uby = build_body_mesh(args.bottom_x, args.bottom_y)
    Xt0, Yt0, Utx, Uty = build_body_mesh(args.top_x, args.top_y)

    # Деформированные координаты
    Xb = Xb0 + args.scale * Ubx
    Yb = Yb0 + args.scale * Uby

    Xt = Xt0 + args.scale * Utx
    Yt = Yt0 + args.scale * Uty

    # Построение
    fig, ax = plt.subplots(figsize=args.figsize)

    # Исходные сетки
    plot_grid(ax, Xb0, Yb0,
              color=bottom_before_color,
              linewidth=args.lw_before,
              alpha=args.alpha_before)

    plot_grid(ax, Xt0, Yt0,
              color=top_before_color,
              linewidth=args.lw_before,
              alpha=args.alpha_before)

    # Деформированные сетки
    plot_grid(ax, Xb, Yb,
              color=bottom_after_color,
              linewidth=args.lw_after,
              alpha=args.alpha_after)

    plot_grid(ax, Xt, Yt,
              color=top_after_color,
              linewidth=args.lw_after,
              alpha=args.alpha_after)

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(args.title)
    ax.grid(False)

    # Автоподбор пределов с небольшим отступом
    all_x = np.concatenate([Xb0.ravel(), Xt0.ravel(), Xb.ravel(), Xt.ravel()])
    all_y = np.concatenate([Yb0.ravel(), Yt0.ravel(), Yb.ravel(), Yt.ravel()])

    dx = all_x.max() - all_x.min()
    dy = all_y.max() - all_y.min()
    pad_x = 0.05 * dx if dx > 0 else 0.1
    pad_y = 0.05 * dy if dy > 0 else 0.1

    ax.set_xlim(all_x.min() - pad_x, all_x.max() + pad_x)
    ax.set_ylim(all_y.min() - pad_y, all_y.max() + pad_y)

    if args.show_legend:
        handles = [
            Line2D([0], [0], color=bottom_before_color, lw=args.lw_before, label="Bottom body: initial"),
            Line2D([0], [0], color=bottom_after_color,  lw=args.lw_after,  label="Bottom body: deformed"),
            Line2D([0], [0], color=top_before_color,    lw=args.lw_before, label="Top body: initial"),
            Line2D([0], [0], color=top_after_color,     lw=args.lw_after,  label="Top body: deformed"),
        ]
        ax.legend(handles=handles)

    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=args.dpi, bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
    main()