import argparse
import json
import gzip
import sys
import numpy as np
import open3d as o3d
import open3d.visualization.gui as gui
import open3d.visualization.rendering as rendering

def visualize_ramachandran(filepath, category, x_angle, y_angle, z_angle):
    points = []
    print(f"Parsing '{category}' category from {filepath}...")
    print(f"Mapping axes -> X: {x_angle}, Y: {y_angle}, Z: {z_angle}")
    
    # Automatically handle both standard and .gz compressed files
    open_func = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'

    try:
        with open_func(filepath, mode) as f:
            for line in f:
                if not line.strip():
                    continue
                    
                data = json.loads(line)
                
                # Filter by category and ensure the requested angles exist
                if (data.get('rama_category') == category and
                    data.get(x_angle) is not None and
                    data.get(y_angle) is not None and
                    data.get(z_angle) is not None):
                    
                    points.append([data[x_angle], data[y_angle], data[z_angle]])
    except FileNotFoundError:
        print(f"Error: Could not find file '{filepath}'.")
        sys.exit(1)

    if not points:
        print(f"No valid points found for category '{category}' with angles {x_angle}, {y_angle}, {z_angle}.")
        sys.exit(0)

    print(f"Successfully extracted {len(points)} points. Building modern UI...")

    # 1. Convert to Open3D format
    points_np = np.array(points)
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(points_np)
    
    # EXPLICIT COLOR: Paint all points a visible deep blue
    pcd.paint_uniform_color([0.1, 0.4, 0.8])

    # 2. Build a strictly demarcated -180 to 180 bounding box
    box_points = [
        [-180, -180, -180], [180, -180, -180], [-180, 180, -180], [180, 180, -180],
        [-180, -180, 180], [180, -180, 180], [-180, 180, 180], [180, 180, 180]
    ]
    box_lines = [
        [0, 1], [0, 2], [1, 3], [2, 3],
        [4, 5], [4, 6], [5, 7], [6, 7],
        [0, 4], [1, 5], [2, 6], [3, 7]
    ]
    # EXPLICIT COLOR: Make the bounding box lines pure black
    box_colors = [[0.0, 0.0, 0.0] for _ in range(len(box_lines))] 
    
    bbox = o3d.geometry.LineSet()
    bbox.points = o3d.utility.Vector3dVector(box_points)
    bbox.lines = o3d.utility.Vector2iVector(box_lines)
    bbox.colors = o3d.utility.Vector3dVector(box_colors)

    # 3. Create small axes (size 50) in the corner
    frame = o3d.geometry.TriangleMesh.create_coordinate_frame(
        size=50.0, 
        origin=[-180.0, -180.0, -180.0]
    )

    # 4. Launch the modern Open3D GUI application
    app = gui.Application.instance
    app.initialize()

    window_title = f"3D Rama Plot: {category}"
    vis = o3d.visualization.O3DVisualizer(window_title, 1024, 768)
    
    # GUI CONTROLS: Turning this to True reveals the side panel with Reset buttons and size dials!
    vis.show_settings = True  

    # Set up rendering materials
    mat_pts = rendering.MaterialRecord()
    mat_pts.shader = "defaultUnlit"
    mat_pts.point_size = 4.0  # Slightly larger default point size
    
    mat_lines = rendering.MaterialRecord()
    mat_lines.shader = "unlitLine"
    mat_lines.line_width = 2.0

    # Add all geometry to the visualizer
    vis.add_geometry("Points", pcd, mat_pts)
    vis.add_geometry("BoundingBox", bbox, mat_lines)
    vis.add_geometry("Axes", frame, mat_pts)

    # Add the 3D Text Labels slightly past the tip of the 50-unit arrows
    vis.add_3d_label([-120, -180, -180], f"X: {x_angle}")
    vis.add_3d_label([-180, -120, -180], f"Y: {y_angle}")
    vis.add_3d_label([-180, -180, -120], f"Z: {z_angle}")

    # Render!
    vis.reset_camera_to_default()
    app.add_window(vis)
    app.run()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract angular data from JSONL files and visualize it in 3D space.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Examples:\n  python3 plot_pydangle.py my_data.jsonl\n  python3 plot_pydangle.py my_data.jsonl.gz -c PrePro\n  python3 plot_pydangle.py my_data.jsonl -x phi -y psi -z omega"
    )
    
    parser.add_argument("filepath", help="Path to the input JSONL or JSONL.gz file")
    parser.add_argument("-c", "--category", default="General", help="Ramachandran category to filter by (default: General)")
    parser.add_argument("-x", "--x-angle", default="phi", dest="x_angle", help="Data field for the X-axis / Red arrow (default: phi)")
    parser.add_argument("-y", "--y-angle", default="psi", dest="y_angle", help="Data field for the Y-axis / Green arrow (default: psi)")
    parser.add_argument("-z", "--z-angle", default="chi1", dest="z_angle", help="Data field for the Z-axis / Blue arrow (default: chi1)")

    args = parser.parse_args()
    visualize_ramachandran(args.filepath, args.category, args.x_angle, args.y_angle, args.z_angle)
