import argparse
import json
import gzip
import sys

def generate_kinemage(input_file, output_file, category, x_angle, y_angle, z_angle):
    points = []
    print(f"Parsing '{category}' category from {input_file}...")
    
    open_func = gzip.open if input_file.endswith('.gz') else open
    mode = 'rt' if input_file.endswith('.gz') else 'r'

    try:
        with open_func(input_file, mode) as f:
            for line in f:
                if not line.strip(): continue
                data = json.loads(line)
                
                if (data.get('rama_category') == category and
                    data.get(x_angle) is not None and
                    data.get(y_angle) is not None and
                    data.get(z_angle) is not None):
                    
                    chain = data.get('chain', '_')
                    resnum = data.get('resnum', '')
                    resname = data.get('resname', 'UNK')
                    label = f"{chain} {resnum} {resname}".strip()
                    points.append((label, data[x_angle], data[y_angle], data[z_angle]))
                    
    except FileNotFoundError:
        print(f"Error: Could not find file '{input_file}'.")
        sys.exit(1)

    if not points:
        print(f"No valid points found for category '{category}'.")
        sys.exit(0)

    print(f"Writing {len(points)} points to {output_file}...")

    with open(output_file, 'w') as out:
        out.write("@text\n")
        out.write(f"Ramachandran Point Cloud\n")
        out.write(f"Category: {category}\n")
        out.write(f"Axes: X={x_angle}, Y={y_angle}, Z={z_angle}\n")
        out.write("@kinemage 1\n")
        out.write(f"@title {{{category} Cloud}}\n")
        
        # Standard XY View
        out.write("@viewid {XY Plane (Standard Rama)}\n")
        out.write("@span 550.0\n")
        out.write("@center 0.0 0.0 0.0\n")
        out.write("@matrix 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0\n")
        
        # Bounding Box
        out.write("@group {Bounding Box} dominant\n")
        out.write("@vectorlist {Box} color= gray\n")
        out.write("{-180 -180 -180} -180.0 -180.0 -180.0\n{ 180 -180 -180}  180.0 -180.0 -180.0\n{ 180  180 -180}  180.0  180.0 -180.0\n{-180  180 -180} -180.0  180.0 -180.0\n{-180 -180 -180} -180.0 -180.0 -180.0\n")
        out.write("@vectorlist {Box Top} color= gray\n")
        out.write("{-180 -180  180} -180.0 -180.0  180.0\n{ 180 -180  180}  180.0 -180.0  180.0\n{ 180  180  180}  180.0  180.0  180.0\n{-180  180  180} -180.0  180.0  180.0\n{-180 -180  180} -180.0 -180.0  180.0\n")
        out.write("@vectorlist {Box Pillars} color= gray\n")
        out.write("P {-180 -180 -180} -180.0 -180.0 -180.0\n  {-180 -180  180} -180.0 -180.0  180.0\n")
        out.write("P { 180 -180 -180}  180.0 -180.0 -180.0\n  { 180 -180  180}  180.0 -180.0  180.0\n")
        out.write("P { 180  180 -180}  180.0  180.0 -180.0\n  { 180  180  180}  180.0  180.0  180.0\n")
        out.write("P {-180  180 -180} -180.0  180.0 -180.0\n  {-180  180  180} -180.0  180.0  180.0\n")

        # Full-length Coordinate Axes
        out.write("@group {Axes & Labels} dominant\n")
        out.write("@vectorlist {X-axis} color= red width= 2\n")
        out.write("P -180.0 -180.0 -180.0\n  180.0 -180.0 -180.0\n")
        out.write("@vectorlist {Y-axis} color= green width= 2\n")
        out.write("P -180.0 -180.0 -180.0\n -180.0  180.0 -180.0\n")
        out.write("@vectorlist {Z-axis} color= blue width= 2\n")
        out.write("P  180.0 -180.0 -180.0\n  180.0 -180.0  180.0\n")

        # Major Tick-marks (every 90 degrees, length 10)
        out.write("@vectorlist {Major Ticks} color= white\n")
        for tick in [-180, -90, 0, 90, 180]:
            out.write(f"P {tick}.0 -180.0 -180.0\n  {tick}.0 -190.0 -180.0\n") # X ticks down
        for tick in [-90, 0, 90, 180]: 
            out.write(f"P -180.0 {tick}.0 -180.0\n  -190.0 {tick}.0 -180.0\n") # Y ticks left
        for tick in [-90, 0, 90, 180]: 
            out.write(f"P 180.0 -180.0 {tick}.0\n  190.0 -180.0 {tick}.0\n") # Z ticks right

        # Minor Tick-marks (every 45 degrees, length 5)
        out.write("@vectorlist {Minor Ticks} color= white\n")
        for tick in [-135, -45, 45, 135]:
            out.write(f"P {tick}.0 -180.0 -180.0\n  {tick}.0 -185.0 -180.0\n") # X minor ticks
            out.write(f"P -180.0 {tick}.0 -180.0\n  -185.0 {tick}.0 -180.0\n") # Y minor ticks
            out.write(f"P 180.0 -180.0 {tick}.0\n  185.0 -180.0 {tick}.0\n") # Z minor ticks

        # Cleanly Spaced 3D Text Labels
        out.write("@labellist {Text Labels} color= white\n")
        
        # Axis Titles
        out.write(f"{{{x_angle}}} 0.0 -225.0 -180.0\n")   # Centered below X
        out.write(f"{{{y_angle}}} -225.0 0.0 -180.0\n")   # Centered left of Y
        out.write(f"{{{z_angle}}} 225.0 -180.0 0.0\n")    # Centered right of Z
        
        # Axis numbers (Omitting 180 and -180)
        for tick in [-90, 0, 90]:
            out.write(f"{{{tick}}} {tick}.0 -205.0 -180.0\n") # X-axis
            out.write(f"{{{tick}}} -205.0 {tick}.0 -180.0\n") # Y-axis
            out.write(f"{{{tick}}} 205.0 -180.0 {tick}.0\n")  # Z-axis

        # High-performance point rendering (Dots only)
        out.write(f"@group {{Data}} dominant\n")
        out.write(f"@subgroup {{{category}}} dominant\n")
        out.write(f"@dotlist {{{category} Points}} color= sky\n")
            
        for p in points:
            out.write(f"{{{p[0]}}} {p[1]:.3f} {p[2]:.3f} {p[3]:.3f}\n")

    print(f"Done! You can now open {output_file} in KiNG.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert JSONL angle data to Kinemage format.")
    
    parser.add_argument("filepath", help="Path to input JSONL or JSONL.gz file")
    parser.add_argument("-o", "--output", default="output.kin", help="Output kinemage filename")
    parser.add_argument("-c", "--category", default="General", help="Ramachandran category")
    parser.add_argument("-x", "--x-angle", default="phi", dest="x_angle", help="X-axis (default: phi)")
    parser.add_argument("-y", "--y-angle", default="psi", dest="y_angle", help="Y-axis (default: psi)")
    parser.add_argument("-z", "--z-angle", default="chi1", dest="z_angle", help="Z-axis (default: chi1)")

    args = parser.parse_args()
    generate_kinemage(args.filepath, args.output, args.category, args.x_angle, args.y_angle, args.z_angle)
