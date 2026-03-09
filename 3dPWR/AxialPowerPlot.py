import re
import os
import matplotlib.pyplot as plt

def extract_and_plot_axial_v2(file_path):
    # Ensure the images directory exists
    output_dir = "images"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(file_path, 'r') as f:
        lines = f.readlines()

    heights = []
    powers = []
    
    # Flag to start capturing once we hit the Axial Power Distribution header
    capture = False
    
    for line in lines:
        if "Axial Power Distribution" in line:
            capture = True
            continue
        
        if capture:
            # Match the specific pattern: Height Power Plane Mesh
            # Example: 10.00 0.316 2 20.00
            match = re.search(r"^\s*([\d\.]+)\s+([\d\.]+)\s+(\d+)\s+([\d\.]+)", line)
            if match:
                heights.append(float(match.group(1)))
                powers.append(float(match.group(2)))
            elif heights and not line.strip():
                # Stop capturing if we hit an empty line after finding data
                break

    if not powers:
        print("Error: Could not extract axial data. Check if the table is present in the file.")
        return

    # Plotting
    plt.figure(figsize=(7, 9))
    
    # Plotting Power (X) vs Height (Y)
    plt.plot(powers, heights, marker='s', linestyle='-', color='#d62728', linewidth=2, label='Axial Power')
    
    # Aesthetics
    plt.xlabel('Relative Planar Power', fontsize=12)
    plt.ylabel('Core Height [cm]', fontsize=12)
    plt.title('IAEA 3D Benchmark: Axial Power Distribution', fontsize=14)
    plt.grid(True, which='both', linestyle='--', alpha=0.6)
    plt.legend()

    # Save to images folder
    save_path = os.path.join(output_dir, "axial_power_distribution.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    print(f"Extraction successful. Plot saved to: {save_path}")
    plt.show()

# Run the script with your specific file
extract_and_plot_axial_v2('IAEA3D.out')