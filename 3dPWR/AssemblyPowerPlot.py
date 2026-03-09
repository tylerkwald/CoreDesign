import re
import matplotlib.pyplot as plt
import numpy as np
import os

def extract_and_plot_power(file_path):
    with open(file_path, 'r') as f:
        content = f.read()

    # Locate the "Average Planar Power Distribution" table
    # This regex looks for the header and then captures the subsequent rows of data
    pattern = re.compile(
        r"Assembly Power Distribution.*?box power\s+1\s+2\s+3\s+4\s+5\s+6\s+7\s+8(.*?)(?=Maximum Pos)", 
        re.DOTALL
    )
    
    match = pattern.search(content)
    if not match:
        print("Could not find the Average Planar Power Distribution table in the file.")
        return

    table_data = match.group(1).strip()
    
    # Initialize an 8x8 grid with NaNs (for non-fuel regions)
    grid = np.full((8, 8), np.nan)
    
    # Parse the table lines
    lines = table_data.split('\n')
    for line in lines:
        parts = line.split()
        if not parts: continue
        
        try:
            row_idx = int(parts[0]) - 1  # Row label (1-8)
            values = [float(v) for v in parts[1:]]
            for col_idx, val in enumerate(values):
                if col_idx < 8:
                    grid[row_idx, col_idx] = val
        except ValueError:
            continue

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(grid, cmap='autumn_r', interpolation='nearest')
    
    # Add labels and colorbar
    plt.colorbar(im, label='Relative Assembly Power')
    ax.set_xticks(np.arange(8))
    ax.set_yticks(np.arange(8))
    ax.set_xticklabels(np.arange(1, 9))
    ax.set_yticklabels(np.arange(1, 9))
    ax.set_xlabel('Assembly Column (I)')
    ax.set_ylabel('Assembly Row (J)')
    ax.set_title('Average Planar Power Distribution (IAEA 3D)')

    # Annotate each cell with its power value
    for i in range(8):
        for j in range(8):
            if not np.isnan(grid[i, j]):
                ax.text(j, i, f'{grid[i, j]:.4f}', ha="center", va="center", color="black", fontsize=8)

    plt.tight_layout()
    output_folder = 'images'
    file_name = 'planar_power_distribution.png'

    # Create the full path
    save_path = os.path.join(output_folder, file_name)

    # Save the figure
    # bbox_inches='tight' ensures the labels aren't cut off
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved successfully to: {save_path}")

# Replace with the path to your file
extract_and_plot_power('IAEA3D.out')