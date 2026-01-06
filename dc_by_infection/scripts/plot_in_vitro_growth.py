import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Data
data = {
    'JW18DOX': [5.313432836, 4.964788732, 5.013303769, 5.291005291, 
                4.542553191, 5.148514851],
    'JW18wMel': [1.660194175, 1.409356725, 2.633027523, 1.926267281, 
                 1.926267281, 2.050847458],
    'S2DOX': [9.77593361, 5.583157895, 3.606493506, 5.98489426, 
              6.804177546, 4.699335548],
    'S2wMel': [2.584183673, 2.517113784, 1.85942029, 2.433839479, 
               3.139255702, 2.508213552]
}

# Color palette
colors = {
    'JW18DOX': '#8fcb84',
    'JW18wMel': '#09aa4b',
    'S2DOX': '#fab280',
    'S2wMel': '#d25727'
}

# Prepare data for plotting
conditions = list(data.keys())
values = [data[cond] for cond in conditions]
color_list = [colors[cond] for cond in conditions]

# Create figure
fig, ax = plt.subplots(figsize=(10, 6))

# Create box plot
bp = ax.boxplot(values, labels=conditions, patch_artist=True,
                widths=0.6, showfliers=False,
                boxprops=dict(linewidth=2),
                medianprops=dict(color='black', linewidth=2),
                whiskerprops=dict(linewidth=2),
                capprops=dict(linewidth=2))

# Color the box outlines and make backgrounds transparent
for patch, color in zip(bp['boxes'], color_list):
    patch.set_facecolor('none')  # No fill
    patch.set_edgecolor(color)   # Colored outline
    patch.set_linewidth(2)

# Color the whiskers and caps
for i, color in enumerate(color_list):
    bp['whiskers'][i*2].set_color(color)
    bp['whiskers'][i*2+1].set_color(color)
    bp['caps'][i*2].set_color(color)
    bp['caps'][i*2+1].set_color(color)
    bp['whiskers'][i*2].set_linewidth(2)
    bp['whiskers'][i*2+1].set_linewidth(2)
    bp['caps'][i*2].set_linewidth(2)
    bp['caps'][i*2+1].set_linewidth(2)

# Add individual points with jitter
np.random.seed(42)  # for reproducible jitter
for i, (cond, vals) in enumerate(zip(conditions, values)):
    # Add jitter to x-coordinates
    x = np.random.normal(i+1, 0.04, size=len(vals))
    ax.scatter(x, vals, alpha=0.8, s=60, color=colors[cond], 
               edgecolors='black', linewidth=0.5, zorder=3)

# Customize plot
ax.set_ylabel('Value', fontsize=12, fontweight='bold')
ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
ax.set_title('Box Plot with Individual Data Points', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3, linestyle='--')
ax.set_axisbelow(True)

# Adjust layout
plt.tight_layout()

# Save figure
plt.savefig('boxplot_with_points.png', dpi=300, bbox_inches='tight')
plt.savefig('boxplot_with_points.pdf', bbox_inches='tight')

# Show plot
plt.show()

print("Plot saved as 'boxplot_with_points.png' and 'boxplot_with_points.pdf'")
