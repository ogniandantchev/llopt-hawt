import matplotlib.pyplot as plt
import re

sections = []
current = {'x': [], 'y': [], 'label': ''}

with open('blade.plt', 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            if current['x']:
                sections.append(current)
                current = {'x': [], 'y': [], 'label': ''}
        elif line.startswith('# r/R'):
            m = re.search(r'r/R\s*=\s*([\d.]+).*twist\s*=\s*([-\d.]+)', line)
            if m:
                current['label'] = f"r/R={m.group(1)}  twist={m.group(2)}°"
        elif not line.startswith('#'):
            parts = line.split()
            if len(parts) == 2:
                current['x'].append(float(parts[0]))
                current['y'].append(float(parts[1]))

if current['x']:
    sections.append(current)

fig, ax = plt.subplots(figsize=(10, 6))
for s in sections:
    ax.plot(s['x'], s['y'], label=s['label'])
ax.set_aspect('equal')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title('Blade Sections (twisted, physical scale)')
ax.legend(fontsize=8, loc='upper right')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
