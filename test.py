import Grid_Operations as GO
import numpy as np
import matplotlib.pyplot as plt

grid, x= GO.create_grid(10, 36, 'linear')
segments = GO.segment_3dgrid(grid)

ax = plt.figure().add_subplot(projection='3d')

for i, (X, Y, Z) in enumerate(segments):
    print(X.shape, Y.shape, Z.shape, i)
    if i == 13:
        continue
    ax.scatter(X, Y, Z, alpha=0.1)
plt.show()