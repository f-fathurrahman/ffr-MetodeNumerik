# Import libraries
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

ax.quiver(0.0, 0.0, 1.0, 2.0, color="blue")
ax.quiver(-1.0, -1.0, -1.0, 2.0, color="red")

ax.set_title("Simple quiver plot")
ax.set_xlim(-10.0, 10.0)
ax.set_ylim(-10.0, 10.0)

ax.set_aspect("equal", "box")
ax.grid(True)

# Show plot
plt.show()
