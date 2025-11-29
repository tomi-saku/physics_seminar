import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


data = np.loadtxt("plummer_initial_position.txt")

# Steps = data.shape[0]
# N = data.shape[1]

x = []
y = []
z = []

print(data.shape)

for i in range(data.shape[1] - 2):
  x.append(data[:,i])
  y.append(data[:,i+1])
  z.append(data[:,i+2])

fig =plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection='3d')

ax.set_xlabel("X position")
ax.set_ylabel("Y position")
ax.set_zlabel("Z position")

ax.set_xlim([-10,10])
ax.set_ylim([-10,10])
ax.set_zlim([-10,10])

plt.title("2-body Orbit using Euler Method")
plt.axis("equal")
plt.grid(True)

print(len(x))
print(len(x[1]))

ims = []
for j in range(len(x)):
  im = ax.plot(x[j],y[j],z[j],linewidth=0,marker=".")
  ims.append(im)

def update(frame):
  ax.cla()
  ax.set_xlabel("X position")
  ax.set_ylabel("Y position")
  ax.set_zlabel("Z position")
  ax.set_xlim([-5,5])
  ax.set_ylim([-5,5])
  ax.set_zlim([-5,5])
  ax.set_title(f"Frame {frame}")
  for i in range(len(x)):
        ax.plot(x[i][frame], y[i][frame], z[i][frame], linewidth=0, marker=".")

ani = animation.FuncAnimation(fig, update, frames=len(x[1]), interval=10)

# mp4ファイルとして保存
ani.save('orbit_animation1000.mp4', writer='ffmpeg', fps=20)

plt.show()