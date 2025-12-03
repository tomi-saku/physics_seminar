import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

data = np.loadtxt("orbit_1000_200_collision.txt")

num_frames = data.shape[0]
num_cols   = data.shape[1]

# 列数から粒子数を計算
num_particles = num_cols // 3

# 粒子ごとに (num_particles, num_frames, 3) にまとめる
points = np.zeros((num_particles, num_frames, 3))

for i in range(num_particles):
    points[i,:,0] = data[:, i*3 + 0]   # x
    points[i,:,1] = data[:, i*3 + 1]   # y
    points[i,:,2] = data[:, i*3 + 2]   # z

# --- 描画設定 ---
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(projection='3d')

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

ax.set_xlim([-5,5])
ax.set_ylim([-5,5])
ax.set_zlim([-5,5])

# unique size (最後の粒子を大きくしたい場合)
sizes = np.full(num_particles, 20)

colors = np.full(num_particles, 'b')
colors[-200:] = 'r'

# 初期位置をプロット
sc = ax.scatter(points[:,0,0], points[:,0,1], points[:,0,2],
                s=sizes, marker=".",c=colors)

def update(frame):
    sc._offsets3d = (points[:,frame,0],
                     points[:,frame,1],
                     points[:,frame,2])
    ax.set_title(f"Frame {frame}")
    return sc,

ani = animation.FuncAnimation(fig, update, frames=range(0, num_frames, 10), interval=10)

ani.save('orbit_1000_200_collision_animation.mp4', writer='ffmpeg', fps=20)

plt.show()
