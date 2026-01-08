import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ======================
# データ読み込み
# ======================
position_data = np.loadtxt("orbit_collision_metric.txt")
velocity_data = np.loadtxt("velocity_collision_metric.txt")

num_frames = position_data.shape[0]
num_cols   = position_data.shape[1]
num_particles = num_cols // 3

# (num_particles, num_frames, 3)
position = np.zeros((num_particles, num_frames, 3))
velocity = np.zeros((num_particles, num_frames, 3))

for i in range(num_particles - 1):
    position[i,:,0] = position_data[:, i*3 + 0]
    position[i,:,1] = position_data[:, i*3 + 1]
    position[i,:,2] = position_data[:, i*3 + 2]

    velocity[i,:,0] = velocity_data[:, i*3 + 0]
    velocity[i,:,1] = velocity_data[:, i*3 + 1]
    velocity[i,:,2] = velocity_data[:, i*3 + 2]

# ======================
# 距離・速度の大きさを事前計算（高速化）
# ======================
r = np.linalg.norm(position, axis=2)   # (num_particles, num_frames)

# 追加：速度と位置の内積 / |r|
vr = np.zeros_like(r)
mask = r > 0
vr[mask] = np.sum(position[mask] * velocity[mask], axis=1) / r[mask]

# 軸範囲を固定
r_max = r.max()
vr_max = vr.max()
vr_min = vr.min()

# ======================
# プロット初期化
# ======================
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlabel(r"$Distance from center |\mathbf{r}|$")
ax.set_ylabel(r"$(\mathbf{v}\cdot\mathbf{r}) / |\mathbf{r}|$")
ax.set_xlim(0, r_max * 1.05)
ax.set_ylim(vr_min * 1.05, vr_max * 1.05)
ax.set_title("Phase Space (r - v)")

# 初期フレーム
scat = ax.scatter(
    r[:, 0],
    vr[:, 0],
    s=5,
    alpha=0.6
)

# ======================
# アニメーション更新関数
# ======================
def update(frame):
    # (N, 2) の形で更新
    offsets = np.column_stack((r[:, frame], vr[:, frame]))
    scat.set_offsets(offsets)
    ax.set_title(f"Phase Space (frame = {frame})")
    return scat,

# ======================
# アニメーション作成
# ======================
ani = animation.FuncAnimation(
    fig,
    update,
    frames=range(0, num_frames, 10),
    interval=30,
    blit=True
)

# ======================
# 保存（必要なら）
# ======================
ani.save("phase_space_r_vr.mp4", writer="ffmpeg",dpi=200, fps=30)

plt.show()
