"""複数のグラフを重ねて描画するプログラム"""
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

data = np.loadtxt("energy_1000_leapfrog_cold.txt")

time = data[:,0]
kinetic_energy = data[:,1]
potential_energy = data[:,2]
whole_energy = data[:,1] + data[:,2]
virial = 2*data[:,1] + data[:,2]

c1,c2,c3 = "blue","green","red"      # 各プロットの色
l1,l2,l3= "K","V","2K+V"  # 各ラベル

ax.set_xlabel('time')  # x軸ラベル
ax.set_ylabel('Energy')  # y軸ラベル
ax.set_title("Kinetic and Potential energy(cold)") # グラフタイトル
# ax.set_aspect('equal') # スケールを揃える
ax.grid()            # 罫線
ax.set_xlim([0, (time[1]-time[0])*len(time)]) # x方向の描画範囲を指定
# ax.set_ylim([0, 1])    # y方向の描画範囲を指定
ax.plot(time, kinetic_energy, color=c1, label=l1)
ax.plot(time, potential_energy, color=c2, label=l2)
ax.plot(time, virial, color=c3, label=l3)


# 目標値
plt.axvline(x=0.035,color='r', ls='--')

ax.legend(loc=0)    # 凡例
fig.tight_layout()  # レイアウトの設定
plt.savefig('energy_1000_subplot_cold_virial.png') # 画像の保存
plt.show()