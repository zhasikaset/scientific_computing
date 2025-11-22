import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d


def draw_all(arrow=False, lvl=20):

    fig, ax = plt.subplots(2, 2, figsize=(8, 6))

    cnt = 0
    for row in range(2):
        for col in range(2):
            axs = ax[row, col]
            cf = axs.contourf(X, Y, results[cnt % 4], levels=lvl)
            fig.colorbar(cf, ax=axs)
            axs.set_title(names[cnt % 4])
            if arrow:
                axs.streamplot(X, Y, U, V, color="black")
                # axs.quiver(X[::1, ::1], Y[::1, ::1], U[::1, ::1], V[::1, ::1]) 

            axs.set_xlabel("x")
            axs.set_ylabel("y")
            cnt += 1

    plt.tight_layout()
    plt.show()
    
def draw_one(Z, name, arrow=False, lvl=20):
    fig, ax = plt.subplots()
    cf = ax.contourf(X, Y, Z, levels=lvl)
    fig.colorbar(cf, ax=ax)
    print(np.min(Z))
    if arrow:
        ax.streamplot(X, Y, U, V, color="black")
        # ax.quiver(X[::1, ::1], Y[::1, ::1], U[::1, ::1], V[::1, ::1]) 

    plt.title(name)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()


x = np.loadtxt(f"x.txt")
y = np.loadtxt(f"y.txt")

xi = np.linspace(x.min(), x.max(), x.size)
yi = np.linspace(y.min(), y.max(), y.size)

X, Y = np.meshgrid(xi, yi)

U = np.loadtxt(f"U.txt", delimiter="\t")
V = np.loadtxt(f"V.txt", delimiter="\t")
P = np.loadtxt(f"P.txt", delimiter="\t")

UV = np.sqrt(U*U + V*V)
results = [U, V, P, UV]
names = ["U", "V", "P", "U + V"]
print(UV.max(), np.abs(U.max()), np.abs(V.max()))

draw_all(arrow=True)

# draw_one(UV, "U + V", arrow=True)

print(X.shape, Y.shape, U.shape, V.shape, P.shape)

big_loaded = np.loadtxt("big.txt", delimiter="\t")


# Load each file (tab-separated)
U0 = np.loadtxt("source/0U.txt", delimiter="\t")
U1 = np.loadtxt("source/1U.txt", delimiter="\t")
U2 = np.loadtxt("source/2U.txt", delimiter="\t")
U3 = np.loadtxt("source/3U.txt", delimiter="\t")

# Combine them into a 2×2 block
big = np.block([
    [U0, U1],
    [U2, U3]
])

print("max difference between sequential and parallel solution is : ", np.max(np.abs(big-U)))