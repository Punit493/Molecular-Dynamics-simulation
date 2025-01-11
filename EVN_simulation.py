import numpy as np
import matplotlib.pyplot as plt

# Define constants
sig = 1
ep = 1
m = 1
lx = 20
ly = lx
lz = lx
num_particles = 512
ni = 1500
dt = 0.005
va = np.sqrt(12)
rc = 2.5 * sig
t2 = (dt**2) / (2 * m)
t1 = dt / (2 * m)

# Function to initialize position
def posinit():
    pl = 2
    dp = 2
    p2 = 16
    kk = 0
    x = np.zeros(num_particles)
    y = np.zeros(num_particles)
    z = np.zeros(num_particles)
    for p in range(pl, p2 + 1, dp):
        for q in range(pl, p2 + 1, dp):
            for r in range(pl, p2 + 1, dp):
                x[kk] = p
                y[kk] = q
                z[kk] = r
                kk += 1
    return x, y, z

# Function to initialize velocity
def velinit():
    vx = va * (np.random.rand(num_particles) - 0.5)
    vy = va * (np.random.rand(num_particles) - 0.5)
    vz = va * (np.random.rand(num_particles) - 0.5)

    # Calculate average velocity along x, y, z directions
    avx = np.mean(vx)
    avy = np.mean(vy)
    avz = np.mean(vz)

    # Adjust velocity components so that their mean is zero
    vx -= avx
    vy -= avy
    vz -= avz

    return vx, vy, vz

# Function to calculate forces using Lennard-Jones potential
def forcecalc(x, y, z):
    fx = np.zeros(num_particles)
    fy = np.zeros(num_particles)
    fz = np.zeros(num_particles)
    PEN = 0
    for i in range(num_particles - 1):
        x1, y1, z1 = x[i], y[i], z[i]
        for j in range(i + 1, num_particles):
            x2, y2, z2 = x[j], y[j], z[j]
            dx, dy, dz = x1 - x2, y1 - y2, z1 - z2
            if abs(dx) > lx / 2:
                dx = (lx - abs(dx)) * (-dx) / abs(dx)
            if abs(dy) > ly / 2:
                dy = (ly - abs(dy)) * (-dy) / abs(dy)
            if abs(dz) > lz / 2:
                dz = (lz - abs(dz)) * (-dz) / abs(dz)
            r = np.sqrt(dx**2 + dy**2 + dz**2)
            if r < rc:
                f = (24 / r**7) * ((2 / r**6) - 1)
                fx[i] += f * (dx / r)
                fx[j] -= f * (dx / r)
                fy[i] += f * (dy / r)
                fy[j] -= f * (dy / r)
                fz[i] += f * (dz / r)
                fz[j] -= f * (dz / r)
                PEN += 4 * ((1 / r**12) - (1 / r**6))
    return fx, fy, fz, PEN

# Function to update positions
def updatepos(x, y, z, vx, vy, vz, fx, fy, fz):
    x = (x + vx * dt + fx * t2) % lx
    y = (y + vy * dt + fy * t2) % lx
    z = (z + vz * dt + fz * t2) % lx
    return x, y, z

# Function to update velocities
def updatevel(vx, vy, vz, fx, fy, fz, fox, foy, foz):
    vx += (fox + fx) * t1
    vy += (foy + fy) * t1
    vz += (foz + fz) * t1
    return vx, vy, vz

# Main Program
x, y, z = posinit()
xl, yl, zl = np.copy(x), np.copy(y), np.copy(z)
vx, vy, vz = velinit()
fx, fy, fz, PEN = forcecalc(x, y, z)

PE = []
KE = []
for k in range(ni):
    x, y, z = updatepos(x, y, z, vx, vy, vz, fx, fy, fz)
    fox, foy, foz = np.copy(fx), np.copy(fy), np.copy(fz)  # Save the current forces
    fx, fy, fz, PEN = forcecalc(x, y, z)
    vx, vy, vz = updatevel(vx, vy, vz, fx, fy, fz, fox, foy, foz)
    KEN = 0.5 * np.sum(vx**2 + vy**2 + vz**2) / num_particles
    PE.append(PEN / num_particles)
    KE.append(KEN)
    print(f'Iterations: {k+1}')

# Plotting results
fig = plt.figure(figsize=(12, 10))

# Initial positions
ax1 = fig.add_subplot(2, 2, 1, projection='3d')
ax1.scatter(xl, yl, zl, c='blue', marker='o')
ax1.set_title('Initial positions of particles')

# Final positions
ax2 = fig.add_subplot(2, 2, 2, projection='3d')
ax2.scatter(x, y, z, c='red', marker='o')
ax2.set_title('Final positions of particles')

# Energy plot
TE = np.array(KE) + np.array(PE)
ax3 = fig.add_subplot(2, 2, 3)
ax3.plot(range(ni), KE, label='Kinetic Energy')
ax3.plot(range(ni), PE, label='Potential Energy')
ax3.plot(range(ni), TE, label='Total Energy')
ax3.legend()
ax3.set_xlabel('Iterations')
ax3.set_ylabel('Energy')

# Speed distribution
velocity = np.sqrt(vx**2 + vy**2 + vz**2)
ax4 = fig.add_subplot(2, 2, 4)
ax4.hist(velocity, bins=20, edgecolor='black')
ax4.set_xlabel('Speed of the particles')

plt.tight_layout()
plt.show()

print("Microcanonical simulation complete.")
