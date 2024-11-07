import numpy as np
import matplotlib.pyplot as plt
from look_locker import curve_fit_wrapper, OptimizeWarning
from common import compute_T1, compute_c

def f(t, x1, x2, x3):
    return np.abs(x1 * (1.0 - (1.1 + x2**2) * np.exp(-(x3**2) * t)))

values = 5
T1 = np.linspace(0.3, 5.0, values**2)
x1 = 1.0
x2 = np.sqrt(0.5)  # todo: what to choose here? not unique with x3
x3 = np.sqrt((0.1 + x2**2) / T1)

num_samples = 14
sequence_length = 2.6  # seconds

fig, axes = plt.subplots(values, values, figsize=(values * 3, values * 3))

for x3_val, ax in zip(x3, axes.flatten()):
    t = np.linspace(0.115, sequence_length, num_samples)
    y = f(t, x1, x2, x3_val)

    # add some Rician noise to the data
    SNR = 25
    max_val = np.max(y)
    a, b = (
        np.random.normal(0, max_val / SNR, len(y)),
        np.random.normal(0, max_val / SNR, len(y)),
    )
    y = y + np.sqrt(a**2 + b**2)
    ax.plot(t, y, "o")

    popt = curve_fit_wrapper(f, t, y, p0=(1.0, 2.0, 1.0))
    t = np.linspace(0.115, sequence_length, 1000)
    ax.plot(t, f(t, *popt))

    ax.set_title(
        f"T1: {(0.1 + popt[1]**2)/(popt[2]**2):.2f} s \n"
        f"Actual T1: {(0.1 + x2**2)/(x3_val**2):.2f} s"
    )

fig.tight_layout()
plt.show()

# do the analysis for a range of T1 values and for each T1 multiple noise realizations
def compute_T1_from_x23(x2, x3):
    return (0.1 + x2**2) / (x3**2)


c_values = np.linspace(0.001, 0.3, 100)
T1 = compute_T1(c_values) * 0.001  # convert to seconds

x1 = 1.0
x2 = np.sqrt(0.5)  # todo: what to choose here? not uniquely defined by T1
x3 = np.sqrt((0.1 + x2**2) / T1)

n = 50
popt = np.zeros((len(T1), n, 3))
for i, x3 in enumerate(x3):
    for j in range(n):
        t = np.linspace(0, sequence_length, num_samples)
        y = f(t, x1, x2, x3)
        SNR = 25
        max_val = np.max(y)
        re, im = (
            np.random.normal(0, max_val / SNR, len(y)),
            np.random.normal(0, max_val / SNR, len(y)),
        )
        y = y + np.sqrt(re**2 + im**2)
        try:
            popt[i, j] = curve_fit_wrapper(f, t, y, (1.0, np.sqrt(0.9), 1.0))
        except (OptimizeWarning, RuntimeError):
            popt[i, j] = np.nan * np.zeros(3)

# calculate the mean and standard deviation of the estimated T1 values
T1_est = compute_T1_from_x23(x2=popt[:, :, 1], x3=popt[:, :, 2]) * 1000  # convert to ms
c_est = compute_c(T1_est)

c_est_mean = np.mean(c_est, axis=1)
c_est_std = np.std(c_est, axis=1)

fig, ax = plt.subplots(1, 3, figsize=(15, 5))
ax[0].plot(c_values, c_values, "k--", lw=2)
ax[0].plot(c_values, c_est_mean, "r-", label="mean")
for i in range(len(c_values)):
    ax[0].scatter([c_values[i]] * n, c_est[i], s=2, c="b", alpha=0.05)
ax[0].set_xlabel("Actual c (mmol/l)")
ax[0].set_ylabel("Estimated c (mmol/l)")
ax[0].set_xlim([0, 0.3])
ax[0].set_ylim([0, 0.5])
ax[0].legend()

# plot estimated T1 values over true concentration
T1_est_mean = np.mean(T1_est, axis=1)
T1_est_std = np.std(T1_est, axis=1)

ax[1].plot(c_values, T1 * 1000, "k--", lw=2)
ax[1].plot(c_values, T1_est_mean, "r-", label="mean")
for i in range(len(T1)):
    ax[1].scatter([c_values[i]] * n, T1_est[i], s=2, c="b", alpha=0.05)
ax[1].set_xlabel("Actual c (mmol/l)")
ax[1].set_ylabel("Estimated T1 (ms)")
ax[1].set_xlim([0, 0.3])
ax[1].set_ylim([0, 5000])
ax[1].legend()

ax[2].plot(T1 * 1000, T1 * 1000, "k--", lw=2)
ax[2].plot(T1 * 1000, T1_est_mean, "r-", label="mean")
for i in range(len(T1)):
    ax[2].scatter([T1[i] * 1000] * n, T1_est[i], s=2, c="b", alpha=0.05)
ax[2].set_xlabel("Actual T1 (ms)")
ax[2].set_ylabel("Estimated T1 (ms)")
ax[2].set_xlim([0, 5000])
ax[2].set_ylim([0, 5000])
ax[2].legend()

fig.tight_layout()
plt.show()

print("T1 estimated from noisy signal (SNR=25) vs actual T1:\n")
print(f"mean stddev for all values: {np.mean(T1_est_std):.2f}ms")
print(f"mean stddev for T1 < 2000ms: {np.mean(T1_est_std[T1<2.0]):.2f}ms")
print(f"mean stddev for T1 > 2000ms: {np.mean(T1_est_std[T1>2.0]):.2f}ms")
print(f"mean stddev/value for all values: {np.mean(T1_est_std/(T1*1000))*100:.2f}%")
print(
    f"mean stddev/value for T1 < 2000ms: {np.mean(T1_est_std[T1<2]/(T1[T1<2]*1000))*100:.2f}%"
)
print(
    f"mean stddev/value for T1 > 2000ms: {np.mean(T1_est_std[T1>2]/(T1[T1>2]*1000))*100:.2f}%"
)
