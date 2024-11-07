#!/usr/bin/env python
# coding: utf-8

# Noise analysis
# Predicted concentration vs estimated concentraion
import matplotlib.pyplot as plt
import numpy as np

from common import compute_T1, compute_c
from mixed import compute_ir_signal, compute_se_signal, LOOKUP_TABLE, extract_mixed_t1

np.random.seed(0)

# test different concentrations for a given noise level
num_samples = 200
c_values = np.linspace(0.01, 0.5, 200)
c = np.vstack([c_values for _ in range(num_samples)])
T1 = compute_T1(c)
SNR = 25
IR = compute_ir_signal(T1, SNR=SNR)
SE = compute_se_signal(T1, SNR=SNR)

T1_volume = extract_mixed_t1(IR=IR, SE=SE, lookup_table=LOOKUP_TABLE)
c_est = compute_c(T1_volume)

# scatter plot of the estimated concentration vs the true concentration
# and add the median, the mean, the identity line and the 95% confidence interval as a shaded area
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
ax[0].scatter(c, c_est, s=2)
c_max = np.max(c)
c_min = np.min(c)
ax[0].plot([c_min, c_max], [c_min, c_max], color="black", label="identity")
ax[0].plot(c_values, np.nanmedian(c_est, axis=0), color="red", label="mean")
ax[0].plot(c_values, np.nanmean(c_est, axis=0), color="blue", label="median")
ax[0].set_xlabel("True concentration (mmol/l)")
ax[0].set_ylabel("Estimated concentration")
ax[0].set_title(f"SNR = {SNR}")
ax[0].legend()

ax[1].scatter(c, T1_volume, s=2)
ax[1].plot(c_values, np.nanmedian(T1_volume, axis=0), color="red", label="mean")
ax[1].plot(c_values, np.nanmean(T1_volume, axis=0), color="blue", label="median")
ax[1].set_xlabel("True concentration (mmol/l)")
ax[1].set_ylabel("T1 estimate")
ax[1].set_title(f"SNR = {SNR}")
ax[1].legend()

plt.show()
