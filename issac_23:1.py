#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 14:39:44 2025

@author: kkhissac
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import ifft, fftfreq
from scipy.integrate import cumtrapz, solve_ivp

# Parameters for road input
kappa = 5e-7  # Road roughness coefficient (m^3/cycle)
V = 25        # Vehicle speed (m/s)
Fs = 1000    # Sampling frequency (Hz)
T = 100       # Simulation time (seconds)
N = int(Fs * T)  # Number of samples
dt = 1 / Fs   # Time step

# Frequency vector
f = fftfreq(N, d=dt)[:N//2+1]
omega = 2 * np.pi * f

# Velocity power spectrum (S_dot_xr)
S_dot_xr = 2 * np.pi * kappa * V  # Constant value for all frequencies

# Generate Gaussian white noise in frequency domain
random_phase = np.exp(1j * 2 * np.pi * np.random.rand(len(f)))
amplitude = np.sqrt(S_dot_xr * Fs / 2)  # Scale by power spectrum
white_noise_freq = amplitude * random_phase

# Mirror spectrum for negative frequencies
white_noise_full = np.concatenate([white_noise_freq, np.conj(white_noise_freq[-2:0:-1])])

# Inverse FFT to get road velocity in time domain
road_velocity = ifft(white_noise_full).real

# Time vector
t = np.linspace(0, T, N, endpoint=False)

# Integrate to get road displacement
road_displacement = cumtrapz(road_velocity, t, initial=0)

# Define system parameters for quarter car model
m1 = 240  # Sprung mass (kg)
m2 = 36   # Unsprung mass (kg)
k = 16e3  # Suspension stiffness (N/m)
kt = 160e3  # Tire stiffness (N/m)
c = 100    # Damping coefficient (N·s/m)
def quarter_car_model(t, y):
    # y = [x1, v1, x2, v2] (displacement and velocity of sprung and unsprung mass)
    x1, v1, x2, v2 = y

    # Road displacement (interpolated from road_displacement array)
    xr = np.interp(t, t_eval, road_displacement)

    # Equations of motion
    dx1dt = v1
    dv1dt = (-c * (v1 - v2) - k * (x1 - x2)) / m1
    dx2dt = v2
    dv2dt = (c * (v1 - v2) + k * (x1 - x2) - kt * (x2 - xr)) / m2

    return [dx1dt, dv1dt, dx2dt, dv2dt]

# Initial conditions: [x1, v1, x2, v2]
y0 = [0, 0, 0, 0]

# Solve the quarter car model using solve_ivp
t_eval = np.linspace(0, T, N)
solution = solve_ivp(quarter_car_model, [0, T], y0, t_eval=t_eval, method='RK45')

# Extract results
x1 = solution.y[0]  # Displacement of sprung mass
v1 = solution.y[1]  # Velocity of sprung mass
a1 = np.gradient(v1, t_eval)  # Acceleration of sprung mass

# Plot results
plt.figure(figsize=(12, 8))

# Road displacement
plt.subplot(3, 1, 1)
plt.plot(t, road_displacement, label='Road Displacement')
plt.title('Road Displacement Input')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.grid(True)
plt.legend()

# Sprung mass displacement
plt.subplot(3, 1, 2)
plt.plot(t_eval, x1, label='Sprung Mass Displacement')
plt.title('Sprung Mass Displacement')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.grid(True)
plt.legend()

# Sprung mass acceleration
plt.subplot(3, 1, 3)
plt.plot(t_eval, a1, label='Sprung Mass Acceleration')
plt.title('Sprung Mass Acceleration')
plt.xlabel('Time (s)')
plt.ylabel('Acceleration (m/s²)')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()


from scipy.fft import fft, fftfreq

# Compute FFT of the acceleration profile
acceleration_fft = fft(a1)  # `a1` is the acceleration of the sprung mass
freqs = fftfreq(len(t_eval), d=(t_eval[1] - t_eval[0]))  # Frequencies in Hz

# Take only the positive frequencies
positive_freqs = freqs[freqs >= 0]
acceleration_magnitude = np.abs(acceleration_fft[freqs >= 0])

# Normalize the acceleration magnitude
acceleration_magnitude_normalized = acceleration_magnitude / len(t_eval)
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.plot(positive_freqs, acceleration_magnitude_normalized, label='Acceleration vs Frequency')
plt.xlim(0,30)
plt.title('Sprung Mass Acceleration vs Frequency')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Acceleration Magnitude (m/s²)')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.legend()
plt.show()
