import numpy as np 

def synthetic_forecast(predefined_series, attract_factor=0.15, noise_strength=0.05, cap_zero=False, random_seed=None):
    # Set random seed
    if random_seed is not None:
        np.random.seed(random_seed)

    num_steps = len(predefined_series)
    t = np.arange(num_steps)  # Time steps 

    # Generate the random walk starting with the first value of predefined_series
    random_walk = np.zeros(num_steps)
    random_walk[0] = predefined_series[0]  # Start random walk at the first value of predefined_series

    for i in range(1, num_steps):
        attraction = attract_factor*(0.8 + 2*(np.exp(-i/50))+i*1/168/3) * (predefined_series[i] - random_walk[i - 1])
        noise = noise_strength * np.random.randn()
        random_walk[i] = random_walk[i - 1] + attraction + noise
        if cap_zero:
            random_walk[i] = max(0, random_walk[i])

    return random_walk