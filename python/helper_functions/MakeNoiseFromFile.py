import numpy as np
import matplotlib.pyplot as plt

def make_noise_from_file(AS):
    t = 0
    n = 0
    tmp = 0
    x = []
    y = []

    while t == 0:
        n += 1
        SS = slice((n - 1) * 1000, n * 1000)
        plt.figure(1)
        plt.plot(AS[SS])
        points = plt.ginput(2, show_clicks=True)

        if not points:
            break

        t1, t2 = zip(*points)
        x.append(t1[0])
        y.append(t2[0])
        tmp += t2[0] - t1[0]

        SSS = slice((n - 1) * 1000 + 1 + round(x[n - 1]), (n - 1) * 1000 + 1 + round(y[n - 1]))
        plt.figure(2)
        plt.plot(AS[SSS])
        plt.title(str(tmp))

    # Anneal bits
    ASnoise = np.concatenate([AS[slice((n - 1) * 1000 + 1 + round(x[n - 1]), (n - 1) * 1000 + round(y[n - 1]))] for n in range(1, len(x))])

    plt.figure()
    plt.plot(ASnoise)
    plt.hold(True)
    plt.plot(AS)

    return ASnoise

# Example usage:
# Replace AS with your actual data
# AS = your_AS_value
# ASnoise_result = make_noise_from_file(AS)
# plt.show()  # Display the plots
