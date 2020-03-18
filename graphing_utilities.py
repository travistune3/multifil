import matplotlib.pyplot as plt


def print_constants(constants, print_address=False):
    """prints all settings in an organized fashion"""
    for f_type, filaments in constants.items():
        print(f_type)
        for address, filament in filaments.items():
            address = "\t" + str(address)
            if not print_address:
                address = ""
            print(address, "\t", end="")

            for constant, value in filament.items():
                print(constant, "=", value, end=" ")
            if len(filaments.keys()) < 50:
                print()
            else:
                print(", ", end="\t")


def plot_input_traces(time, length, ap, title=None):
    """plots the experimental traces"""
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(16, 9))
    axes[0].plot(time, length)
    if title is None:
        title = str(len(time)) + " timestep simulation"
    axes[0].set(ylabel='hs length (nm)', title=title)
    axes[1].plot(time, ap)
    axes[1].set(ylabel='actin permissiveness',
                xlabel='millis (ms)')
    plt.tight_layout()
    plt.show()
    
def plot_data(data, l_key='axial_force', r_key='actin_permissiveness',
              title="multifilament model of muscle contraction", save_dir=None):
    
    fs = 12  # font size

    # recreate time trace in milliseconds
    time_trace = data['timestep'].copy()
    for i in range(len(time_trace)):
        time_trace[i] *= data['timestep_length']

    # plot
    fig, axes = plt.subplots(figsize=(8, 4.5))
    axes.plot(time_trace, data[l_key], color='black')

    title = title
    plt.title(title, fontsize=fs*1.5)
    plt.xlabel("time (ms)", fontsize=fs)
    plt.ylabel(l_key, fontsize=fs)

    ax2 = plt.twinx()
    ax2.plot(time_trace, data[r_key]),
    ax2.set(ylabel=r_key)
    
    if save_dir is not None:
        plt.savefig(save_dir + "\\" + data['name'] + ".png", dpi=450)
        print("saved to:", save_dir + "\\" + data['name'] + ".png")
