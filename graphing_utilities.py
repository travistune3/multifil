import matplotlib.pyplot as plt

def print_constants(constants, print_address=False):
    """prints all settings in an organized fashion"""
    for f_type, filaments in constants.items():
        print(f_type)
        for address, filament in filaments.items():
            address = "\t" + str(address)
            if not print_address:
                address = ""
            print(address, "\t", end = "")
            
            for constant, value in filament.items():
                print(constant, "=", value, end=" ")
            if len(filaments.keys()) < 50:
                print()
            else:
                print(", ", end="\t")
                
def plot_input_traces(time, length, ap):
    """plots the experimental traces"""
    fig, axes = plt.subplots(2,1,sharex=True, figsize=(16,9))
    axes[0].plot(time, length)
    axes[0].set(ylabel='hs length (nm)', title=str(len(time)) + " timestep simulation")
    axes[1].plot(time, ap)
    axes[1].set(ylabel='actin permissiveness',
               xlabel='time (ms)')
    plt.tight_layout()
    plt.show()