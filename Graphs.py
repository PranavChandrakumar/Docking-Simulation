import numpy as np
import matplotlib.pyplot as plt

v_x = []
v_y = []
v_z = []
r = []
v = []
t = []

def main():
    #Plots for each variable
    fig, axs = plt.subplots(5, 1, figsize=(10, 25))

    # Plot v_x against t
    axs[0].plot(t, v_x, label='v_x')
    axs[0].set_xlabel('t')
    axs[0].set_ylabel('v_x')
    axs[0].set_title('v_x vs t')
    axs[0].legend()
    axs[0].grid(True)

    # Plot v_y against t
    axs[1].plot(t, v_y, label='v_y')
    axs[1].set_xlabel('t')
    axs[1].set_ylabel('v_y')
    axs[1].set_title('v_y vs t')
    axs[1].legend()
    axs[1].grid(True)

    # Plot v_z against t
    axs[2].plot(t, v_z, label='v_z')
    axs[2].set_xlabel('t')
    axs[2].set_ylabel('v_z')
    axs[2].set_title('v_z vs t')
    axs[2].legend()
    axs[2].grid(True)

    # Plot r against t
    axs[3].plot(t, r, label='r')
    axs[3].set_xlabel('t')
    axs[3].set_ylabel('r')
    axs[3].set_title('r vs t')
    axs[3].legend()
    axs[3].grid(True)

    # Plot v against t
    axs[4].plot(t, v, label='v')
    axs[4].set_xlabel('t')
    axs[4].set_ylabel('v')
    axs[4].set_title('v vs t')
    axs[4].legend()
    axs[4].grid(True)

    plt.tight_layout()
    
    #Display plots
    plt.show()

    years = int(np.floor(t[-1]//3.154e7))
    months = int(np.floor((t[-1]-years*3.154e7)/2.62e6))
    weeks = int(np.floor((t[-1]-years*3.154e7-months*2.62e6)/604800))
    days = int(np.floor((t[-1]-years*3.154e7-months*2.62e6-weeks*604800)/86400))
    hours = int(np.floor((t[-1]-years*3.154e7-months*2.62e6-weeks*604800-days*86400)/3600))
    minutes = int(np.floor((t[-1]-years*3.154e7-months*2.62e6-weeks*6048800-days*86400-hours*3600)/60))
    seconds = int((t[-1]-years*3.154e7-months*2.62e6-weeks*6048800-days*86400-hours*3600-minutes*60))
    print("Total Simulation Time:", years,"years",months,"months",weeks,"weeks",days,"days",hours,"hours",minutes,"minutes",seconds,"seconds")
    return