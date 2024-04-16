import csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def load_C_mean(filename):
    with open(filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            C = np.array([float(i) for i in row])
    return C

def load_C_list(filename):
    C_list = []
    with open(filename, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            C = np.array([float(i) for i in row])
            C_list.append(C)
    csvfile.close()
    return C_list


# cosine similarity between vectors
def cos_sim(v1,v2):
    return np.sum(v1*v2) / (np.sqrt((np.sum(v1*v1))*np.sum(v2*v2)))


def plot_sim(sim, mu, xy_range, ax_label, save, figname):
    fig,ax = plt.subplots(1,1,figsize=(4.5,4))

    minValue = min([min(sim[n][np.nonzero(sim[n])]) for n in range(len(sim)-1)])
    maxValue = max([max(sim[n]) for n in range(len(sim))])


    plt.rcParams.update({'font.size': 13})

    sns.heatmap(sim.transpose(),
                xticklabels=xy_range,
                yticklabels=xy_range,
                annot=False,
                mask=np.triu(np.ones_like(sim, dtype=bool)),
                vmin=minValue,
                vmax=maxValue,
                cmap="viridis")
    plt.xlabel(ax_label)
    plt.ylabel(ax_label)
    plt.tight_layout()
    if save:
        figname += '.pdf'
        print(figname)
        plt.savefig(figname)
    plt.show()


def plot_sim_different_axes(sim, x_range, y_range, x_label, y_label, figname, save):
    fig,ax = plt.subplots(1,1,figsize=(4.5,4))
    maxValue = max([max(sim[n]) for n in range(len(sim))])
    minValue = min([min(sim[n][np.nonzero(sim[n])]) for n in range(len(sim)-1)])
    print(minValue,maxValue)
    plt.rcParams.update({'font.size': 13})
    sns.heatmap(sim,
                xticklabels=x_range,
                yticklabels=y_range,
                annot=False,
                vmin=minValue,
                vmax=maxValue,
                cmap="viridis")
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    #plt.title(title)
    #ax.set_title(title)
    plt.tight_layout()
    if save:
        plt.savefig(figname)
    plt.show()


def save_on_csv(filename, variable_list, writing_operation):
    with open(filename, writing_operation) as csvfile:
        writer = csv.writer(csvfile)
        try:
            [writer.writerow(s) for s in variable_list]
        except:
            writer.writerow(variable_list)
