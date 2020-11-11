from glob import glob
import pandas as pd
import seaborn as sns

class Analysis:
    def __init__(self, folder):
        """folder(str): path to the simulation folder"""
        self.folder = folder
        self.dataframe = None
        self.distance = None
        self.distribution = None

    def filtering(self):
        reports = []
        for files in glob("{}/output/0/report_*".format(self.folder)):
            reports.append(pd.read_csv(files, sep="    ", engine="python"))

        self.dataframe = pd.concat(reports)
        self.dataframe.sort_values(by="Binding Energy", inplace=True)
        self.dataframe.reset_index(drop=True, inplace=True)
        self.dataframe = self.dataframe.head(len(self.dataframe)*10/100)
        self.distance = self.dataframe["distance0.5"].copy()
        self.distance.sort_values(inplace=True)
        self.distance.reset_index(drop=True, inplace=True)
        self.distance = self.distance.head(min(50, len(self.distance)))

        if "original" in self.folder:
            self.distance = self.distance[0].copy()

        return self.dataframe, self.distance

    def set_distribution(self, num):
        self.distribution = self.distance - num

def analyse_all(folders="PELE_*"):
    """folders (str): path to the different PELE simulation folders to be analyzed"""
    dictio = {}
    for folder in glob("{}".format(folders)):
        simulation = Analysis(folder)
        simulation.filtering()
        dictio[folder[5:]] = simulation
        if "original" in folder:
            distance = simulation.distance
        else:
            distance = 3

    plot_dictio = {}
    for key, value in dictio.items():
        if "original" not in key:
            value.set_distribution(distance)
            plot_dictio[key] = value.distribution

    plot_dataframe = pd.DataFrame(plot_dictio)
    plot_dataframe.columns.name = "mutations"

    return dictio, plot_dataframe

def box_plot(dataframe):
    """dataframe (tabular data): Any tabular structure that corresponds to the wide-form data accepted by seaborn"""
    ax = sns.catplot(data=dataframe, kind="box", palette="Set2")
    ax.set(xlabel="Mutations", ylabel="Distance", title="Distance variantion with respect to original")
    ax.savefig("distance.pdf" )

    return ax
