# -*- coding: utf-8 -*-

# Global imports

from __future__ import absolute_import
import os,re,sys
import glob
import math
import mdtraj as md
import multiprocessing as mp
from functools import partial


# Script information
__author__ = "Sergi Roda"
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Sergi Roda"
__email__ = "sergi.rodallordes@bsc.es"

def GetResults(trajectory_list, n_processors,atom_index,residue, metric, res_name, topology):

    Metric_results = []
    print(trajectory_list, n_processors,atom_index,residue, metric, res_name, topology)
    pool = mp.Pool(n_processors)
    if len(atom_index) == 2 or len(residue) == 2:
        AddDistance_func = partial(AddDistance, atom_index = atom_index, residue = residue, metric = metric, res_name = res_name, topology = topology)
        Metric_results.append(pool.map(AddDistance_func,trajectory_list))
    elif len(atom_index) == 3 or len(residue) == 3:
        AddAngle_func = partial(AddAngle, atom_index = atom_index, residue = residue, metric = metric, res_name = res_name, topology = topology)
        Metric_results.append(pool.map(AddAngle_func,trajectory_list))
    elif len(atom_index) == 4 or len(residue) == 4:
        AddDihedral_func = partial(AddDihedral, atom_index = atom_index, residue = residue, metric = metric, res_name = res_name, topology = topology)
        Metric_results.append(pool.map(AddDihedral_func,trajectory_list))
    pool.terminate()

    return Metric_results

def GetLigResults(trajectory_list, n_processors,atom_index,residue, metric, res_name, topology):

    Metric_results = []

    pool = mp.Pool(n_processors)
    if len(atom_index) == 2 or len(residue) == 2:
        AddLigDistance_func = partial(AddLigDistance, atom_index = atom_index, residue = residue, metric = metric, res_name = res_name, topology = topology)
        Metric_results.append(pool.map(AddLigDistance_func,trajectory_list))
    elif len(atom_index) == 3 or len(residue) == 3:
        AddAngle_func = partial(AddAngle, atom_index = atom_index, residue = residue, metric = metric, res_name = res_name, topology = topology)
        Metric_results.append(pool.map(AddAngle_func,trajectory_list))
    elif len(atom_index) == 4 or len(residue) == 4:
        AddDihedral_func = partial(AddDihedral, atom_index = atom_index, residue = residue, metric = metric, res_name = res_name, topology = topology)
        Metric_results.append(pool.map(AddDihedral_func,trajectory_list))
    pool.terminate()

    return Metric_results

def AddLigDistance(trajectory, atom_index, residue, metric, res_name, topology):
    """
    Take the PELE simulation trajectory files and returns the list of values of the desired distance metric

    RETURNS
    -------
    metric_list: list
            List of values of the desired distance metric
    """

    if ".xtc" in trajectory:
        traj = md.load_xtc(trajectory, topology)
    else:
        traj = md.load_pdb(trajectory)

    if len(atom_index) != 0:
        metric_list = 10*md.compute_distances(traj,[[int(atom_index[0]),int(atom_index[1])]])
    else:
        Atom_pair_1 = int(traj.topology.select("name {} and resname {}".format(metric[0].strip("_"),res_name[0])))
        Atom_pair_2 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(residue[1],metric[1].strip("_"),res_name[1])))
        metric_list = 10*md.compute_distances(traj,[[Atom_pair_1,Atom_pair_2]])

    return metric_list

def AddDistance(trajectory, atom_index, residue, metric, res_name, topology):
    """
    Take the PELE simulation trajectory files and returns the list of values of the desired distance metric

    RETURNS
    -------
    metric_list: list
            List of values of the desired distance metric
    """

    if ".xtc" in trajectory:
        traj = md.load_xtc(trajectory, topology)
    else:
        traj = md.load_pdb(trajectory)

    if len(atom_index) != 0:
        metric_list = 10*md.compute_distances(traj,[[int(atom_index[0]),int(atom_index[1])]])
    else:
        Atom_pair_1 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(residue[0],metric[0].strip("_"),res_name[0])))
        Atom_pair_2 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(residue[1],metric[1].strip("_"),res_name[1])))
        print(Atom_pair_1, Atom_pair_2)
        metric_list = 10*md.compute_distances(traj,[[Atom_pair_1,Atom_pair_2]])

    return metric_list

def AddAngle(trajectory, atom_index, residue, metric, res_name, topology):
    """
    Take the PELE simulation trajectory files and returns the list of values of the desired angle metric

    RETURNS
    -------
    metric_list: list
            List of values of the desired angle metric
    """

    if ".xtc" in trajectory:
        traj = md.load_xtc(trajectory, topology)
    else:
        traj = md.load_pdb(trajectory)

    if len(atom_index) != 0:
        metric_list = md.compute_angles(traj,[[int(atom_index[0]),int(atom_index[1]), int(atom_index[2])]])
    else:
        Atom_pair_1 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(residue[0],metric[0].strip("_"),res_name[0])))
        Atom_pair_2 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(residue[1],metric[1].strip("_"),res_name[1])))
        Atom_pair_3 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(residue[2],metric[2].strip("_"),res_name[2])))
        metric_list = md.compute_angles(traj,[[Atom_pair_1,Atom_pair_2, Atom_pair_3]])

    return metric_list

def AddDihedral(trajectory, atom_index, residue, metric, res_name, topology):
    """
    Take the PELE simulation trajectory files and returns the list of values of the desired dihedral metric

    RETURNS
    -------
    metric_list: list
            List of values of the desired dihedral metric
    """

    if ".xtc" in trajectory:
        traj = md.load_xtc(trajectory, topology)
    else:
        traj = md.load_pdb(trajectory)

    if len(atom_index) != 0:
        metric_list = md.compute_dihedrals(traj,[[int(atom_index[0]),int(atom_index[1]), int(atom_index[2]), int(atom_index[3])]])
    else:
        Atom_pair_1 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(residue[0],metric[0].strip("_"),res_name[0])))
        Atom_pair_2 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(residue[1],metric[1].strip("_"),res_name[1])))
        Atom_pair_3 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(residue[2],metric[2].strip("_"),res_name[2])))
        Atom_pair_4 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(residue[3],metric[3].strip("_"),res_name[3])))
        metric_list = md.compute_dihedrals(traj,[[Atom_pair_1,Atom_pair_2, Atom_pair_3, Atom_pair_4]])

    return metric_list

class MetricCalculator(): 
    
    def __init__(self, Input, list_of_residue_numbers, list_of_atom_names, list_of_residue_names, column_name, report_format, trajectory_format, topology, output_name, n_processors, atom_index = "", verbose = False):

        self.input, self.atom_index, self.residue, self.metric, self.res_name, self.column_name, \
        self.report_format, self.trajectory_format, self.topology, \
        self.output_name, self.verbose, self.n_processors = Input, atom_index, list_of_residue_numbers, list_of_atom_names, \
        list_of_residue_names, column_name, report_format, trajectory_format, topology, output_name, verbose, n_processors

    def DecompressList(self, l_of_lists):
        """
        This method decompress a 
        list of lists into a list

        PARAMETERS
        ----------
        l_of_lists: list of lists
                List of lists

        RETURNS
        -------
        new_list: list
                New list    
        """

        new_list = []
        for sublist in l_of_lists:
            for item in sublist:
                new_list.append(item)

        return new_list

    def AddDistance(self, trajectory):
        """
        Take the PELE simulation trajectory files and returns the list of values of the desired distance metric

        RETURNS
        -------
        metric_list: list
                List of values of the desired distance metric
        """

        if ".xtc" in trajectory:
            traj = md.load_xtc(trajectory, self.topology)
        else:
            traj = md.load_pdb(trajectory)

        if len(self.atom_index) != 0:
            metric_list = 10*md.compute_distances(traj,[[int(self.atom_index[0]),int(self.atom_index[1])]])
        else:
            Atom_pair_1 = int(traj.topology.select("name {} and resname {}".format(self.metric[0].strip("_"),self.res_name[0])))
            Atom_pair_2 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[1],self.metric[1].strip("_"),self.res_name[1])))
            metric_list = 10*md.compute_distances(traj,[[Atom_pair_1,Atom_pair_2]])

        return metric_list

    def AddAngle(self, trajectory):
        """
        Take the PELE simulation trajectory files and returns the list of values of the desired angle metric

        RETURNS
        -------
        metric_list: list
                List of values of the desired angle metric
        """

        if ".xtc" in trajectory:
            traj = md.load_xtc(trajectory, self.topology)
        else:
            traj = md.load_pdb(trajectory)

        if len(self.atom_index) != 0:
            metric_list = md.compute_angles(traj,[[int(self.atom_index[0]),int(self.atom_index[1]), int(self.atom_index[2])]])
        else:
            Atom_pair_1 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[0],self.metric[0].strip("_"),self.res_name[0])))
            Atom_pair_2 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[1],self.metric[1].strip("_"),self.res_name[1])))
            Atom_pair_3 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[2],self.metric[2].strip("_"),self.res_name[2])))
            metric_list = md.compute_angles(traj,[[Atom_pair_1,Atom_pair_2, Atom_pair_3]])

        return metric_list

    def AddDihedral(self, trajectory):
        """
        Take the PELE simulation trajectory files and returns the list of values of the desired dihedral metric

        RETURNS
        -------
        metric_list: list
                List of values of the desired dihedral metric
        """

        if ".xtc" in trajectory:
            traj = md.load_xtc(trajectory, self.topology)
        else:
            traj = md.load_pdb(trajectory)

        if len(self.atom_index) != 0:
            metric_list = md.compute_dihedrals(traj,[[int(self.atom_index[0]),int(self.atom_index[1]), int(self.atom_index[2]), int(self.atom_index[3])]])
        else:
            Atom_pair_1 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[0],self.metric[0].strip("_"),self.res_name[0])))
            Atom_pair_2 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[1],self.metric[1].strip("_"),self.res_name[1])))
            Atom_pair_3 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[2],self.metric[2].strip("_"),self.res_name[2])))
            Atom_pair_4 = int(traj.topology.select("resSeq {} and name {} and resname {}".format(self.residue[3],self.metric[3].strip("_"),self.res_name[3])))
            metric_list = md.compute_dihedrals(traj,[[Atom_pair_1,Atom_pair_2, Atom_pair_3, Atom_pair_4]])

        return metric_list

    def AddMetric(self):
        """
        Take the PELE simulation trajectory files and returns the report files with the desired metric
        
        OUTPUT
        ------
        Report files with the desired metric added.
        """
     
        def atoi(text):
            return int(text) if text.isdigit() else text
        
        def natural_keys(text):
            return [atoi(c) for c in re.split(r'(\d+)', text)]

        trajectory_list, report_list = [],[]


        for path in glob.glob(self.input): 
            trajectory_list += glob.glob(os.path.join(path,"*"+self.trajectory_format))
            report_list += glob.glob(os.path.join(path, self.report_format))
            
        report_list.sort(key=natural_keys)
        trajectory_list.sort(key=natural_keys)

        Metric_results = GetResults(trajectory_list, self.n_processors, self.atom_index, self.residue, self.metric, self.res_name, self.topology)

        Metric_results = self.DecompressList(Metric_results)

        j = -1
        for report in report_list:
            if self.verbose:
                print(report)

            j+=1

            if self.verbose:
                print(j)

            with open(report, 'r') as report_file:
                if report.endswith(".out"):
                    out_report = open(report.split(".out")[0]+self.output_name,'w')
                else:
                    out_report = open(report+self.output_name,'w')
                for i,line in enumerate(report_file):
                    if i==0:
                        out_report.write(line.strip("\n")+'    '+self.column_name+"\n")
                    else:
                        if len(self.atom_index) > 2 or len(self.residue) > 2:
                            out_report.write(line.strip("\n")+'    '+str(math.degrees(Metric_results[j][i-1][0]))+"\n")
                        else:
                            out_report.write(line.strip("\n")+'    '+str(Metric_results[j][i-1][0])+"\n")

        if self.verbose:
            print("{} report files have the desired metric added".format(j+1))

    def AddLigMetric(self):
        """
        Take the PELE simulation trajectory files and returns the report files with the desired metric
        
        OUTPUT
        ------
        Report files with the desired metric added.
        """

        def atoi(text):
            return int(text) if text.isdigit() else text

        def natural_keys(text):
            return [atoi(c) for c in re.split(r'(\d+)', text)]

        trajectory_list, report_list = [],[]


        for path in glob.glob(self.input):
            trajectory_list += glob.glob(os.path.join(path,"*"+self.trajectory_format))
            report_list += glob.glob(os.path.join(path, self.report_format))

        report_list.sort(key=natural_keys)
        trajectory_list.sort(key=natural_keys)

        Metric_results = GetLigResults(trajectory_list, self.n_processors, self.atom_index, self.residue, self.metric, self.res_name, self.topology)

        Metric_results = self.DecompressList(Metric_results)

        j = -1
        for report in report_list:
            if self.verbose:
                print(report)

            j+=1

            if self.verbose:
                print(j)

            with open(report, 'r') as report_file:
                if report.endswith(".out"):
                    out_report = open(report.split(".out")[0]+self.output_name,'w')
                else:
                    out_report = open(report+self.output_name,'w')
                for i,line in enumerate(report_file):
                    if i==0:
                        out_report.write(line.strip("\n")+'    '+self.column_name+"\n")
                    else:
                        if len(self.atom_index) > 2 or len(self.residue) > 2:
                            out_report.write(line.strip("\n")+'    '+str(math.degrees(Metric_results[j][i-1][0]))+"\n")
                        else:
                            try:
                                out_report.write(line.strip("\n")+'    '+str(Metric_results[j][i-1][0])+"\n")
                            except IndexError:
                                import ipdb; ipdb.set_trace()

        if self.verbose:
            print("{} report files have the desired metric added".format(j+1))

if __name__ == "__main__":
#    #"""Call the main function"""
    Metriccalculator = MetricCalculator(sys.argv[1],[int(i) for i in sys.argv[2:4]],sys.argv[4:6],sys.argv[6:8],sys.argv[8],sys.argv[9],sys.argv[10],sys.argv[11],sys.argv[12],int(sys.argv[13]))
    Metriccalculator.AddMetric()
