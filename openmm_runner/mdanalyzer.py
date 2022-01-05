import pandas as pd
import numpy as np
import nglview as nv
import MDAnalysis as mda
from MDAnalysis.analysis import rms, diffusionmap, align
from MDAnalysis.analysis.distances import dist
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA


def align_trajectroy(data_path, topologyfile='topology.pdb', trajectroyfile='trajectrory.xtc'):
    md_universe = mda.Universe(str(data_path / topologyfile), str(data_path / trajectroyfile))
    # Set trajectory pointer to first frame as reference
    md_universe.trajectory[0]
    
    # in_memory=True is needed to actually modify the loaded trajectory
    # if the trajectory is too big for memory, one can write out the aligned trajectory and reload it into a new universe
    alignment = align.AlignTraj(
        mobile=md_universe, reference=md_universe, select="protein", in_memory=True
    )
    alignment.run()
    return md_universe


def rmsd_for_atomgroups(universe, selection1, selection2=None):
    """Calulate the RMSD for selected atom groups.

    Parameters
    ----------
    universe: MDAnalysis.core.universe.Universe
        MDAnalysis universe.
    selection1: str
        Selection string for main atom group, also used during alignment.
    selection2: list of str, optional
        Selection strings for additional atom groups.

    Returns
    -------
    rmsd_df: pandas.core.frame.DataFrame
        DataFrame containing RMSD of the selected atom groups over time.
    """

    universe.trajectory[0]
    ref = universe
    rmsd_analysis = rms.RMSD(universe, ref, select=selection1, groupselections=selection2)
    rmsd_analysis.run()
    columns = [selection1, *selection2] if selection2 else [selection1]
    rmsd_df = pd.DataFrame(np.round(rmsd_analysis.rmsd[:, 2:], 2), columns=columns)
    rmsd_df.index.name = "frame"
    return rmsd_df


