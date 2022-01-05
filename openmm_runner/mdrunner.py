import yaml
import sys
import os
import argparse
from openmm.app import PDBFile
from . import mdutils
from . import mdanalyzer
import mdtraj as md
import openmm as mm
import openmm.app as app
from openmm import unit
from openff.toolkit.topology import Molecule, Topology
from openmmforcefields.generators import GAFFTemplateGenerator
from pathlib import Path


def getParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('configs')
    return parser


def generate_forcefield(
    rdkit_mol=None, protein_ff="amber14-all.xml", solvent_ff="amber14/tip3pfb.xml"
):
    """
    Generate an OpenMM Forcefield object and register a small molecule.

    Parameters
    ----------
    rdkit_mol: rdkit.Chem.rdchem.Mol
        Small molecule to register in the force field.
    protein_ff: string
        Name of the force field.
    solvent_ff: string
        Name of the solvent force field.

    Returns
    -------
    forcefield: simtk.openmm.app.Forcefield
        Forcefield with registered small molecule.
    """
    forcefield = app.ForceField(protein_ff, solvent_ff)

    if rdkit_mol is not None:
        gaff = GAFFTemplateGenerator(
            molecules=Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
        )
        forcefield.registerTemplateGenerator(gaff.generator)

    return forcefield


def runmd():
    print('load config')
    config = getParser().parse_args()
    with open(config.configs, 'r') as file:
        md_config = yaml.load(file)
    
    data = md_config['data_path']
    #HERE = Path(_dh[-1])
    HERE = Path(os.getcwd())
    DATA = HERE / data
    DATA.mkdir(exist_ok=True)

    print('Prepare Protein....')
    prepared_protein = mdutils.prepare_protein(md_config['pdb_file'],
    ignore_missing_residues=md_config['ignore_missing_residues'],
    ignore_terminal_missing_residues=md_config['ignore_terminal_missing_residues'],
    ph=md_config['ph'])

    print('Prepare Ligand....')
    rdkit_ligand = mdutils.prepare_ligand(md_config['pdb_file'],
    resname=md_config['resname'],
    smiles=md_config['lig_smiles'])

    omm_ligand = mdutils.rdkit_to_openmm(rdkit_ligand, md_config['resname'])
    
    print('Merge Protein and Ligand....')
    complex_topology, complex_positions = mdutils.merge_protein_and_ligand(prepared_protein, omm_ligand)
    print("Complex topology has", complex_topology.getNumAtoms(), "atoms.")

    forcefield = generate_forcefield(rdkit_ligand, md_config['protein_ff'], md_config['solvent_ff'])

    modeller = app.Modeller(complex_topology, complex_positions)
    modeller.addSolvent(forcefield, padding=1.0 * unit.nanometers, ionicStrength=0.15 * unit.molar)

    print('Setup simulation...')
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME)
    
    # this part should be refactaring!!!
    integrator = mm.LangevinIntegrator(
        md_config['integrator_settings']['temperature'] * unit.kelvin, 
        md_config['integrator_settings']['frictionCoeff'] / unit.picoseconds, 
        md_config['integrator_settings']['stepSize']  * unit.femtoseconds
    )
    if md_config['platform'] != "":
        platform = mm.Platform.getPlatformByName(md_config['platform']['name'])
        properties = md_config['platform']['properties']
    else:
        platform = None
        properties = None
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)
    
    print('Run MD')
    simulation.minimizeEnergy()
    with open(DATA / "topology.pdb", "w") as pdb_file:
        app.PDBFile.writeFile(
            simulation.topology,
            simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(),
            file=pdb_file,
            keepIds=True,
        )

    simulation.reporters.append(
    md.reporters.XTCReporter(file=str(DATA / "trajectory.xtc"), reportInterval=md_config['write_interval'])
    )
    simulation.reporters.append(
        app.StateDataReporter(
            sys.stdout,
            md_config['log_interval'],
            step=True,
            potentialEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=md_config['md_steps'],
            separator="\t",
        )
    )
    
    simulation.context.setVelocitiesToTemperature(md_config['integrator_settings']['temperature'] * unit.kelvin)
    simulation.step(md_config['md_steps'])  # perform the simulation
    (DATA / "trajectory.xtc").stat().st_size > 0


if __name__=='__main__':
    runmd()