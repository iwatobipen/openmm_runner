import yaml
from collections import OrderedDict


config = {
    'pdb_file':'',
    'resname':'',
    'lig_smiles':'',
    'ignore_missing_residues':False,
    'ignore_terminal_missing_residues':False,
    'ph':7.0,
    'data_path': 'data',
    'protein_ff':"amber14-all.xml",
    'solvent_ff':"amber14/tip3pfb.xml",
    'integrator' : 'LangevinIntegrator',
    'integrator_settings':  {'temperature':300,
                            'frictionCoeff':1.0,
                            'stepSize':2.0},
    'platform': {'name':'CUDA',
                 'properties': {'DeviceIndex':'0'}
                },
    'md_steps' : 1000,
    'write_interval': 1,
    'log_interval': 1
          }

def makeconfig():
    with open('config_template.yml', 'w') as writer:
        yaml.dump(config, writer,default_flow_style=False, sort_keys=False)

if __name__=='__main__':
    makeconfig()