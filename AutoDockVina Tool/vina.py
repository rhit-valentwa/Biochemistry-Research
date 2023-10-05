# vina.py

print('\n\n\nVina Terminal Tool\n\n')

# Required Modules
import os
import subprocess
import signal
import csv
import statistics
import math
from datetime import datetime

# Ensure necessary software are installed
def checkInstalled(software, args):
    output = ''
    installed = True
    try:
        output = subprocess.run(args, capture_output=True).stdout.decode('utf-8').strip().lower()
    except subprocess.CalledProcessError as e:
        installed = False
    except Exception as e:
        installed = False
    
    if 'command not found' in output or 'error' in output:
        installed = False
    
    if installed == False:
        print(f'The piece of software {software} has not been installed.')
        print(f'Here is a link to help install {software}.')
        print(f'https://google.com/search?q=how+to+install+{software.lower().replace(" ", "+")}')
        print('Closing program.\n')
        os._exit(1)

    print(f'"{software}" is installed.')

checkInstalled('Python module: requests', ['python', '-c', 'import requests; print(True)'])
checkInstalled('Python module: tqdm', ['python', '-c', 'import tqdm; print(True)'])
checkInstalled('Open Babel', ['obabel'])
checkInstalled('AutoDock Vina 1.2', ['vina1.2'])

print('All required software have been installed.\n')

# More Required Modules
import requests
from tqdm import tqdm

date = datetime.now().strftime("%m-%d-%Y")

# Ease-of-Use Methods
def prompt(string):
    output = input(string + '\n')
    if output:
        return (output).lower().strip()
    else:
        return ''

def command(string):
    return os.system(string)
    
def writeFile(file, data):
    with open(file, "w") as file:
        file.write(data)

ligand_names = []
temp_ligand_names = []

def get_ligand(name):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/SDF?record_type=3d'
    response = requests.get(url)
    if response.status_code == 200:
        writeFile('input.sdf', response.text)
        command(f"""
                    obabel input.sdf -O ./data/pdbqt/{name}.pdbqt -p 7.4 >/dev/null 2>&1 
                    rm input.sdf
                """)
        temp_ligand_names.append(name)
    else:
        print(f"Unable to get {name}.")
        

def exit_handler(signal, frame):
    print('\n\nCtrl+C detected. Ending Program.')
    os._exit(1)

signal.signal(signal.SIGINT, exit_handler)

protein_name = prompt('Enter name of protein (default: protein_6.pdbqt):')
ligand_name = prompt('Enter names of ligands (one ligand on each new line):')

if protein_name == '':
    protein_name = 'protein_6.pdbqt'

if ligand_name != '':
    ligand_names.append(ligand_name.replace(' ', '-'))
else:
    print('No ligand entered. Exiting program.')
    os._exit(1)

while True:
    ligand_grammar = 's' if len(ligand_names) > 1 else ''
    output = prompt(f"""You have entered {len(ligand_names)} ligand{ligand_grammar} (blank line to end):""")

    if output == "":
        break
    ligand_names.append(output.replace(' ', '-'))

command("""
        echo ""
        mkdir -p data
        mkdir -p data/pdbqt
        mkdir -p results
        """)

for ligand in ligand_names:
    folder_name = f'{protein_name.replace(".pdbqt", "")}_{ligand}_{date}'
    command(f"""
        mkdir -p "./data/{folder_name}"
        mkdir -p "./data/{folder_name}/Individual Experiments"
        mkdir -p "./data/{folder_name}/Batch Experiments"
        """)
    get_ligand(ligand)

ligand_names = temp_ligand_names

config = {
    'receptor': protein_name,
    'center_x': -10.159,
    'center_y': 3.755,
    'center_z': -11.291,
    'size_x': 82,
    'size_y': 100,
    'size_z': 92,
    'energy_range': 4,
    'exhaustiveness': 64, 
}

if prompt('Would you like to make changes to the configuration? (Y/N, default: N)') == 'y':
    for key, value in config.items():
        output = prompt(f'Input value for {key} (default: {value}):')
        if output:
            config[key] = output

for ligand in tqdm(ligand_names, desc='Docking in progress...'):
    writeFile('config.txt', f"""
    receptor = {os.getcwd()}/{config['receptor']}
    ligand = {os.getcwd()}/data/pdbqt/{ligand}.pdbqt

    center_x = {config['center_x']}
    center_y = {config['center_y']}
    center_z = {config['center_z']}

    size_x = {config['size_x']}
    size_y = {config['size_y']}
    size_z = {config['size_z']}

    energy_range = {config['energy_range']}

    exhaustiveness = {config['exhaustiveness']}
    """)
    folder_name = f'{protein_name.replace(".pdbqt", "")}_{ligand}_{date}/Batch Experiments' if len(ligand_names) > 1 else f'{protein_name.replace(".pdbqt", "")}_{ligand}_{date}/Individual Experiments'

    
    affinities = []

    affinities_mins = []

    for run in tqdm(range(0, 8), desc = f'Docking simulation {protein_name.replace(".pdbqt", "")}-{ligand}'):
        
        current_affinities = []
        
        command(f"""vina1.2 --config config.txt --out "./data/{folder_name}/{protein_name.replace(".pdbqt", "")}_{ligand}_run-{run + 1}_output.pdbqt" > log.txt""")

        raw_data = open('log.txt', 'r').read().split('-----+------------+----------+----------')[1].split('Writing')[0]
        log_data = [list(map(float, row.split())) for row in raw_data.strip().split('\n')]

        best_docking = float('-inf')
        for row in log_data:
            current_affinities.append(float(row[1]))
        
        affinities_mins.append(min(current_affinities))
        best_docking = min(affinities_mins)

        raw_data = """mode |   affinity | dist from best mode | (kcal/mol) | rmsd l.b.| rmsd u.b.\n""" + raw_data

        with open(f'./results/{protein_name.replace(".pdbqt", "")}_ligands-{len(ligand_names)}_{date}_results.csv', 'a', newline='') as file:
            writer = csv.writer(file)
            if ligand == ligand_names[0] and run == 0:
                writer.writerow(['protein', 'ligand', 'run', 'center_x', 'center_y', 'center_z', 'size_x', 'size_y', 'size_z', 'energy_range', 'exhaustiveness', 'log_results', 'best_docking'])
            
            writer.writerow([protein_name.replace(".pdbqt", ""), ligand, run + 1, config['center_x'], config['center_y'], config['center_z'], 
                             config['size_x'], config['size_y'], config['size_z'], config['energy_range'], 
                             config['exhaustiveness'], raw_data, best_docking])

    with open(f'./results/{protein_name.replace(".pdbqt", "")}_ligands-{len(ligand_names)}_{date}_Analysed_results.csv', 'a', newline='') as file:
            writer = csv.writer(file)
            if sum(1 for _ in open(f'./results/{protein_name.replace(".pdbqt", "")}_ligands-{len(ligand_names)}_{date}_Analysed_results.csv')) == 0:
                writer.writerow(['receptor', 'ligand', 'max', 'min', 'avg', 'st_dev', 'CI_lower_95', 'CI_upper_95'])

            lower_bound = (sum(affinities_mins) / len(affinities_mins)) - (2.262 * (statistics.stdev(affinities_mins) / math.sqrt(len(affinities_mins))))
            upper_bound = (sum(affinities_mins) / len(affinities_mins)) + (2.262 * (statistics.stdev(affinities_mins) / math.sqrt(len(affinities_mins))))
            writer.writerow([protein_name.replace(".pdbqt", ""), ligand, max(affinities_mins), min(affinities_mins), sum(affinities_mins) / len(affinities_mins), statistics.stdev(affinities_mins), lower_bound, upper_bound])

# Deleting undesired files
command("""rm config.txt
           rm log.txt""")