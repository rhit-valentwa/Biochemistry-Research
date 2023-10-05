import subprocess
import csv
import os
import signal
from scipy import signal as scipy_signal
import numpy
from pathlib import Path
from datetime import datetime

base_ions = ['a', 'b', 'c', 'x', 'y', 'z']
modifications = ['', '-NH3', '-H2O', '+H2O']

n_ions = [f'{ion}{mod}' for ion in base_ions[:3] for mod in modifications]
c_ions = [f'{ion}{mod}' for ion in base_ions[3:] for mod in modifications]
row_names = []

failed_matches = [[], [], [], []]

color_script = ''

def exit_handler(signal, frame):
    print('\n\nCtrl+C detected. Ending Program.')
    os._exit(1)

signal.signal(signal.SIGINT, exit_handler)

def prompt(string):
    return input(string + '\n').strip().lower() or ''

def determineTerminal(name):
    if name in n_ions:
        return 'n'
    elif name in c_ions:
        return 'c'

def prepareData(file_name):
    with open(f'./{file_name}') as file:

        terminals = [[], [], [], []]
        data = []
        results = []
        intensities = []

        reader = csv.DictReader(file, delimiter="\t")
        for line in reader:
            parsed_line = [
                float(line['m/z']) if line['m/z'] else None,
                int(line['z']) if line['z'] else None,
                float(line['Intensity']) if line['Intensity'] else None,
                int(line['Peptide #']) if line['Peptide #'] else None,
                str(line['Ion Type']) if line['Ion Type'] else None,
                int(line['Index']) if line['Index'] else None,
                int(line['Charge']) if line['Charge'] else None,
                str(line['Elem Comp']) if line['Elem Comp'] else None,
                float(line['Error']) if line['Error'] else None,
            ]

            if line['Intensity'] != None and line['Intensity'] != '' and float(line['Intensity']) != 0:
                intensities.append(float(line['Intensity']))
                data.append(parsed_line)

        y_coordinates = numpy.array(intensities)
        peak_indices = scipy_signal.find_peaks_cwt(y_coordinates, numpy.arange(1, 2))

        for peak in peak_indices:
            results.append(data[peak])

        # remove duplicates
        # results = set(tuple(result) for result in results)

        for row in results:
            if row[3] != None and row[4] != None:
                terminal = determineTerminal(row[4])
                peptide_num = row[3]
                if terminal == 'n':
                    if peptide_num == 1:
                        terminals[0].append(row)
                    elif peptide_num == 2:
                        terminals[1].append(row)
                elif terminal == 'c':
                    if peptide_num == 1:
                        terminals[2].append(row)
                    elif peptide_num == 2:
                        terminals[3].append(row)

        return [numpy.array(sublist) for sublist in terminals]

def processData(terminals, sequence):
    output = [[], [], []]

    fragment_label_fractions = [['n-terminal labeled fraction'],
                                ['n-terminal denominator'],
                                ['c-terminal labeled fraction'],
                                ['c-terminal denominator'],
                                ['flipped c-terminal fraction']]

    output[0].extend([''] + ['n-terminal'] + ([''] * (len(n_ions)) * 2) + ['c-terminal'])
    output[1].extend(([''] + ['unlabeled'] + ([''] * (len(n_ions) - 1)) + ['labeled'] + ([''] * (len(n_ions) - 1))) * 2)
    output[2].extend([''] + (n_ions * 2) + [''] + (c_ions * 2))

    for i in range(0, len(failed_matches)):
        for j in range(0, len(sequence)):
            failed_matches[i].append({})
            for ion in (n_ions if (i == 0 or i == 1) else c_ions):
                failed_matches[i][len(failed_matches[i]) - 1][ion] = False

    for index in range(0, len(sequence)):
        output.append([f'{index + 1} ({sequence[index]})'])

    for l in range(len(terminals)):
        i = 0
        for ion in (n_ions if (l <= 1) else c_ions):
            if l == 2 and ion == c_ions[0]:
                for j in range(len(sequence), 0, -1):
                    output[i + 3].append(f'{j} ({sequence[i]})')
                    i += 1
            for a in range(len(sequence)) if (l <= 1) else range(len(sequence) - 1, 0, -1):
                value = 0
                for row in terminals[l]:
                    if row[4] == ion and int(row[5]) == a:
                        value += float(row[2])
                
                if value != 0 and failed_matches[l][a][ion] == False:
                    output[(a + 3) if (l <= 1) else a].append(value)
                else:
                    output[(a + 3) if (l <= 1) else a].append('')
                    failed_matches[l][a][ion] = True

    return output

def sheetGenerator(data, sequence):
    global color_script
    global row_names
    # prepare data storage arrays
    sheet = [[], [], [], [], [], [], [], [], [], [], []]

    i = 0
    for row in data[3:]:
        n_terminal_unlabeled = [x for x in row[1:len(n_ions)] if x != '']
        n_terminal_labeled = [x for x in row[len(n_ions): len(n_ions) * 2 + 1] if x != '']

        c_terminal_unlabeled = [x for x in row[len(n_ions) * 2 + 2: len(n_ions) * 2 + 1 + len(c_ions)] if x != '']
        c_terminal_labeled = [x for x in row[len(n_ions) * 2 + 2 + len(c_ions):] if x != '']

        sheet[0].append(sum(n_terminal_unlabeled))
        sheet[1].append(sum(n_terminal_labeled))

        if sum(n_terminal_unlabeled) + sum(n_terminal_labeled) == 0:
            sheet[2].append(-1)
        else:
            sheet[2].append(sum(n_terminal_unlabeled) / (sum(n_terminal_unlabeled) + sum(n_terminal_labeled)))

        sheet[3].append(sum(c_terminal_unlabeled))
        sheet[4].append(sum(c_terminal_labeled))

        if sum(c_terminal_unlabeled) + sum(c_terminal_labeled) == 0:
            sheet[5].append(-1)
        else:
             sheet[5].append(sum(c_terminal_unlabeled) / (sum(c_terminal_unlabeled) + sum(c_terminal_labeled)))

        i += 1
        row_names = [f'{index + 1} ({sequence[index]})' for index in range(0, len(sequence))]
    
    # c-terminal value flipper
    for value in sheet[5]:
        if value != -1:
            sheet[6].append(1 - value)
        else:
            sheet[6].append(-1)

    # weighted average
    for index in range(len(sheet[2])):
        N = sheet[2][index]
        Nd = sheet[1][index]
        Cd = sheet[4][index]
        FlippedC = sheet[6][index] 

        # if c-terminal fraction flipped is not found
        if N != -1 and FlippedC == -1:
            sheet[7].append(N)
        elif N == -1 and FlippedC != -1:
            sheet[7].append(FlippedC)
        elif N == -1 and FlippedC == -1:
            sheet[7].append('---')
        else:
            sheet[7].append(
                N * (Nd / (Nd + Cd)) + (FlippedC * (Cd / (Nd + Cd)))
            )

    # denominator sums
    for index in range(len(sheet[1])):
        value = 0
        if sheet[1][index] != -1:
            value += sheet[1][index]
        if sheet[4][index] != -1:
            value += sheet[4][index]
        if value == 0:
            sheet[8].append('---')
        else:
            sheet[8].append(value)

    # relative label calculator
    for index in range(len(sheet[7])):
        if index == 0 or index == len(sheet[7]) - 1:
            sheet[9].append('---')
        elif index == 1:
            sheet[9].append(sheet[7][index])
        else:
            sheet[9].append(sheet[7][index] - sheet[7][index - 1])

    # rgb calculator
    values_to_normalize = [float(val) for val in sheet[8] if isinstance(val, (int, float))]
    for index in range(len(sheet[8])):
        if sheet[8][index] != '---':
            normalized_value = (sheet[8][index] - min(values_to_normalize)) * (255 / (max(values_to_normalize) - min(values_to_normalize)))
            sheet[10].append(f'color [255, {normalized_value}, {normalized_value}]')
        else:
            sheet[10].append('---')

    # color script generator
    for index in range(len(sheet[10])):
        color_script += f'select {index+1}\n{sheet[10][index]}\n'

    # replaces placeholder -1 with ---
    for index in range(len(sheet[2])):
        if sheet[2][index] == -1:
            sheet[2][index] = '---'

    for index in range(len(sheet[5])):
        if sheet[5][index] == -1:
            sheet[5][index] = '---'
    
    for index in range(len(sheet[5])):
        if sheet[6][index] == -1:
            sheet[6][index] = '---'

    return list(zip(*[
        [''] + row_names,
        ['n-terminal numerators'] + [value for value in sheet[0]],
        ['n-terminal denominators'] + [value for value in sheet[1]],
        ['n-terminal fraction'] + [value for value in sheet[2]],
        ['c-terminal numerators'] + [value for value in sheet[3]],
        ['c-terminal denominators'] + [value for value in sheet[4]],
        ['c-terminal fraction'] + [value for value in sheet[5]],
        ['c-terminal fraction flipped'] + [value for value in sheet[6]],
        ['weighted averages'] + [value for value in sheet[7]],
        ['denominator sum'] + [value for value in sheet[8]],
        ['relative label'] + [value for value in sheet[9]],
        ['rgb value'] + [value for value in sheet[10]]
    ]))

file_name = prompt('Please enter the name of the spectroscopy results file from Protein Prospector:')
peptide_sequence = prompt('Please enter the peptide sequence (default: *DAEFRHDSGYEVHHQK*): ') or '*DAEFRHDSGYEVHHQK*'

data = prepareData(file_name)
processed_data = processData(data, peptide_sequence)

generated_sheet = sheetGenerator(processed_data, peptide_sequence)

current_date = datetime.now().strftime('%m:%d:%Y')

with open(f'fragment-intensities-{current_date}.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(processed_data)

with open(f'processed-fragment-intensities-{current_date}.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(generated_sheet)

with open(f'color-script-{current_date}.txt', 'w', newline='') as file:
    file.write(color_script)