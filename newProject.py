#!/usr/bin/python3
import subprocess
from os import walk
# import sys

# sys.path.append('../src')
# import readHEFTConfig


def generateHEFTConfig(projName, n_bare, n_ch):
    """ Generate the main config file for HEFT """
    ch = []
    ch_pWaves = []
    partialWaves = ['S', 'P', 'D', 'F', 'G', 'H']
    dirName = f'{projName}_{n_bare}b{n_ch}c'

    # -------------------------Get channel data-------------------------
    for i in range(n_ch):
        particle1 = input(f'Enter the first particle for channel {i+1}: ')
        particle2 = input(f'Enter the second particle for channel {i+1}: ')
        while True:
            ch_pWave = input('Enter the partial wave which this channel undergoes: ')
            if ch_pWave.upper() in partialWaves:
                break
            else:
                print(f'Enter a valid partial wave ({partialWaves})')
                print('')
        ch.append([particle1, particle2])
        ch_pWaves.append(ch_pWave.upper())

    if (n_ch>1):
        while True:
            onCh = int(input('Which number channel is on-shell? '))
            if (onCh<=n_ch and onCh>0):
                break
            else:
                print(f'Enter a number between 1 and {n_ch}')
    else:
        onCh = 1

    # -----------------------Write to config file-----------------------
    with open(f'{dirName}/HEFT.config', 'w') as f:
        f.write(f'# {dirName} config\n\n')
        f.write(f'n_channels   {n_ch}\n')
        f.write(f'n_bare       {n_bare}\n\n')
        f.write('# Channel                   partialWave     latex\n')
        for i in range(n_ch):
            f.write(f'{ch[i][0]:<12}  {ch[i][1]:<12}  {ch_pWaves[i]}\n')
        f.write(f'\nOnshellChannel    {onCh}\n')
        f.write(f'useCustomMasses   F\n')
        f.write('particles.in\n\n')
        f.write('# Bare state labels\n')
        for i in range(n_bare):
            f.write(f'Bare{i+1}\n')
        f.write('\n# g (2-1 potential), v (2-2),   u (regulator)\n')
        f.write('A                    A          A')



def generateAllFitsFile(projName, n_bare, n_ch):
    """ Generates the file which will store all sets of fit params """
    dirName = f'{projName}_{n_bare}b{n_ch}c'
    # Default values for parameters
    m_bare_def =  1.6
    g_def      =  0.1
    Lam_def    =  0.8
    v_def      = -0.1
    Lam_v_def  =  0.8

    # Parameter end points
    nParamsTotal = int(n_bare + 2*n_ch*n_bare + (n_ch*(n_ch+1))/2 + n_ch)
    bare_end     = int(n_bare)
    g_end        = int(bare_end + n_ch*n_bare)
    Lam_end      = int(g_end + n_ch*n_bare)
    v_end        = int(Lam_end + (n_ch*(n_ch+1))/2)
    Lam_v_end    = int(v_end + n_ch)

    # Create the header, first amd zero string
    defStr   = '    n    '
    firstStr = '    1    '
    zeroStr  = '    0    '
    for nb in range(n_bare):
        defStr = defStr + f'm_b{nb+1}          '
        firstStr = firstStr + f'{m_bare_def:.8f}    '
        zeroStr  = zeroStr + f'{0:.8f}    '
    for nb in range(n_bare):
        for nc in range(n_ch):
            defStr = defStr + f'g_b{nb+1}c{nc+1}        '
            firstStr = firstStr + f'{g_def:.8f}    '
            zeroStr  = zeroStr + f'{0:.8f}    '
    for nb in range(n_bare):
        for nc in range(n_ch):
            defStr = defStr + f'Lam_b{nb+1}c{nc+1}      '
            firstStr = firstStr + f'{Lam_def:.8f}    '
            zeroStr  = zeroStr + f'{0:.8f}    '
    for nc1 in range(n_ch):
        for nc2 in range(nc1, n_ch):
            defStr = defStr + f'v_c{nc1+1}c{nc2+1}       '
            firstStr = firstStr + f'{v_def:.8f}   '
            zeroStr  = zeroStr + f' {0:.8f}   '
    for nc in range(n_ch):
        defStr = defStr + f'Lam_v_c{nc+1}      '
        firstStr = firstStr + f'{Lam_v_def:.8f}    '
        zeroStr  = zeroStr + f'{0:.8f}    '
        # firstStr = firstStr + '{:.8f}'.format(Lam_v_def) + '    '
    defStr = defStr + '          chi2     Notes\n'
    firstStr = firstStr + '          0.0      Default\n'
    zeroStr =  zeroStr  + '          0.0\n'


    # Create the allFits file based on the HEFT.config info
    with open(f'{dirName}/allFits.params', 'w') as f:
        f.write(f'    nParam     {nParamsTotal}\n')
        f.write('    iChoice    1\n')
        f.write(f'    paramEnds  {bare_end}  {g_end}  '
                + f'{Lam_end}  {v_end}  {Lam_v_end}  \n')
        f.write(defStr)
        f.write(firstStr)
        f.write(zeroStr)


def generateFittingConfig(projName, n_bare, n_ch):
    """ Generates the file which will store the fitting config """
    dirName = f'{projName}_{n_bare}b{n_ch}c'
    # Create the string of param names/isActive states
    paramStr   = '# '
    isActiveStr = '  '
    for nb in range(n_bare):
        paramStr = paramStr + f'm_b{nb+1}      '
        isActiveStr = isActiveStr + 'T         '
    for nb in range(n_bare):
        for nc in range(n_ch):
            paramStr = paramStr + f'g_b{nb+1}c{nc+1}    '
            isActiveStr = isActiveStr + 'T         '
    for nb in range(n_bare):
        for nc in range(n_ch):
            paramStr = paramStr + f'Lam_b{nb+1}c{nc+1}  '
            isActiveStr = isActiveStr + 'T         '
    for nc1 in range(n_ch):
        for nc2 in range(nc1, n_ch):
            paramStr = paramStr + f' v_c{nc1+1}c{nc2+1}    '
            isActiveStr = isActiveStr + ' T         '
    for nc in range(n_ch):
        paramStr = paramStr + f'Lam_v_c{nc+1}  '
        isActiveStr = isActiveStr + 'T         '
    paramStr = paramStr + '\n'
    isActiveStr = isActiveStr + '\n'

    with open(f'{dirName}/HEFTFitting.config', 'w') as f:
        f.write(f'#      Bounds                  Errors\n')
        for nb in range(n_bare):
            f.write(f'm{nb+1}     1.2         3.0         1e-5\n')
        f.write('g     -1.0         1.0	       1e-4\n')
        f.write('Lam    0.6         1.6	       1e-5\n')
        f.write('v     -1.0         1.0        1e-4\n')
        f.write('Lamv   0.6         1.6	       1e-5\n\n')

        f.write('# Which parameters are active\n')
        f.write(paramStr)
        f.write(isActiveStr)


def main():
    # ---------------------Setup new file directory---------------------
    projName = input('Enter the name of the new project: ')
    n_bare = int(input('Enter the number of bare states: '))
    n_ch = int(input('Enter the number of scattering channels: '))

    dirName = f'{projName}_{n_bare}b{n_ch}c'
    print(f'Creating directory {dirName}')
    proc = subprocess.run(f'mkdir {dirName}', shell=True)
    print('Copying these files:')
    proc = subprocess.run(f'ls Templates/*', shell=True)
    proc = subprocess.run(f'cp Templates/* {dirName}/', shell=True)
    proc = subprocess.run(f'ln src/makefile {dirName}/makefile', shell=True)
    proc = subprocess.run(f'mkdir {dirName}/figs', shell=True)
    proc = subprocess.run(f'mkdir {dirName}/data', shell=True)
    print()

    # ----------------------Generate required files---------------------
    print('Generating general HEFT Config for this system')
    generateHEFTConfig(projName, n_bare, n_ch)
    print()
    print('Generating allFits.params')
    generateAllFitsFile(projName, n_bare, n_ch)
    print('Generating the fitting config file')
    generateFittingConfig(projName, n_bare, n_ch)
    print('Done!')

    # ----------------------Select scattering data----------------------
    print('\nSelect the desired scattering data:')
    files = []
    for (_, _, filenames) in walk(f'ScatteringData/'):
        files.extend(filenames)
        nFiles = len(files)
    for i,filename in enumerate(files):
        print(f'{i+1}: {filename}')
    print('0: No scattering data')

    # Force a correct user input
    while True:
        choice = abs(int(input()))
        if (choice in range(0,nFiles+1)):
            break
        else:
            print('Enter a valid input')

    if (choice>0):
        res = subprocess.run(f'cp ScatteringData/{files[choice-1]} {dirName}/dataInf.in'
                              , shell=True, capture_output=True)
        stderr = res.stderr.decode('utf-8') # check copy was successful
        if (stderr==''):
            print(f'{files[choice-1]} copied')
            res = subprocess.run(f'cp {dirName}/dataInf.in {dirName}/dataInf_orig.in'
                                 , shell=True)
            res = subprocess.run(f'cp {dirName}/dataInf.in {dirName}/dataInf_pseudo.in'
                                 , shell=True)
        else:
            print(f'Error in copying:')
            print(f'   stderr: {stderr}')


if __name__ == '__main__':
    main()
