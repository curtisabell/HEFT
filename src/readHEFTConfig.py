#!/usr/bin/python3


# -------------Class to store details about each channel------------
class channel:
    def __init__(self):
        self._particle = ['', '']
        self._partialWave = ''
        self._isOnshell = False
        self._latexLabel = ''
        self._mass = [0.0, 0.0]

    @property
    def particle(self):
        return self._particle
    @particle.setter
    def particle(self, value):
        self._particle = value

    @property
    def partialWave(self):
        return self._partialWave
    @partialWave.setter
    def partialWave(self, value):
        self._partialWave = value

    @property
    def isOnshell(self):
        return self._isOnshell
    @isOnshell.setter
    def isOnshell(self, value):
        self._isOnshell = value

    @property
    def latexLabel(self):
        return self._latexLabel
    @latexLabel.setter
    def latexLabel(self, value):
        self._latexLabel = value

    @property
    def mass(self):
        return self._mass
    @mass.setter
    def mass(self, value):
        self._mass = value

    # makes the print output compatible with older scripts
    def __str__(self):
        return f'{self.particle[0]} {self.particle[1]}'



# -----Class to contain all info stored in the HEFT.config file-----
class HEFTConfig:
    def __init__(self):
        self._n_ch = 0
        self._n_bare = 0
        self._chs = []
        self._pWaves = []
        self._useCustomMasses = False
        self._customMassFile = ''
        self._bareStates = []

    @property
    def n_ch(self):
        return self._n_ch
    @n_ch.setter
    def n_ch(self, value):
        self._n_ch = value

    @property
    def chs(self):
        return self._chs
    @chs.setter
    def chs(self, value):
        self._chs = value

    @property
    def pWaves(self):
        return self._pWaves
    @pWaves.setter
    def pWaves(self, value):
        self._pWaves = value

    @property
    def useCustomMasses(self):
        return self._useCustomMasses
    @useCustomMasses.setter
    def useCustomMasses(self, value):
        self._useCustomMasses = value

    @property
    def customMassFile(self):
        return self._customMassFile
    @customMassFile.setter
    def customMassFile(self, value):
        self._customMassFile = value

    @property
    def bareStates(self):
        return self._bareStates
    @bareStates.setter
    def bareStates(self, value):
        self._bareStates = value

    @property
    def n_bare(self):
        return self._n_bare
    @n_bare.setter
    def n_bare(self, value):
        self._n_bare = value

    @property
    def projectName(self):
        return self._projectName
    @projectName.setter
    def projectName(self, value):
        self._projectName = value

    def readHEFTConfigFile(self):
        with open('HEFT.config', 'r') as f:
            # Read n_b and n_c info for this system
            self.projectName = f.readline().split()[1]
            f.readline()
            line = f.readline()
            self.n_ch = int(line[-3:])
            line = f.readline()
            self.n_bare = int(line[-3:])

            # Read info about the channel(s)
            f.readline()
            f.readline()
            self.chs = []
            for i_ch in range(self.n_ch):
                chLine = f.readline().split()
                self.chs.append(channel())
                self.chs[i_ch].particle[0] = chLine[0]
                self.chs[i_ch].particle[1] = chLine[1]
                self.chs[i_ch].partialWave = chLine[2]
                # self.chs[i_ch].latexLabel = chLine[3]

            f.readline()

            # get which channel is onshell
            line = f.readline()
            onshellChannel = int(line.split()[1])
            self.chs[onshellChannel-1].isOnshell = True

            line = f.readline()
            self.useCustomMasses = True if line.split()[1]=='T' else False
            line = f.readline()
            self.customMassFile = line.split()[0]

            f.readline()
            f.readline()

            # Read info about the bare state(s)
            for i_bare in range(self.n_bare):
                line = f.readline()
                self.bareStates.append(line.strip())

            f.readline()
            f.readline()

            line = f.readline()
            potChoice = line.replace(" ", "")

    # def establishParticleMasses(self):
    #     if self.useCustomMasses:
    #         particleMassFile = self.customMassFile
    #     else:
    #         particleMassFile = '../src/particles.in'

    #     with open(particleMassFile, 'r') as f:
    #         # TODO read this file, find the particles in each channel and set themasses

    def printHEFTInfo(self):
        nbmc = f'{self.n_bare}b{self.n_ch}c'
        print(f'--------------------HEFT {nbmc}--------------------')
        print('Bare States:')
        for bare in self.bareStates:
            print(f'\t{bare}')
        print('Scattering Channels:')
        for i in range(self.n_ch):
            print(f'\t{self.chs[i]} ({self.chs[i].partialWave}-wave)')
        print(f'-------------------------------------------------')
        print()

def main():
    HEFT = HEFTConfig()
    HEFT.readHEFTConfigFile()
    HEFT.printHEFTInfo()


if __name__ == '__main__':
    main()
