#!/usr/bin/python3


class HEFTConfig:
    def __init__(self):
        self._n_ch = 0
        self._n_bare = 0
        self._chs = []
        self._pWaves = []
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

    def readHEFTConfigFile(self):
        with open('HEFT.config', 'r') as f:
            f.readline()
            f.readline()
            line = f.readline()
            self.n_ch = int(line[-3:])
            line = f.readline()
            self.n_bare = int(line[-3:])

            f.readline()
            f.readline()
            self.chs = []
            for i_ch in range(self.n_ch):
                line = f.readline()
                self.chs.append(line[:-2])
                self.pWaves.append(line[-2:])

            f.readline()
            f.readline()
            f.readline()
            f.readline()
            f.readline()

            for i_bare in range(self.n_bare):
                line = f.readline()
                self.bareStates.append(line.strip())

            f.readline()
            f.readline()

            line = f.readline()
            potChoice = line.replace(" ", "")

    def printHEFTInfo(self):
        nbmc = f'{self.n_bare}b{self.n_ch}c'
        print(f'--------------------HEFT {nbmc}--------------------')
        print('Bare States:')
        for bare in self.bareStates:
            print(f'\t{bare}')
        print('Scattering Channels:')
        for i in range(len(self.chs)):
            print(f'\t{self.chs[i]} {self.pWaves[i].strip()}-wave')
        print(f'-------------------------------------------------')
        print()

def main():
    HEFT = HEFTConfig()
    HEFT.readHEFTConfigFile()
    HEFT.printHEFTInfo()


if __name__ == '__main__':
    main()
