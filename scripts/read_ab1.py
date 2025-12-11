from Bio import SeqIO
import numpy as np
from abifpy import Trace


def main():
    fname = "/data1/ab1-debug/S22509070002-Epi5A-1.ab1"
    yummy = Trace(fname)
    print(np.array(yummy.get_data("PLOC1")).astype(np.int16).view(np.uint16))
    print(np.array(yummy.get_data("DATA9"))[:10])
    print(np.array(yummy.get_data("DATA10"))[:10])
    print(np.array(yummy.get_data("DATA11"))[:10])
    print(np.array(yummy.get_data("DATA12"))[:10])
    print(yummy)
    pass


if __name__ == "__main__":
    main()
    pass
