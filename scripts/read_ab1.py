from Bio import SeqIO
import numpy as np
from abifpy import Trace


def main():
    # fname = "/data1/ab1-debug/kangwei-ab1-data/Sequence/S22512022101-A1_1.ab1"
    fname = ""
    yummy = Trace(fname)
    print(np.array(yummy.get_data("PLOC1")).astype(np.int16).view(np.uint16))
    print(np.array(yummy.get_data("DATA9"))[:50])
    print(np.array(yummy.get_data("DATA10"))[:50])
    print(np.array(yummy.get_data("DATA11"))[:50])
    print(np.array(yummy.get_data("DATA12"))[:50])
    print(yummy)
    pass


if __name__ == "__main__":
    main()
    pass
