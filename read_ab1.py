from Bio import SeqIO

from abifpy import Trace

def main():
    fname = "/data-slow/kangwei-deliver/kangwei-deliver/S22509070002-Epi5A-1.ab1"
    yummy = Trace(fname)
    print(yummy)
    pass


if __name__ == "__main__":
    main()
    pass
