import matplotlib.pyplot as plt
import sys
from abifpy import Trace
import numpy as np

if len(sys.argv) < 2:
    print("Usage: python inspect_ab1_abifpy.py input.ab1")
    sys.exit(1)

fn = sys.argv[1]
abi = Trace(fn)


def get_tag(tag):
    """Safe getter; some AB1 fields may not exist or be empty."""
    try:
        return abi.get_data(tag)
    except Exception:
        return None


# ---- PBAS1 ----
PBAS1 = get_tag("PBAS1")
PBAS1_len = len(PBAS1) if PBAS1 is not None else 0

print("PBAS1 seq len:", PBAS1_len)
print("PBAS1 (head):", PBAS1[:100] if PBAS1 else "<missing>")

# ---- Key tags ----
for k in ["PLOC1", "DATA9", "DATA10", "DATA11", "DATA12",  "PCON1"]:
    v = get_tag(k)
    if v is None:
        print(f"{k}: <MISSING>")
        continue

    try:
        ln = len(v)
        head = v[:100]
        tail = v[-100:]
    except Exception:
        ln = None
        head = tail = "<unprintable>"

    print(f"{k}: type={type(v).__name__}, len={ln}")
    print(f"  head={head}")
    print(f"  tail={tail}")

# ---- PLOC1 details ----
PLOC1 = get_tag("PLOC1")
if PLOC1 is not None:
    try:
        print("PLOC1[0..9]:", PLOC1[:100])
        print("PLOC1[-10:]:", PLOC1[-100:])
        print("PLOC1 MAX: ", max(PLOC1))
        print("PLOC1 RANGE", np.array(PLOC1[3110:3130]))
        print("PLOC1 RANGE", np.array(PLOC1[3119:3122]).astype(np.int16).view(np.uint16) )
    except Exception as e:
        print("Could not slice PLOC1:", e)

# ---- combined trace stats ----
channels = []
for k in ["DATA9", "DATA10", "DATA11", "DATA12"]:
    v = get_tag(k)
    if v is not None:
        channels.append(np.array(v, dtype=float))
    else:
        print(f"{k} not found. ")

if channels:
    s = np.sum(channels, axis=0)
    print("combined trace length:", len(s))
    print("combined trace sample:", s[:100])
    print("combined mean/std:", float(s.mean()), float(s.std()))


def load_ab1():
    seq = get_tag("PBAS1")
    ploc = np.array(get_tag("PLOC1"))
    G = np.array(get_tag("DATA9"), dtype=float)
    A = np.array(get_tag("DATA10"), dtype=float)
    T = np.array(get_tag("DATA11"), dtype=float)
    C = np.array(get_tag("DATA12"), dtype=float)

    return seq, ploc, A, C, G, T


def plot_ab1(start=None, end=None):
    seq, ploc, A, C, G, T = load_ab1()

    if start is None:
        start = 0
    if end is None or end > len(seq):
        end = len(seq)

    region_idx = np.arange(start, end)
    region_peaks = ploc[start:end]

    # ---- NEW: refine peaks using actual local maxima ----
    refined_peaks = []
    for p in region_peaks:
        win_left = max(0, p - 12)
        win_right = min(len(A), p + 12)

        # Sum trace to find strongest peak among A,C,G,T
        window_sum = (
            A[win_left:win_right] +
            C[win_left:win_right] +
            G[win_left:win_right] +
            T[win_left:win_right]
        )

        # true peak index inside trace
        local_max_idx = np.argmax(window_sum)
        refined_peaks.append(win_left + local_max_idx)

    refined_peaks = np.array(refined_peaks)

    # trace plotting range, based on refined peaks
    x_min = max(refined_peaks[0] - 50, 0)
    x_max = min(refined_peaks[-1] + 50, len(A))
    x = np.arange(x_min, x_max)

    plt.figure(figsize=(18, 6))

    plt.plot(x, A[x], label="A")
    plt.plot(x, C[x], label="C")
    plt.plot(x, G[x], label="G")
    plt.plot(x, T[x], label="T")

    # ---- mark refined peaks ----
    for peak_x, base in zip(refined_peaks, seq[start:end]):
        height = max(A[peak_x], C[peak_x], G[peak_x], T[peak_x])
        plt.text(peak_x, height + 20, base,
                 ha='center', va='bottom', fontsize=9)
        plt.axvline(peak_x, color='gray', linewidth=0.4, alpha=0.3)

    plt.xlim(x_min, x_max)
    plt.xlabel("Trace index")
    plt.ylabel("Signal intensity")
    plt.title(f"Basecalls {start}..{end}")
    plt.legend()
    plt.tight_layout()
    plt.savefig("show.png", dpi=150)
    plt.show()


plot_ab1(start=3120, end=3150)
