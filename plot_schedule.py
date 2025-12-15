import argparse
import matplotlib.pyplot as plt


def read_schedule(schedule_file):
    temps = []
    current_km = []
    best_km = []

    with open(schedule_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            temps.append(float(parts[0]))
            current_km.append(float(parts[1]))
            best_km.append(float(parts[2]))

    return temps, current_km, best_km


def main():
    parser = argparse.ArgumentParser(description="Plot annealing schedule")
    parser.add_argument("schedule_file")
    parser.add_argument("--out", default=None)
    parser.add_argument("--no_show", action="store_true")
    args = parser.parse_args()

    temps, current_km, best_km = read_schedule(args.schedule_file)
    if len(temps) == 0:
        raise RuntimeError("No schedule points found")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(temps, current_km, lw=1, label="current")
    ax.plot(temps, best_km, lw=2, label="best")
    ax.set_xscale("log")
    ax.set_xlabel("Temperature")
    ax.set_ylabel("Total distance (km)")
    ax.set_title("Annealing schedule")
    ax.legend()

    out_file = args.out
    if out_file is None:
        base = args.schedule_file.rsplit(".", 1)[0]
        out_file = base + ".png"
    fig.savefig(out_file, format="png", dpi=150, facecolor="white")

    if not args.no_show:
        print('close plot or "^C" to exit')
        try:
            plt.show()
        except KeyboardInterrupt:
            plt.close("all")


if __name__ == "__main__":
    main()
