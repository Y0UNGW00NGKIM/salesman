import sys

# simple coordinate "struct"
class Coord:
    def __init__(self, lon, lat):
        self.lon = lon
        self.lat = lat

def get_data(filename):
    cities = []
    with open(filename, "r") as fp:
        for line in fp:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # read first two numbers only (matches C++ "sscanf")
            parts = line.split()
            try:
                lon = float(parts[0])
                lat = float(parts[1])
            except (ValueError, IndexError):
                continue  # skip malformed lines

            cities.append(Coord(lon, lat))
    return cities

def main():
    if len(sys.argv) < 2:
        print("Please provide a data file path as argument")
        return

    filename = sys.argv[1]
    cities = get_data(filename)

    print(f"Read {len(cities)} cities from data file")
    print("Longitude  Latitude")

    for c in cities:
        print(f"{c.lon:.6f} {c.lat:.6f}")

if __name__ == "__main__":
    main()

