#include <algorithm>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#define EARTH_RADIUS_KM 6371.0
#define PI 3.141592653589793238462643383279502884

typedef struct {
    double lon_deg;
    double lat_deg;
    double lon_rad;
    double lat_rad;
} city_t;

typedef struct {
    std::size_t n;
    std::vector<double> d;
} dist_matrix_t;

typedef struct {
    std::string input_file;
    std::string init_route_file;
    std::string opt_route_file;
    std::string schedule_file;

    bool shuffle_init;
    std::uint64_t seed;

    std::string move;
    int max_seg;

    double tmax;
    double tmin;
    double alpha;
    int melt_sweeps;
    int sweeps_per_temp;
    int tmax_trials;
} args_t;

static double deg_to_rad(double deg) {
    return deg * PI / 180.0;
}

static bool is_comment_or_blank(const std::string &line) {
    bool saw_non_space = false;
    for (char ch : line) {
        if (!std::isspace(static_cast<unsigned char>(ch))) {
            saw_non_space = true;
            return ch == '#';
        }
    }
    return !saw_non_space;
}

static std::vector<city_t> read_cities(const std::string &filename) {
    std::ifstream fin(filename);
    if (!fin) {
        throw std::runtime_error("failed to open input file: " + filename);
    }

    std::vector<city_t> cities;
    std::string line;
    while (std::getline(fin, line)) {
        if (is_comment_or_blank(line)) {
            continue;
        }

        std::istringstream iss(line);
        double lon_deg = 0.0;
        double lat_deg = 0.0;
        if (!(iss >> lon_deg >> lat_deg)) {
            continue;
        }

        city_t c;
        c.lon_deg = lon_deg;
        c.lat_deg = lat_deg;
        c.lon_rad = deg_to_rad(lon_deg);
        c.lat_rad = deg_to_rad(lat_deg);
        cities.push_back(c);
    }

    if (cities.size() < 2) {
        throw std::runtime_error("need at least 2 cities");
    }

    return cities;
}

static dist_matrix_t dist_matrix_create(std::size_t n) {
    dist_matrix_t m;
    m.n = n;
    m.d.assign(n * n, 0.0);
    return m;
}

static double dist_matrix_get(const dist_matrix_t &m, std::size_t i, std::size_t j) {
    return m.d[i * m.n + j];
}

static void dist_matrix_set(dist_matrix_t &m, std::size_t i, std::size_t j, double v) {
    m.d[i * m.n + j] = v;
}

static double haversine_km(double lon1_rad, double lat1_rad, double lon2_rad, double lat2_rad) {
    double dlat = lat2_rad - lat1_rad;
    double dlon = lon2_rad - lon1_rad;

    double sin_dlat = std::sin(dlat / 2.0);
    double sin_dlon = std::sin(dlon / 2.0);
    double a = (sin_dlat * sin_dlat) + std::cos(lat1_rad) * std::cos(lat2_rad) * (sin_dlon * sin_dlon);
    if (a < 0.0) {
        a = 0.0;
    }
    if (a > 1.0) {
        a = 1.0;
    }
    double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
    return EARTH_RADIUS_KM * c;
}

static dist_matrix_t build_dist_matrix(const std::vector<city_t> &cities) {
    std::size_t n = cities.size();
    dist_matrix_t dist = dist_matrix_create(n);

    for (std::size_t i = 0; i < n; i += 1) {
        dist_matrix_set(dist, i, i, 0.0);
        for (std::size_t j = i + 1; j < n; j += 1) {
            double dij = haversine_km(cities[i].lon_rad, cities[i].lat_rad, cities[j].lon_rad, cities[j].lat_rad);
            dist_matrix_set(dist, i, j, dij);
            dist_matrix_set(dist, j, i, dij);
        }
    }

    return dist;
}

static double tour_length_km(const dist_matrix_t &dist, const std::vector<int> &route) {
    std::size_t n = route.size();
    double total = 0.0;
    for (std::size_t i = 0; i < n; i += 1) {
        std::size_t j = (i + 1) % n;
        total += dist_matrix_get(
            dist,
            static_cast<std::size_t>(route[i]),
            static_cast<std::size_t>(route[j])
        );
    }
    return total;
}

static void write_route_lonlat(const std::string &filename, const std::vector<city_t> &cities, const std::vector<int> &route) {
    std::ofstream fout(filename);
    if (!fout) {
        throw std::runtime_error("failed to open output file: " + filename);
    }

    fout << std::setprecision(10);
    for (int idx : route) {
        const city_t &c = cities[static_cast<std::size_t>(idx)];
        fout << c.lon_deg << " " << c.lat_deg << "\n";
    }
}

static double delta_two_opt(const dist_matrix_t &dist, const std::vector<int> &route, int i, int j) {
    int n = static_cast<int>(route.size());
    int ip1 = (i + 1) % n;
    int jp1 = (j + 1) % n;

    int a = route[static_cast<std::size_t>(i)];
    int b = route[static_cast<std::size_t>(ip1)];
    int c = route[static_cast<std::size_t>(j)];
    int d = route[static_cast<std::size_t>(jp1)];

    double old_len = dist_matrix_get(dist, static_cast<std::size_t>(a), static_cast<std::size_t>(b))
        + dist_matrix_get(dist, static_cast<std::size_t>(c), static_cast<std::size_t>(d));
    double new_len = dist_matrix_get(dist, static_cast<std::size_t>(a), static_cast<std::size_t>(c))
        + dist_matrix_get(dist, static_cast<std::size_t>(b), static_cast<std::size_t>(d));

    return new_len - old_len;
}

static void apply_two_opt(std::vector<int> &route, int i, int j) {
    std::reverse(route.begin() + (i + 1), route.begin() + (j + 1));
}

static double delta_swap(const dist_matrix_t &dist, const std::vector<int> &route, int i, int j) {
    int n = static_cast<int>(route.size());
    if (i == j) {
        return 0.0;
    }
    if (i > j) {
        std::swap(i, j);
    }

    int im1 = (i - 1 + n) % n;
    int ip1 = (i + 1) % n;
    int jm1 = (j - 1 + n) % n;
    int jp1 = (j + 1) % n;

    int a = route[static_cast<std::size_t>(im1)];
    int b = route[static_cast<std::size_t>(i)];
    int c = route[static_cast<std::size_t>(ip1)];

    int d = route[static_cast<std::size_t>(jm1)];
    int e = route[static_cast<std::size_t>(j)];
    int f = route[static_cast<std::size_t>(jp1)];

    if (j == i + 1) {
        double old_len = dist_matrix_get(dist, static_cast<std::size_t>(a), static_cast<std::size_t>(b))
            + dist_matrix_get(dist, static_cast<std::size_t>(b), static_cast<std::size_t>(e))
            + dist_matrix_get(dist, static_cast<std::size_t>(e), static_cast<std::size_t>(f));
        double new_len = dist_matrix_get(dist, static_cast<std::size_t>(a), static_cast<std::size_t>(e))
            + dist_matrix_get(dist, static_cast<std::size_t>(e), static_cast<std::size_t>(b))
            + dist_matrix_get(dist, static_cast<std::size_t>(b), static_cast<std::size_t>(f));
        return new_len - old_len;
    }

    double old_len = dist_matrix_get(dist, static_cast<std::size_t>(a), static_cast<std::size_t>(b))
        + dist_matrix_get(dist, static_cast<std::size_t>(b), static_cast<std::size_t>(c))
        + dist_matrix_get(dist, static_cast<std::size_t>(d), static_cast<std::size_t>(e))
        + dist_matrix_get(dist, static_cast<std::size_t>(e), static_cast<std::size_t>(f));

    double new_len = dist_matrix_get(dist, static_cast<std::size_t>(a), static_cast<std::size_t>(e))
        + dist_matrix_get(dist, static_cast<std::size_t>(e), static_cast<std::size_t>(c))
        + dist_matrix_get(dist, static_cast<std::size_t>(d), static_cast<std::size_t>(b))
        + dist_matrix_get(dist, static_cast<std::size_t>(b), static_cast<std::size_t>(f));

    return new_len - old_len;
}

static void apply_swap(std::vector<int> &route, int i, int j) {
    std::swap(route[static_cast<std::size_t>(i)], route[static_cast<std::size_t>(j)]);
}

static void usage(const char *prog) {
    std::cerr << "usage: " << prog << " -i cities.dat [options]\n";
    std::cerr << "options:\n";
    std::cerr << "  -i <file>                  input cities file (required)\n";
    std::cerr << "  --init <file>              write initial route (lon lat) to file\n";
    std::cerr << "  -o <file>                  write optimized route (lon lat) to file\n";
    std::cerr << "  --schedule <file>          write annealing schedule (T current_km best_km)\n";
    std::cerr << "  --seed <u64>               RNG seed (default 1)\n";
    std::cerr << "  --no_shuffle               do not shuffle initial route\n";
    std::cerr << "  --move <two_opt|swap>      move type (default two_opt)\n";
    std::cerr << "  --max_seg <int>            max 2-opt segment length (default 200; <=0 means full)\n";
    std::cerr << "  --tmax <float>             starting temperature (0 => auto)\n";
    std::cerr << "  --tmin <float>             stopping temperature (default 1e-3)\n";
    std::cerr << "  --alpha <float>            cooling factor (default 0.95)\n";
    std::cerr << "  --melt <int>               melt sweeps (proposals = sweeps * N, default 20)\n";
    std::cerr << "  --sweeps <int>             sweeps per temperature (default 50)\n";
}

static args_t default_args(void) {
    args_t a;
    a.input_file = "";
    a.init_route_file = "";
    a.opt_route_file = "";
    a.schedule_file = "";

    a.shuffle_init = true;
    a.seed = 1;

    a.move = "two_opt";
    a.max_seg = 200;

    a.tmax = 0.0;
    a.tmin = 1e-3;
    a.alpha = 0.95;
    a.melt_sweeps = 20;
    a.sweeps_per_temp = 50;
    a.tmax_trials = 2000;

    return a;
}

static args_t parse_args(int argc, char **argv) {
    args_t args = default_args();

    for (int i = 1; i < argc; i += 1) {
        std::string a = argv[i];

        if (a == "-i" && i + 1 < argc) {
            args.input_file = argv[i + 1];
            i += 1;
        } else if (a == "--init" && i + 1 < argc) {
            args.init_route_file = argv[i + 1];
            i += 1;
        } else if (a == "-o" && i + 1 < argc) {
            args.opt_route_file = argv[i + 1];
            i += 1;
        } else if (a == "--schedule" && i + 1 < argc) {
            args.schedule_file = argv[i + 1];
            i += 1;
        } else if (a == "--seed" && i + 1 < argc) {
            args.seed = static_cast<std::uint64_t>(std::stoull(argv[i + 1]));
            i += 1;
        } else if (a == "--no_shuffle") {
            args.shuffle_init = false;
        } else if (a == "--move" && i + 1 < argc) {
            args.move = argv[i + 1];
            i += 1;
        } else if (a == "--max_seg" && i + 1 < argc) {
            args.max_seg = std::stoi(argv[i + 1]);
            i += 1;
        } else if (a == "--tmax" && i + 1 < argc) {
            args.tmax = std::stod(argv[i + 1]);
            i += 1;
        } else if (a == "--tmin" && i + 1 < argc) {
            args.tmin = std::stod(argv[i + 1]);
            i += 1;
        } else if (a == "--alpha" && i + 1 < argc) {
            args.alpha = std::stod(argv[i + 1]);
            i += 1;
        } else if (a == "--melt" && i + 1 < argc) {
            args.melt_sweeps = std::stoi(argv[i + 1]);
            i += 1;
        } else if (a == "--sweeps" && i + 1 < argc) {
            args.sweeps_per_temp = std::stoi(argv[i + 1]);
            i += 1;
        } else {
            usage(argv[0]);
            throw std::runtime_error("unknown or incomplete argument: " + a);
        }
    }

    if (args.input_file.empty()) {
        usage(argv[0]);
        throw std::runtime_error("missing required -i <file>");
    }

    auto infer_tag = [&](const std::string &infile) -> std::string {
        std::string base = infile;
        std::size_t slash = base.find_last_of("/\\");
        if (slash != std::string::npos) {
            base = base.substr(slash + 1);
        }
        std::size_t dot = base.find_last_of('.');
        if (dot != std::string::npos) {
            base = base.substr(0, dot);
        }
        return base;
    };

    std::string tag = infer_tag(args.input_file);
    if (args.init_route_file.empty()) {
        args.init_route_file = tag + "_init.dat";
    }
    if (args.opt_route_file.empty()) {
        args.opt_route_file = tag + "_opt.dat";
    }
    if (args.schedule_file.empty()) {
        if (tag.rfind("cities", 0) == 0) {
            args.schedule_file = "an" + tag.substr(6) + ".dat";
        } else {
            args.schedule_file = "an" + tag + ".dat";
        }
    }

    if (args.alpha <= 0.0 || args.alpha >= 1.0) {
        throw std::runtime_error("alpha must be in (0,1)");
    }
    if (args.tmin <= 0.0) {
        throw std::runtime_error("tmin must be > 0");
    }
    if (args.melt_sweeps < 0 || args.sweeps_per_temp <= 0) {
        throw std::runtime_error("melt must be >=0 and sweeps must be >0");
    }
    if (!(args.move == "two_opt" || args.move == "swap")) {
        throw std::runtime_error("move must be two_opt or swap");
    }

    return args;
}

static bool propose_two_opt_indices(int n, int max_seg, std::mt19937_64 &rng, int *i_out, int *j_out) {
    if (n < 4) {
        return false;
    }

    int max_len = n - 1;
    if (max_seg > 0 && max_seg < max_len) {
        max_len = max_seg;
    }
    if (max_len < 2) {
        return false;
    }

    std::uniform_int_distribution<int> uni_i(0, n - 3);
    int i = uni_i(rng);

    int max_j = i + max_len;
    if (max_j > n - 1) {
        max_j = n - 1;
    }
    if (max_j < i + 2) {
        return false;
    }

    std::uniform_int_distribution<int> uni_j(i + 2, max_j);
    int j = uni_j(rng);

    *i_out = i;
    *j_out = j;
    return true;
}

static double estimate_tmax(const dist_matrix_t &dist, const std::vector<int> &route, std::mt19937_64 &rng, const args_t &args) {
    int n = static_cast<int>(route.size());
    if (n < 4) {
        return 1.0;
    }

    std::uniform_int_distribution<int> uni_idx(0, n - 1);
    double max_delta = 0.0;

    for (int t = 0; t < args.tmax_trials; t += 1) {
        int i = 0;
        int j = 0;

        if (args.move == "two_opt") {
            if (!propose_two_opt_indices(n, args.max_seg, rng, &i, &j)) {
                continue;
            }
            double delta = delta_two_opt(dist, route, i, j);
            if (delta > max_delta) {
                max_delta = delta;
            }
        } else {
            i = uni_idx(rng);
            j = uni_idx(rng);
            if (i == j) {
                continue;
            }
            double delta = delta_swap(dist, route, i, j);
            if (delta > max_delta) {
                max_delta = delta;
            }
        }
    }

    if (max_delta <= 0.0) {
        max_delta = 1.0;
    }

    return 1.1 * max_delta;
}

int main(int argc, char **argv) {
    try {
        args_t args = parse_args(argc, argv);

        std::vector<city_t> cities = read_cities(args.input_file);
        dist_matrix_t dist = build_dist_matrix(cities);

        std::size_t n_sz = cities.size();
        int n = static_cast<int>(n_sz);

        std::vector<int> route(n_sz);
        for (std::size_t i = 0; i < n_sz; i += 1) {
            route[i] = static_cast<int>(i);
        }

        std::mt19937_64 rng(args.seed);
        if (args.shuffle_init) {
            std::shuffle(route.begin(), route.end(), rng);
        }

        double initial_km = tour_length_km(dist, route);
        double current_km = initial_km;

        std::vector<int> best_route = route;
        double best_km = initial_km;

        if (!args.init_route_file.empty()) {
            write_route_lonlat(args.init_route_file, cities, route);
        }

        double tmax = args.tmax;
        if (tmax <= 0.0) {
            tmax = estimate_tmax(dist, route, rng, args);
        }

        std::ofstream schedule(args.schedule_file);
        if (!schedule) {
            throw std::runtime_error("failed to open schedule file: " + args.schedule_file);
        }
        schedule << "# T current_km best_km\n";
        schedule << std::setprecision(10);

        std::uniform_real_distribution<double> uni01(0.0, 1.0);
        std::uniform_int_distribution<int> uni_idx(0, n - 1);

        auto t_start = std::chrono::high_resolution_clock::now();

        int melt_proposals = args.melt_sweeps * n;
        for (int step = 0; step < melt_proposals; step += 1) {
            double delta = 0.0;
            bool have_move = false;
            int i = 0;
            int j = 0;

            if (args.move == "two_opt") {
                have_move = propose_two_opt_indices(n, args.max_seg, rng, &i, &j);
                if (have_move) {
                    delta = delta_two_opt(dist, route, i, j);
                }
            } else {
                i = uni_idx(rng);
                j = uni_idx(rng);
                if (i != j) {
                    have_move = true;
                    delta = delta_swap(dist, route, i, j);
                }
            }

            if (!have_move) {
                continue;
            }

            double accept_p = 1.0;
            if (delta > 0.0) {
                accept_p = std::exp(-delta / tmax);
            }

            if (uni01(rng) < accept_p) {
                if (args.move == "two_opt") {
                    apply_two_opt(route, i, j);
                } else {
                    apply_swap(route, i, j);
                }
                current_km += delta;

                if (current_km < best_km) {
                    best_km = current_km;
                    best_route = route;
                }
            }
        }

        double T = tmax;
        while (T > args.tmin) {
            int proposals = args.sweeps_per_temp * n;

            for (int step = 0; step < proposals; step += 1) {
                double delta = 0.0;
                bool have_move = false;
                int i = 0;
                int j = 0;

                if (args.move == "two_opt") {
                    have_move = propose_two_opt_indices(n, args.max_seg, rng, &i, &j);
                    if (have_move) {
                        delta = delta_two_opt(dist, route, i, j);
                    }
                } else {
                    i = uni_idx(rng);
                    j = uni_idx(rng);
                    if (i != j) {
                        have_move = true;
                        delta = delta_swap(dist, route, i, j);
                    }
                }

                if (!have_move) {
                    continue;
                }

                double accept_p = 1.0;
                if (delta > 0.0) {
                    accept_p = std::exp(-delta / T);
                }

                if (uni01(rng) < accept_p) {
                    if (args.move == "two_opt") {
                        apply_two_opt(route, i, j);
                    } else {
                        apply_swap(route, i, j);
                    }
                    current_km += delta;

                    if (current_km < best_km) {
                        best_km = current_km;
                        best_route = route;
                    }
                }
            }

            schedule << T << " " << current_km << " " << best_km << "\n";
            T *= args.alpha;
        }

        auto t_end = std::chrono::high_resolution_clock::now();
        double elapsed_s = std::chrono::duration<double>(t_end - t_start).count();

        write_route_lonlat(args.opt_route_file, cities, best_route);

        std::cout << std::setprecision(10);
        std::cout << "n_cities " << n_sz << "\n";
        std::cout << "initial_km " << initial_km << "\n";
        std::cout << "best_km " << best_km << "\n";
        std::cout << "time_s " << elapsed_s << "\n";
        std::cout << "init_route_file " << args.init_route_file << "\n";
        std::cout << "opt_route_file " << args.opt_route_file << "\n";
        std::cout << "schedule_file " << args.schedule_file << "\n";

        return 0;
    } catch (const std::exception &e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
