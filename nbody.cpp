#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

struct Simulation {
    int n = 0;

    // Per-particle arrays
    std::vector<double> m;
    std::vector<double> x, y, z;
    std::vector<double> vx, vy, vz;
    std::vector<double> fx, fy, fz;

    void resize(int count) {
        n = count;
        m.assign(n, 0.0);
        x.assign(n, 0.0); y.assign(n, 0.0); z.assign(n, 0.0);
        vx.assign(n, 0.0); vy.assign(n, 0.0); vz.assign(n, 0.0);
        fx.assign(n, 0.0); fy.assign(n, 0.0); fz.assign(n, 0.0);
    }

    void reset_forces() {
        std::fill(fx.begin(), fx.end(), 0.0);
        std::fill(fy.begin(), fy.end(), 0.0);
        std::fill(fz.begin(), fz.end(), 0.0);
    }

    // O(n^2) force calculation with softening, symmetric updates
    void compute_forces(double G, double softening) {
        reset_forces();

        const double eps2 = softening * softening;

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double dx = x[j] - x[i];
                double dy = y[j] - y[i];
                double dz = z[j] - z[i];

                double r2 = dx*dx + dy*dy + dz*dz + eps2;
                double inv_r = 1.0 / std::sqrt(r2);
                double inv_r3 = inv_r * inv_r * inv_r;

                double scale = G * m[i] * m[j] * inv_r3;

                double fxi = scale * dx;
                double fyi = scale * dy;
                double fzi = scale * dz;

                fx[i] += fxi;  fy[i] += fyi;  fz[i] += fzi;
                fx[j] -= fxi;  fy[j] -= fyi;  fz[j] -= fzi;
            }
        }
    }

    // Semi-implicit Euler like the handout
    void step(double dt) {
        for (int i = 0; i < n; i++) {
            double ax = fx[i] / m[i];
            double ay = fy[i] / m[i];
            double az = fz[i] / m[i];

            vx[i] += ax * dt;
            vy[i] += ay * dt;
            vz[i] += az * dt;

            x[i] += vx[i] * dt;
            y[i] += vy[i] * dt;
            z[i] += vz[i] * dt;
        }
    }

    // Output format: one state per line, tab-separated.
    // Each line: n, then for each particle:
    // mass x y z vx vy vz fx fy fz
    void dump(std::ostream& out) const {
        out << n;
        out << std::setprecision(15);
        for (int i = 0; i < n; i++) {
            out << '\t' << m[i]
                << '\t' << x[i] << '\t' << y[i] << '\t' << z[i]
                << '\t' << vx[i] << '\t' << vy[i] << '\t' << vz[i]
                << '\t' << fx[i] << '\t' << fy[i] << '\t' << fz[i];
        }
        out << '\n';
    }
};

static bool is_integer(const std::string& s) {
    if (s.empty()) return false;
    size_t i = 0;
    if (s[0] == '-' || s[0] == '+') i = 1;
    for (; i < s.size(); i++) {
        if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
    }
    return true;
}

static void init_random(Simulation& sim, int n, unsigned int seed = 42) {
    sim.resize(n);

    std::mt19937 rng(seed);

    // These ranges are "reasonable" for a simple demo
    std::uniform_real_distribution<double> mass_dist(1e22, 1e26);
    std::uniform_real_distribution<double> pos_dist(-1e11, 1e11);
    std::uniform_real_distribution<double> vel_dist(-2e4, 2e4);

    for (int i = 0; i < n; i++) {
        sim.m[i]  = mass_dist(rng);
        sim.x[i]  = pos_dist(rng);
        sim.y[i]  = pos_dist(rng);
        sim.z[i]  = pos_dist(rng);
        sim.vx[i] = vel_dist(rng);
        sim.vy[i] = vel_dist(rng);
        sim.vz[i] = vel_dist(rng);
    }
}

static void init_preset_sem(Simulation& sim) {
    // Very simple Sun-Earth-Moon-ish starting setup
    sim.resize(3);

    // Masses (kg)
    sim.m[0] = 1.9891e30; // Sun
    sim.m[1] = 5.972e24;  // Earth
    sim.m[2] = 7.342e22;  // Moon

    // Positions (m) on x-axis for simplicity
    sim.x[0] = 0.0;       sim.y[0] = 0.0;       sim.z[0] = 0.0;
    sim.x[1] = 1.496e11;  sim.y[1] = 0.0;       sim.z[1] = 0.0;   // 1 AU
    sim.x[2] = 1.496e11 + 3.844e8; sim.y[2] = 0.0; sim.z[2] = 0.0; // Earth + Moon distance

    // Velocities (m/s) roughly circular in +y
    sim.vx[0] = 0.0; sim.vy[0] = 0.0;   sim.vz[0] = 0.0;
    sim.vx[1] = 0.0; sim.vy[1] = 29780; sim.vz[1] = 0.0;
    sim.vx[2] = 0.0; sim.vy[2] = 29780 + 1022; sim.vz[2] = 0.0;

    sim.reset_forces();
}

static bool load_from_tsv_first_line(Simulation& sim, const std::string& path) {
    std::ifstream in(path);
    if (!in) return false;

    std::string line;
    if (!std::getline(in, line)) return false;

    std::istringstream ss(line);
    int n;
    ss >> n;
    if (!ss || n <= 0) return false;

    sim.resize(n);

    // The recommended format contains mass, pos(3), vel(3), force(3) per particle
    for (int i = 0; i < n; i++) {
        double mass, px, py, pz, vxx, vyy, vzz, fxx, fyy, fzz;
        ss >> mass >> px >> py >> pz >> vxx >> vyy >> vzz;

        if (!ss) return false;

        // forces might not exist in some files
        if (!(ss >> fxx >> fyy >> fzz)) {
            fxx = fyy = fzz = 0.0;
            ss.clear();
        }

        sim.m[i]  = mass;
        sim.x[i]  = px;  sim.y[i]  = py;  sim.z[i]  = pz;
        sim.vx[i] = vxx; sim.vy[i] = vyy; sim.vz[i] = vzz;
        sim.fx[i] = fxx; sim.fy[i] = fyy; sim.fz[i] = fzz;
    }

    return true;
}

static void print_usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " <n_or_file_or_preset> <dt> <steps> <dumpEvery> [output.tsv]\n\n"
        << "Where:\n"
        << "  <n_or_file_or_preset> is either:\n"
        << "    - an integer N: random simulation with N particles\n"
        << "    - a file path: load initial state from first line of a TSV file\n"
        << "    - the string 'sem': built-in Sun/Earth/Moon preset\n"
        << "  <dt> is the time step (double)\n"
        << "  <steps> is number of iterations (int)\n"
        << "  <dumpEvery> is how often to output (int). 1 = every step.\n"
        << "  [output.tsv] optional output file (default: stdout)\n\n"
        << "Example:\n"
        << "  " << prog << " 100 1 10000 10 out.tsv\n"
        << "  " << prog << " solar.tsv 200 5000 10 solar_out.tsv\n";
}

int main(int argc, char** argv) {
    if (argc < 5) {
        print_usage(argv[0]);
        return 1;
    }

    std::string mode = argv[1];
    double dt = std::atof(argv[2]);
    int steps = std::atoi(argv[3]);
    int dumpEvery = std::atoi(argv[4]);

    if (dt <= 0 || steps < 0 || dumpEvery <= 0) {
        std::cerr << "Error: dt must be > 0, steps >= 0, dumpEvery > 0\n";
        return 1;
    }

    std::ofstream fout;
    std::ostream* out = &std::cout;
    if (argc >= 6) {
        fout.open(argv[5]);
        if (!fout) {
            std::cerr << "Error: could not open output file: " << argv[5] << "\n";
            return 1;
        }
        out = &fout;
    }

    Simulation sim;

    if (mode == "sem") {
        init_preset_sem(sim);
    } else if (is_integer(mode)) {
        int n = std::stoi(mode);
        if (n <= 0) {
            std::cerr << "Error: N must be > 0\n";
            return 1;
        }
        init_random(sim, n, 42);
    } else {
        if (!load_from_tsv_first_line(sim, mode)) {
            std::cerr << "Error: could not load initial state from file: " << mode << "\n";
            return 1;
        }
    }

    // Constants
    const double G = 6.674e-11;     // gravitational constant
    const double softening = 1e3;   // meters (prevents divide-by-zero)

    // Output initial state (step 0)
    sim.compute_forces(G, softening);
    sim.dump(*out);

    for (int t = 1; t <= steps; t++) {
        sim.compute_forces(G, softening);
        sim.step(dt);

        if (t % dumpEvery == 0) {
            sim.dump(*out);
        }
    }

    return 0;
}
