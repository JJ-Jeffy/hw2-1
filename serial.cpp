#include "common.h"
#include <cmath>
#include <vector>

// Global declaration of bins
std::vector<std::vector<int>> bins;
int binCountX, binCountY; // Number of bins in each dimension
int binSize; 

// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

// Integrate the ODE
void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}


void init_simulation(particle_t* parts, int num_parts, double size) {
    binSize = std::max(static_cast<int>(std::ceil(cutoff * 0.5)), 1); // problem 
    binCountX = std::ceil(size / binSize);
    binCountY = std::ceil(size / binSize);
    bins.resize(binCountX * binCountY);
}

// Assign particles to bins 
void assign_particles_to_bins(particle_t* parts, int num_parts, double size) {
    // Clear bins
    for (auto& bin : bins) {
        bin.clear();
    }

    // Assign particles to bins
    for (int i = 0; i < num_parts; ++i) {
        int binX = parts[i].x / binSize; 
        int binY = parts[i].y / binSize;
        int binIndex = binY * binCountX + binX; // Calculate linear index
        bins[binIndex].push_back(i);
    }
}

std::vector<std::pair<int, int>> get_neighboring_bins(int binX, int binY, int maxX, int maxY) {
    std::vector<std::pair<int,int>> neighbors;
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            int newX = binX + dx, newY = binY + dy;
            if (newX >= 0 && newX < maxX && newY >= 0 && newY < maxY) {
                neighbors.emplace_back(newX, newY); // Store as pairs for compatibility
            }
        }
    }
    return neighbors;
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Re-assign particles to bins to account for movement
    assign_particles_to_bins(parts, num_parts, size);

    // Compute forces with binning
    for (int i = 0; i < num_parts; ++i) {
        parts[i].ax = parts[i].ay = 0;
        int binX = parts[i].x / binSize;
        int binY = parts[i].y / binSize;
        int binIndex = binY * binCountX + binX;

        auto neighbors = get_neighboring_bins(binX, binY, binCountX, binCountY);
        for (auto &[nbX, nbY] : neighbors) {
            int nbIndex = nbY * binCountX + nbX;
            for (int j : bins[nbIndex]) {
                apply_force(parts[i], parts[j]);
            }
        }
    }

    // Move particles
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }
}
