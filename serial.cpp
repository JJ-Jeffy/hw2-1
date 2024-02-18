#include "common.h"
#include <cmath>
#include <vector>

// Global declaration of bins
std::vector<std::vector<std::vector<int>>> bins;
int binSize; // This will be initialized based on cutoff

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
	// You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any particle simulation here
    // Initialize binSize based on global cutoff value
    binSize = static_cast<int>(cutoff * 2); 

    // Calculate the number of bins based on the simulation size and binSize
    int numBinsX = ceil(size / binSize);
    int numBinsY = ceil(size / binSize);
    bins.resize(numBinsX, std::vector<std::vector<int>>(numBinsY));
}

// Assign particles to bins 
void assign_particles_to_bins(particle_t* parts, int num_parts, double size){
    int numBinsX = ceil(size/binSize);
    int numBinsY = ceil(size/binSize); 

    // clear bins 
    for (auto& row : bins){
        for (auto& bin: row){
            bin.clear(); 
        }
    }

    // Assign particles to bins
    for (int i = 0; i < num_parts; ++i){
        int binX = parts[i].x / binSize; 
        int binY = parts[i].y / binSize; 
        bins[binX][binY].push_back(i); 
    }
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Assign particles to bins 
    assign_particles_to_bins(parts, num_parts, size); 

    // Compute forces with binning 
    for(int i = 0; i < num_parts; ++i){
        parts[i].ax = parts[i].ay = 0; 
        int binX = parts[i].x/binSize; 
        int binY = parts[i].y/binSize; 

        // get neighboring bins 
        auto neighbors = get_neighboring_bins(binX, binY, bins.size(), bins[0].size()); 

        for (auto &nb: neighbors){
            for (int j : bins[nb.first][nb.second]) {
                apply_force(parts[i], parts[j]); 
            }
        }
    }

    // Move particles 
    for (int i = 0; i < num_parts; ++i){
        move(parts[i], size); 
    }

    // Re-assign particles to bins after movement 
    assign_particles_to_bins(parts, num_parts, size); 
}


std::vector<std::pair<int, int>> get_neighboring_bins(int binX, int binY, int maxX, int maxY){
    std::vector<std::pair<int,int>> neighbors; 
    for (int dx = -1; dx <= 1; ++dx){
        for (int dy = -1; dy <= 1; ++dy){
            int newX = binX + dx, newY = binY + dy; 
            if (newX >= 0 && newX < maxX && newY >= 0 && newY < maxY){
                neighbors.emplace_back(newX, newY); 
            }
        }
    }

    return neighbors; 
}
