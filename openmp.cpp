#include "common.h"
#include <omp.h>
#include <cmath>
#include <array> 

constexpr int MAX_PARTICLES = 1000000;
// constexpr int MAX_BINS = 10000; // Choose an appropriate value

// Define a struct for linked list node
struct ListNode {
    int index;
    ListNode* next;
    ListNode(int idx) : index(idx), next(nullptr) {}
};

// Global declaration of bins
std::array<std::array<ListNode*, MAX_PARTICLES>, MAX_PARTICLES> bins; // Use fixed-size array instead of vector
int binCountX, binCountY; // Number of bins in each dimension
int binSize; 
// omp_lock_t locks[MAX_BINS][MAX_BINS]; // an array to store the locks
// omp_lock_t** locks; // an array to store the locks 

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
    binSize = std::max(static_cast<int>(std::ceil(cutoff * 0.5)), 1); 
    binCountX = std::ceil(size / binSize);
    binCountY = std::ceil(size / binSize);

    // Initialize locks for bins
    // for (int i = 0; i < binCountX; ++i) {
    //     for (int j = 0; j < binCountY; ++j) {
    //         omp_init_lock(&locks[i][j]);
    //     }
    // }
}

// Assign particles to bins 
void assign_particles_to_bins(particle_t* parts, int num_parts, double size) {
    // Clear bins
    for (int i = 0; i < binCountX; ++i) {
        for (int j = 0; j < binCountY; ++j) {
            bins[i][j] = nullptr; // Initialize all bins to nullptr
        }
    }

    // Assign particles to bins
    // #pragma omp parallel for 
    for (int i = 0; i < num_parts; ++i) {
        int binX = parts[i].x / binSize;
        int binY = parts[i].y / binSize;
        // omp_set_lock(&locks[binX][binY]);
        bins[binX][binY] = new ListNode(i); // Assign a new node to the bin
        // omp_unset_lock(&locks[binX][binY]);
    }
}


void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Re-assign particles to bins to account for movement

    // Compute forces with binning
    #pragma omp for
    for (int i = 0; i < num_parts; ++i) {
        parts[i].ax = parts[i].ay = 0;
        int binX = parts[i].x / binSize;
        int binY = parts[i].y / binSize;

        // Iterate through neighboring bins
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                int newX = binX + dx, newY = binY + dy;
                if (newX >= 0 && newX < binCountX && newY >= 0 && newY < binCountY) {
                    ListNode* node = bins[newX][newY]; // Get the head of the linked list
                    while (node != nullptr) {
                        apply_force(parts[i], parts[node->index]);
                        node = node->next;
                    }
                }
            }
        }
    }

    // Move particles
    #pragma omp for
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }

    // assign_particles_to_bins(parts, num_parts, size);

    // Recalculate bins for particles that moved in the previous time step
    #pragma omp master
    {
        assign_particles_to_bins(parts, num_parts, size);
    }
    // Synchronize threads
    #pragma omp barrier

     // Clean up bins after each simulation step --> nvm makes it slower
    //cleanup_bins();
}