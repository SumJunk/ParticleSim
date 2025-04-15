#include <iostream>
#include <omp.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <string>

struct Particle {
    double mass;
    double posX, posY, posZ;
    double velX, velY, velZ;
    double forceX, forceY, forceZ;
};

class Simulation {
public:
    std::vector<Particle> particles;
    double G = 6.67430e-11;

    void initializeRandomParticles(int numParticles) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-1.0e7, 1.0e7);

        particles.clear();
        for (int i = 0; i < numParticles; ++i) {
            Particle p;
            p.mass = dis(gen);
            p.posX = dis(gen);
            p.posY = dis(gen);
            p.posZ = dis(gen);
            p.velX = dis(gen);
            p.velY = dis(gen);
            p.velZ = dis(gen);
            p.forceX = p.forceY = p.forceZ = 0.0;
            particles.push_back(p);
        }
    }

    void loadFromFile(const std::string &filename) {
        std::ifstream file(filename);
        int numParticles;
        file >> numParticles;
        particles.clear();
        for (int i = 0; i < numParticles; ++i) {
            Particle p;
            file >> p.mass >> p.posX >> p.posY >> p.posZ >> p.velX >> p.velY >> p.velZ >> p.forceX >> p.forceY >> p.forceZ;
            particles.push_back(p);
        }
    }

    void calculateGravitationalForces() {
    
        #pragma omp parallel for
        for (int i = 0; i < particles.size(); ++i) {
            particles[i].forceX = particles[i].forceY = particles[i].forceZ = 0.0;
        }
    
        #pragma omp parallel
        {
            std::vector<double> localForceX(particles.size(), 0.0);
            std::vector<double> localForceY(particles.size(), 0.0);
            std::vector<double> localForceZ(particles.size(), 0.0);
    
            #pragma omp for schedule(dynamic)
            for (int i = 0; i < particles.size(); ++i) {
                for (int j = i + 1; j < particles.size(); ++j) {
                    Particle &p1 = particles[i];
                    Particle &p2 = particles[j];
    
                    double dx = p2.posX - p1.posX;
                    double dy = p2.posY - p1.posY;
                    double dz = p2.posZ - p1.posZ;
                    double distSqr = dx * dx + dy * dy + dz * dz;
                    double dist = std::sqrt(distSqr);
                    double force = (G * p1.mass * p2.mass) / distSqr;
    
                    double fx = force * dx / dist;
                    double fy = force * dy / dist;
                    double fz = force * dz / dist;
    
                    localForceX[i] += fx;
                    localForceY[i] += fy;
                    localForceZ[i] += fz;
    
                    localForceX[j] -= fx;
                    localForceY[j] -= fy;
                    localForceZ[j] -= fz;
                }
            }
            #pragma omp critical
            {
                for (int i = 0; i < particles.size(); ++i) {
                    particles[i].forceX += localForceX[i];
                    particles[i].forceY += localForceY[i];
                    particles[i].forceZ += localForceZ[i];
                }
            }
        }
    }

    void applyForces(double dt) {
        #pragma omp parallel for
        for (int i = 0; i < particles.size(); ++i) {
            particles[i].velX += (particles[i].forceX / particles[i].mass) * dt;
            particles[i].velY += (particles[i].forceY / particles[i].mass) * dt;
            particles[i].velZ += (particles[i].forceZ / particles[i].mass) * dt;
        }
    }
    
    void moveParticles(double dt) {
        #pragma omp parallel for
        for (int i = 0; i < particles.size(); ++i) {
            particles[i].posX += particles[i].velX * dt;
            particles[i].posY += particles[i].velY * dt;
            particles[i].posZ += particles[i].velZ * dt;
        }
    }

    void logSimulationState(std::ofstream &outputFile) {
        outputFile << particles.size() << "\t";
        for (const auto &p : particles) {
            outputFile << p.mass << "\t" << p.posX << "\t" << p.posY << "\t" << p.posZ << "\t"
                       << p.velX << "\t" << p.velY << "\t" << p.velZ << "\t"
                       << p.forceX << "\t" << p.forceY << "\t" << p.forceZ << "\t";
        }
        outputFile << std::endl;
    }
};

int main(int argc, char *argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <num_particles> <dt> <num_steps> <output_interval>" << std::endl;
        return 1;
    }

    int numParticles = std::stoi(argv[1]);
    double dt = std::stod(argv[2]);
    int numSteps = std::stoi(argv[3]);
    int outputInterval = std::stoi(argv[4]);

    Simulation sim;
    sim.initializeRandomParticles(numParticles);
    std::ofstream outputFile("output.log");

    for (int step = 0; step < numSteps; ++step) {
        sim.calculateGravitationalForces();
        sim.applyForces(dt);
        sim.moveParticles(dt);

        if (step % outputInterval == 0) {
            sim.logSimulationState(outputFile);
        }
    }

    outputFile.close();
    return 0;
}
