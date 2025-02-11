# Particle

# Overview  
Simulates particle movement under gravitational forces using Gravitational Force and Equations of Motion to update the position of a particle over time. You may add the number of particles, time step per update, number of simulation steps, and the output intervals to modify your experience. "output.log" is then created to check each interval.
 

# Features  
- Random particle generation, calculation of gravitational forces between each pair of particles.
- Velocity updates based on forces. 
- Runs on Centaurus computing nodes, able to use SLURM batch jobs.
 
 
# Clone the Repository  
To download and use the project on Centaurus ensure you are on a computing node, run:  
```bash
git clone https://github.com/SumJunk/ParticleSim.git
cd ParticleSim
git checkout master
g++ -O3 Particle.cpp -o Particle
./Particle <num_particles> <dt> <num_steps> <output_interval>


