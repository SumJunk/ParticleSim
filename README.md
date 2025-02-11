# MergeSort

# Overview  
Implements a mergesort algorithm by sorting an array of a given size in order while producing the time it takes to sort.
 

# Features  
- Show execution times for 10^1 to 10^9.  
- Runs on Centaurus computing nodes, able to use SLURM batch jobs
- Chart to also show execution times.
 
 
# Clone the Repository  
To download and use the project on Centaurus ensure you are on a computing node, run:  
```bash
git clone https://github.com/SumJunk/MergeSort.git
cd MergeSort
git checkout master
g++ -O3 mergesort.cpp -o mergesort
./mergesort (value)


