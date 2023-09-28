#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <thread>
#include <mutex>

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
typedef int64_t PageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
typedef float PageRankType;
#endif

std::mutex mutexLock;
const uintV vertices_per_mutex = 1000;

void pageRankParallelWorker(uintV start, uintV end, Graph &g, int max_iters, PageRankType* pr_curr, PageRankType* pr_next, std::mutex* locks, CustomBarrier& my_barrier, double &time_taken) {
    timer t1;
    t1.start();

   
    
   for (int iter = 0; iter < max_iters; iter++) {
        for (uintV u = start; u < end; u++) {
            uintE out_degree = g.vertices_[u].getOutDegree();
            for (uintE i = 0; i < out_degree; i++) {
                uintV v = g.vertices_[u].getOutNeighbor(i);
                uintV mutex_index = v / vertices_per_mutex;  // Find the mutex associated with vertex v
                locks[mutex_index].lock();
                pr_next[v] += (pr_curr[u] / out_degree);
                locks[mutex_index].unlock();
            }
        }
        
        my_barrier.wait();  // Barrier to ensure all threads finish the first loop
        
        for (uintV v = start; v < end; v++) {
            pr_next[v] = PAGE_RANK(pr_next[v]);
            pr_curr[v] = pr_next[v];
            pr_next[v] = 0.0;
        }
        
        my_barrier.wait();  // Barrier to ensure all threads finish the second loop before next iteration
    }
    time_taken = t1.stop();
}

void pageRankParallel(Graph &g, int max_iters, uint n_workers) {
    uintV n = g.n_;
    PageRankType *pr_curr = new PageRankType[n];
    PageRankType *pr_next = new PageRankType[n];
    
    
    uintV n_mutexes = (n + vertices_per_mutex - 1) / vertices_per_mutex;
    std::mutex *locks = new std::mutex[n_mutexes];

    CustomBarrier my_barrier(n_workers);  // Barrier object

    for (uintV i = 0; i < n; i++) {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    std::vector<std::thread> workers;
    std::vector<double> times(n_workers, 0.0);
    uintV vertices_per_thread = n / n_workers;

    for (uint t = 0; t < n_workers; t++) {
        uintV start = t * vertices_per_thread;
        uintV end = (t == n_workers - 1) ? n : (t + 1) * vertices_per_thread;
        workers.push_back(std::thread(pageRankParallelWorker, start, end, std::ref(g), max_iters, pr_curr, pr_next, locks, std::ref(my_barrier), std::ref(times[t])));
    }

    for (std::thread &worker : workers) {
        worker.join();
    }

    std::cout << "thread_id, time_taken\n";
    for (uint t = 0; t < n_workers; t++) {
        std::cout << t << ", " << times[t] << "\n";
    }

    PageRankType sum_of_page_ranks = 0;
    for (uintV u = 0; u < n; u++) {
        sum_of_page_ranks += pr_curr[u];
    }
    std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
    std::cout << "Time taken (in seconds) : " << *max_element(times.begin(), times.end()) << "\n";

    delete[] pr_curr;
    delete[] pr_next;
    delete[] locks;
}

int main(int argc, char *argv[]) {
      cxxopts::Options options(
      "page_rank_push",
      "Calculate page_rank using serial and parallel execution");
  options.add_options(
      "",
      {
          {"nWorkers", "Number of workers",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
          {"nIterations", "Maximum number of iterations",
           cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_workers = cl_options["nWorkers"].as<uint>();
  uint max_iterations = cl_options["nIterations"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();

#ifdef USE_INT
  std::cout << "Using INT\n";
#else
  std::cout << "Using FLOAT\n";
#endif
  std::cout << std::fixed;
  std::cout << "Number of workers : " << n_workers << "\n";

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
    std::cout << "Created graph\n";
    pageRankParallel(g, max_iterations, n_workers);

    return 0;
}
