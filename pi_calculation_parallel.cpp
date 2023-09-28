#include "core/utils.h"
#include <iostream>
#include <vector>
#include <thread>
#include <atomic>
#include <iomanip>
#include <stdlib.h>
#include <numeric>

#define sqr(x) ((x) * (x))
#define DEFAULT_NUMBER_OF_POINTS "1000"
#define DEFAULT_NUMBER_OF_WORKERS "4"

uint c_const = (uint)RAND_MAX + (uint)1;
inline double get_random_coordinate(uint *random_seed) {
  return ((double)rand_r(random_seed)) / c_const;
}

uint get_points_in_circle(uint n, uint random_seed) {
  uint circle_count = 0;
  double x_coord, y_coord;
  for (uint i = 0; i < n; i++) {
    x_coord = (2.0 * get_random_coordinate(&random_seed)) - 1.0;
    y_coord = (2.0 * get_random_coordinate(&random_seed)) - 1.0;
    if ((sqr(x_coord) + sqr(y_coord)) <= 1.0)
      circle_count++;
  }
  return circle_count;
}

void threadPiCalculation(uint n, uint tid, std::atomic<uint> &global_circle_count, std::vector<uint> &local_circle_count, std::vector<double> &time_per_thread) {
  timer thread_timer;
  thread_timer.start();

  uint circle_count = get_points_in_circle(n, tid);
  
  local_circle_count[tid]  = circle_count;
  global_circle_count += circle_count;

  time_per_thread[tid] = thread_timer.stop();
 

}

int main(int argc, char *argv[]) {
     cxxopts::Options options("pi_calculation",
                           "Calculate pi using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"nPoints", "Number of points",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_POINTS)},
          {"nWorkers", "Number of workers",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
      });

  auto cl_options = options.parse(argc, argv);
  
    uint n_points = cl_options["nPoints"].as<uint>();
    uint n_workers = cl_options["nWorkers"].as<uint>();
    std::cout << std::fixed;
    std::cout << "Number of points : " << n_points << "\n";
    std::cout << "Number of workers : " << n_workers << "\n";

    std::vector<std::thread> threads(n_workers);
    std::vector<uint> local_circle_count(n_workers, 0);
    std::vector<double> time_per_thread(n_workers);

    std::atomic<uint> global_circle_count(0);

    uint base_points_per_thread = n_points / n_workers;
    uint remainder = n_points % n_workers;



    
    for (uint i = 0; i < n_workers; i++) {
       uint current_points = base_points_per_thread + (i < remainder ? 1 : 0);
        threads[i] = std::thread(threadPiCalculation, current_points, i, std::ref(global_circle_count), std::ref(local_circle_count), std::ref(time_per_thread));
    }

    for (uint i = 0; i < n_workers; i++) {
        threads[i].join();
    }

    uint total_circle_count = std::accumulate(local_circle_count.begin(), local_circle_count.end(), 0u);
    double pi_value = 4.0 * (double)global_circle_count.load() / (double)n_points;

    std::cout << "thread_id, points_generated, circle_points, time_taken\n";
    for (uint i = 0; i < n_workers; i++) {
        std::cout << i << ", " << (base_points_per_thread + (i < remainder ? 1 : 0)) << ", " << local_circle_count[i] << ", " << time_per_thread[i] << "\n";
    }

    std::cout << "Total points generated : " << n_points << "\n";
    std::cout << "Total points in circle : " << global_circle_count.load() << "\n";
    std::cout << "Result : " << std::setprecision(VAL_PRECISION) << pi_value << "\n";

    double max_time = *std::max_element(time_per_thread.begin(), time_per_thread.end()); // Get max time
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION) << max_time << "\n";

    return 0;
}
