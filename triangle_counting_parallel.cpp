#include "core/graph.h"
#include "core/utils.h"
#include <future>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>

std::mutex mtx;

struct ReturnStruct{
    long triangleCount;
    double time;
};

long countTriangles(uintV *array1, uintE len1, uintV *array2, uintE len2,
                     uintV u, uintV v) {

  uintE i = 0, j = 0; // indexes for array1 and array2
  long count = 0;

  if (u == v)
    return count;

  while ((i < len1) && (j < len2)) {
    if (array1[i] == array2[j]) {
      if ((array1[i] != u) && (array1[i] != v)) {
        count++;
      }
      i++;
      j++;
    } else if (array1[i] < array2[j]) {
      i++;
    } else {
      j++;
    }
  }
  return count;
}

void threadFunction(Graph& graph, uintV beginVertex, uintV endVertex, long& totalTriangleCount, ReturnStruct& result) {
    timer threadTimer;
    threadTimer.start();

    long threadTriangleCount = 0;

    // Process each edge <u,v>
    for (uintV u = beginVertex; u < endVertex; u++) {
        uintE outDegree = graph.vertices_[u].getOutDegree();

        // For each outNeighbor v, find the intersection of inNeighbor(u) and outNeighbor(v)
        for (uintE idx = 0; idx < outDegree; idx++) {
            uintV v = graph.vertices_[u].getOutNeighbor(idx);
            threadTriangleCount += countTriangles(graph.vertices_[u].getInNeighbors(),
                                                  graph.vertices_[u].getInDegree(),
                                                  graph.vertices_[v].getOutNeighbors(),
                                                  graph.vertices_[v].getOutDegree(), u, v);
        }
    }

    {
        mtx.lock();
        totalTriangleCount += threadTriangleCount;
        mtx.unlock();
    }

    result.triangleCount = threadTriangleCount;
    result.time = threadTimer.stop();
}


void triangleCountParallel(Graph &graph, uint numThreads) {
    uintV totalVertices = graph.n_;
    long overallTriangleCount = 0;
    double overallTimeTaken = 0.0;
    timer mainTimer;

    mainTimer.start();

    std::thread workers[numThreads];
    ReturnStruct results[numThreads];

    auto verticesPerThread = totalVertices / numThreads;
    uintV startVertex = 0;
    uintV endVertex = startVertex + verticesPerThread;

    for(uint idx = 0; idx < numThreads; idx++) {
        if(idx == numThreads - 1)
            endVertex = totalVertices;
        workers[idx] = std::thread(threadFunction, std::ref(graph), startVertex, endVertex, std::ref(overallTriangleCount), std::ref(results[idx]));
        startVertex += verticesPerThread;
        endVertex += verticesPerThread;
    }

    for(uint idx = 0; idx < numThreads; idx++) {
        workers[idx].join();
    }
    
    overallTimeTaken = mainTimer.stop();

    std::cout << "thread_id, triangle_count, time_taken\n"; 
    for(uint idx = 0; idx < numThreads; idx++) {
        std::cout << idx << ", " << results[idx].triangleCount << ", " << results[idx].time << "\n";
    }

    std::cout << "Number of triangles : " << overallTriangleCount << "\n";
    std::cout << "Number of unique triangles : " << overallTriangleCount / 3 << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << overallTimeTaken << "\n";
}


int main(int argc, char *argv[]) {
  cxxopts::Options options(
      "triangle_counting_serial",
      "Count the number of triangles using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"nWorkers", "Number of workers",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_workers = cl_options["nWorkers"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();
  std::cout << std::fixed;
  std::cout << "Number of workers : " << n_workers << "\n";

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";

  triangleCountParallel(g, n_workers);

  return 0;
}
