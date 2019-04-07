#include <gflags/gflags.h>
#include "experiments.hpp"

DEFINE_int32(N, 50, "Number of particles.");

DEFINE_int32(M, 100, "Number of propagation steps.");
DEFINE_int32(steps, 1, "Number invisible steps inbetween.");
DEFINE_double(z, .03, "Propagation step coefficient.");
DEFINE_int32(scheme, 2, "Number of split-step scheme as S{scheme} (eg. 2 for S2).");

DEFINE_string(file, "../data/a.txt", "Path for the output file.");

DEFINE_int32(mode, 0, "Operation mode: 0-phaseSum, 1-localSpin, 2-spinFlux");

DEFINE_int32(j, 0, "Spin index, for local spin correlation");

int main(int argc, char* argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    gflags::ParseCommandLineFlags(&argc, &argv, true);

    std::vector<double> output;
    switch(FLAGS_mode)
    {
        case 0: // mode
            std::cout << "Calculating ..." << std::endl;
           break;
        default:
            std::cerr << "Invalid mode." << std::endl;
            return -1;
    }
    FLAGS_file.append(".bin");
    std::cout << "Writting to file: " << FLAGS_file << std::endl;
    writeBinary(output, FLAGS_file);

    auto finish = std::chrono::high_resolution_clock::now();
    std::cout   << "Finished in "
                << std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count()
                << "ms" << std::endl;
    return 0;
}