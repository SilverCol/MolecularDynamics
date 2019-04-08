#include <gflags/gflags.h>
#include "experiments.hpp"

DEFINE_double(relax, 1.0, "Inverse relaxation time.");
DEFINE_double(tl, 1.0, "Left thermostat.");
DEFINE_double(tr, 2.0, "Right thermostat.");
DEFINE_double(lambda, 1.0, "Anharmonic constant.");

DEFINE_int32(reads, 1000, "Number of reads for output file.");
DEFINE_int32(samples, 1000, "Number of tau steps between reads (maxwell).");
DEFINE_int32(steps, 1, "Number of steps between reads (or inside sample for maxwell).");
DEFINE_double(step, .03, "Time step (a.k.a. tau for maxwell).");

DEFINE_string(file, "../data/a.txt", "Path for the output file.");

DEFINE_int32(mode, 0, "Operation mode: 0-TProfile, 1-Flux, 2-MTProfile");

int main(int argc, char* argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    gflags::ParseCommandLineFlags(&argc, &argv, true);

    double params[4] = {FLAGS_relax, FLAGS_tl, FLAGS_tr, FLAGS_lambda};

    std::normal_distribution<double> tl(0.0, sqrt(FLAGS_tl));
    std::normal_distribution<double> tr(0.0, sqrt(FLAGS_tr));

    // TODO move those inside function
    double y[DIM];
    stateInit(y);

    std::vector<double> output;
    switch(FLAGS_mode)
    {
        case 0: // T profile mode
            std::cout << "Calculating T profile" << std::endl;
            FLAGS_file.append(std::to_string(N));
            makeTProfile(FLAGS_step, FLAGS_steps, params, y, output);
           break;
        case 1: // flux mode
            std::cout << "Calculating T profile" << std::endl;
            FLAGS_file.append(std::to_string(N));
            makeFlux(FLAGS_step, FLAGS_steps, params, y, output);
           break;
        case 2: // maxwell T profile mode
            std::cout << "Calculating T profile (maxwell)" << std::endl;
            FLAGS_file.append(std::to_string(N));
            maxwelTProfile(FLAGS_step, FLAGS_reads, FLAGS_samples, FLAGS_steps, FLAGS_lambda, tl, tr, output);
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