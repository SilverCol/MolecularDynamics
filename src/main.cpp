#include <gflags/gflags.h>
#include "experiments.hpp"

DEFINE_double(relax, 1.0, "Inverse relaxation time.");
DEFINE_double(tl, 1.0, "Left thermostat.");
DEFINE_double(tr, 2.0, "Right thermostat.");
DEFINE_double(lambda, 1.0, "Anharmonic constant.");

DEFINE_int32(reads, 1000, "Number of reads for output file.");
DEFINE_int32(samples, 1000, "Number of tau steps (samples) between reads (maxwell).");
DEFINE_int32(steps, 1, "Number of steps inside a sample (maxwell).");
DEFINE_double(step, .03, "Time step (a.k.a. tau for maxwell).");
DEFINE_double(cutoff, 10000., "Initial cutoff.");

DEFINE_string(file, "../data/a.txt", "Path for the output file.");

DEFINE_int32(mode, 0, "Operation mode: 0-TProfile, 1-Flux, 2-MTProfile, 3-MFlux");

int main(int argc, char* argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    gflags::ParseCommandLineFlags(&argc, &argv, true);

    double params[4] = {FLAGS_relax, FLAGS_tl, FLAGS_tr, FLAGS_lambda};

    std::normal_distribution<double> tl(0.0, sqrt(FLAGS_tl));
    std::normal_distribution<double> tr(0.0, sqrt(FLAGS_tr));

    std::vector<double> x;
    std::vector<double> output;
    switch(FLAGS_mode)
    {
        case 0: // T profile mode
            std::cout << "Calculating T profile" << std::endl;
            FLAGS_file.append(std::to_string(N));
            x.resize(2*N + 2, 0);
            stateInit(x);
            makeTProfile(FLAGS_step * FLAGS_samples, FLAGS_reads, params, output, x, FLAGS_cutoff);
           break;
        case 1: // flux mode
            std::cout << "Calculating flux" << std::endl;
            FLAGS_file.append(std::to_string(N));
            x.resize(2*N + 2, 0);
            stateInit(x);
            makeFlux(FLAGS_step * FLAGS_samples, FLAGS_reads, params, output, x);
           break;
        case 2: // maxwell T profile mode
            std::cout << "Calculating T profile (maxwell)" << std::endl;
            FLAGS_file.append(std::to_string(N));
            x.resize(2*N, 0);
            maxwellInit(x);
            maxwelTProfile(FLAGS_step, FLAGS_reads, FLAGS_samples, FLAGS_steps, FLAGS_lambda, tl, tr, output, x,
                    FLAGS_cutoff);
            break;
        case 3: // maxwell flux mode
            std::cout << "Calculating flux (maxwell)" << std::endl;
            FLAGS_file.append(std::to_string(N));
            x.resize(2*N, 0);
            maxwellInit(x);
            maxwelFlux(FLAGS_step, FLAGS_reads, FLAGS_samples, FLAGS_steps, FLAGS_lambda, tl, tr, output, x);
            break;
        default:
            std::cerr << "Invalid mode." << std::endl;
            return -1;
    }
    if ((int)FLAGS_lambda != 0) FLAGS_file.append("a");
    FLAGS_file.append(".bin");
    std::cout << "Writting to file: " << FLAGS_file << std::endl;
    writeBinary(output, FLAGS_file);

    auto finish = std::chrono::high_resolution_clock::now();
    std::cout   << "Finished in "
                << std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count()
                << "ms" << std::endl;
    return 0;
}