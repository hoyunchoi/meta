#include <iostream>
#include <chrono>
#include <random>
#include <filesystem>

#include "../library-Git/stringFormat.hpp"
#include "meta.hpp"

namespace fs = std::filesystem;

int main(int argc, char *argv[]){
    /*  Input variables
        networkSize = std::stoul(argv[1]);
        meanDegree = std::stod(argv[2]);
        SI_II = std::stod(argv[3]);
        I_S = std::stod(argv[4]);
    */

    //* Weighted meta-population network parameter
    const Size networkSize = std::stoul(argv[1]);               //* Number of nodes of  meta-population network
    const double meanDegree = std::stod(argv[2]);               //* Mean degree of meta-population network
    const Size linkSize = (Size)networkSize * meanDegree / 2;   //* Number of total links of meta-population network
    constexpr Size meanPopulation = 1000;                       //* Average size of node(sub-population)
    const double degreeExponent = 3.0;                          //* Exponent of degree distribution (scale-free)
    const double linkWeightExponent = 0.5;                      //* Link weight = (degree1*degree2)^linkWeightExponent

    //* SIS process parameter
    const double SI_II = std::stod(argv[3]);
    const double I_S = std::stod(argv[4]);
    const double diffusion = 0.5;                               //* Portion of each nodes to diffuse
    SIS_meta::seedSize = 1;                                     //* Number of initial seed(I)
    SIS_meta::maxIteration = 50;
    SIS_meta::deltaT = 0.1;
    const std::vector<double> rates = {SI_II, I_S, diffusion};

    //* Ensemble and write parameters
    constexpr Size ensembleSize = 2;
    const std::string seperater = "\t";
    const std::string secondSeperater = "\n";
    const bool append = true;

    //* Directory and File Name
    const std::string rootPath = "../data/epidemics/meta_SIS/";
    const std::string networkName = "N" + std::to_string(networkSize) + ",M" + to_stringWithPrecision(meanDegree, 0) + ",DE" + to_stringWithPrecision(degreeExponent, 2);
    const std::string rateName = "SIII" + to_stringWithPrecision(SI_II, 2) + ",IS" + to_stringWithPrecision(I_S,2) + ",sD" + to_stringWithPrecision(diffusion, 2);
    const std::string directory = rootPath + networkName + "/" + rateName + "/";

    if (!fs::exists(directory)){
        fs::create_directories(directory);
    }

    //*------------------------------------------------------------------------------------
    auto start = std::chrono::system_clock::now();

    for (unsigned ensemble=0; ensemble<ensembleSize; ++ensemble){
        //* Generate CL network and write
        pcg32 randomEngine(ensemble);
        const WNetwork network = WCL::generate(networkSize, meanPopulation, meanDegree, degreeExponent, linkWeightExponent, randomEngine);
        network.printAdjacency(directory + "network_bool.txt", seperater, secondSeperater, append);
        network.printWAdjacency(directory + "network_weight.txt", seperater, secondSeperater, append);

        //* Do SIS_meta and write
        SIS_meta::syncRun(network, rates, directory + "nodewise.txt", directory + "total.txt");
    }

    std::chrono::duration<double> sec = std::chrono::system_clock::now()-start;
    printf("%0.10f second\n", sec.count());

    return 0;
}